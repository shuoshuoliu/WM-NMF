function [fout,xout,u]=statkscompute(ftype,xi,xispecified,m,u,L,U,weight,cutoff,...
                                     kernelname,ty,yData,foldpoint,maxp,varargin)
%STATKSCOMPUTE Perform computations for kernel smoothing density
%estimation.

%   Copyright 2008-2017 The MathWorks, Inc.

if nargin < 15
    d = 1;
    isXChunk = false;
    bdycorr = 'log';
else
    d = varargin{1};
    isXChunk = varargin{2};
    bdycorr = varargin{3};
    if isempty(bdycorr)
        bdycorr = 'log';
    end
end
[kernel,iscdf,kernelcutoff,kernelname,ftype] = statkskernelinfo(ftype,kernelname,d);

if isempty(cutoff)
    cutoff = kernelcutoff;
end

% Inverse cdf is special, so deal with it here
if isequal(ftype,'icdf') % d==1
    % Put the foldpoint and yData on the transformed scale
    yData = yData(:);
    yData = transform(yData,L,U,1,bdycorr);
    foldpoint = transform(foldpoint,L,U,1,bdycorr);

    % Compute on that scale
    [fout,xout]=compute_icdf(xi,xispecified,m,u,L,U,weight,cutoff,...
                             kernelname,ty,yData,foldpoint,maxp,1,bdycorr);

    % Transform back - f(p) is the icdf at p
    if isequal(bdycorr, 'log')
        fout = untransform(fout,L,U);
    end
else
    [fout,xout,u] = compute_pdf_cdf(xi,xispecified,m,L,U,weight,kernel,...
                                    cutoff,iscdf,u,ty,foldpoint,d,isXChunk,bdycorr);

    % If another function based on the cdf, compute it now
    if isequal(ftype,'survivor')
        fout = 1-fout;
    elseif isequal(ftype,'cumhazard')
        fout = 1-fout;
        t = (fout>0);
        fout(~t) = NaN;
        fout(t) = -log(fout(t));
    end
end

% -----------------------------
function [fout,xout,u]=compute_pdf_cdf(xi,xispecified,m,L,U,weight,...
                          kernel,cutoff,iscdf,u,ty,foldpoint,d,isXChunk,bdycorr)

foldwidth = min(cutoff,3);
issubdist = isfinite(foldpoint);
if ~xispecified
    xi = compute_default_xi(ty,foldwidth,issubdist,m,u,U,L,d,bdycorr);
end

% Compute transformed values of evaluation points that are in bounds
if d == 1
    fout = zeros(size(xi)); % f has the same orientation as xi
    xisize = numel(xi);
else
    xisize = size(xi,1);
    fout = zeros(xisize,1);
end
if iscdf && all(isfinite(U))
    fout(all(bsxfun(@ge,xi,U),2)) = sum(weight);
end
xout = xi;
if d == 1
    xi = xi(:);
end

if all(L==-Inf) && all(U==Inf)   % unbounded support
    inbounds = true(xisize,1);
elseif all(L==0) && all(U==Inf)  % positive support
    inbounds = all(xi>0,2);
    xi = xi(inbounds,:);
else % finite support [L, U]
    inbounds = all(bsxfun(@gt,xi,L),2) & all(bsxfun(@lt,xi,U),2);
    xi = xi(inbounds,:);
end
txi = transform(xi,L,U,d,bdycorr);
if d == 1
    foldpoint = transform(foldpoint,L,U,1,bdycorr);
end

% If the density is censored at the end, add new points so that we can fold
% them back across the censoring point as a crude adjustment for bias.
if issubdist
    needfold = (txi >= foldpoint - foldwidth*u);
    txifold = (2*foldpoint) - txi(needfold);
    nfold = sum(needfold);
else
    nfold = 0;
end

if isempty(xi)
    f = xi(:);
else
    % Compute kernel estimate at the requested points
    f = dokernel(iscdf,txi,ty,u,weight,kernel,cutoff,d,L,U,xi,bdycorr);
    
    % If we need extra points for folding, do that now
    if nfold>0
        % Compute the kernel estimate at these extra points
        ffold = dokernel(iscdf,txifold,ty,u,weight,kernel,cutoff,d,L,U,xi(needfold),bdycorr);
        if iscdf
            % Need to use upper tail for cdf at folded points
            ffold = sum(weight) - ffold;
        end
        
        % Fold back over the censoring point
        f(needfold) = f(needfold) + ffold;
        
        if iscdf
            % For cdf, extend last value horizontally
            maxf = max(f(txi<=foldpoint));
            f(txi>foldpoint) = maxf;
        else
            % For density, define a crisp upper limit with vertical line
            f(txi>foldpoint) = 0;
            if ~xispecified
                xi(end+1) = xi(end);
                f(end+1) = 0;
                inbounds(end+1) = true;
            end
        end
    end
end

if iscdf
    if ~isXChunk
        % Guard against roundoff.  Lower boundary of 0 should be no problem.
        f = min(1,f);
    end
else
    f = f(:) ./ prod(u);
end
fout(inbounds) = f;
if d==1
    xout(inbounds) = xi;
end

% -----------------------------
function xi = compute_default_xi(ty,foldwidth,issubdist,m,u,U,L,d,bdycorr)
% Get XI values at which to evaluate the density

% Compute untransformed values of lower and upper evaluation points
ximin = min(ty,[],1) - foldwidth*u;
if issubdist
    ximax = max(ty,[],1);
else
    ximax = max(ty,[],1) + foldwidth*u;
end

if isequal(bdycorr,'log')
    for i =1:d
        ximin(i) = untransform(ximin(i),L(i),U(i));
        ximax(i) = untransform(ximax(i),L(i),U(i));
    end
else
    for i =1:d
        ximin(i) = max(ximin(i),L(i));
        ximax(i) = min(ximax(i),U(i));
    end
end

if d == 1
    xi = linspace(ximin, ximax, m);
else
    x1 = linspace(ximin(1), ximax(1), m);
    x2 = linspace(ximin(2), ximax(2), m);
    [x1,x2] = meshgrid(x1,x2);
    x1 = x1(:);
    x2 = x2(:);
    xi = [x1 x2];
end

% -----------------------------
function f = dokernel(iscdf,txi,ty,u,weight,kernel,cutoff,d,L,U,xi,bdycorr)
% Now compute density estimate at selected points
blocksize = 3e4;
if d == 1
    m = length(txi);
    n = length(ty);
else
    m = size(txi,1);
    n = size(ty,1);
end

needUntransform = (~(isinf(L)&L<0)|~isinf(U)) & ~iscdf & isequal(bdycorr,'log');
reflectionPDF = isequal(bdycorr,'reflection') && ~iscdf;
reflectionCDF = isequal(bdycorr,'reflection') && iscdf;

if n*m<=blocksize && ~iscdf
    % For small problems, compute kernel density estimate in one operation
    ftemp = ones(n,m);   
    for i = 1:d
        z = (repmat(txi(:,i)',n,1)-repmat(ty(:,i),1,m))/u(i);
        if reflectionPDF
            zleft = (repmat(txi(:,i)',n,1)+ repmat(ty(:,i),1,m) - 2*L(i))/u(i);
            zright = (repmat(txi(:,i)',n,1) + repmat(ty(:,i),1,m) - 2*U(i))/u(i);
            f =  feval(kernel, z)+ feval(kernel, zleft) + feval(kernel, zright);
        else
            f =  feval(kernel, z);
        end
        % Apply reverse transformation and create return value of proper size
        if needUntransform(i)
            f = untransform_f(f,L(i),U(i),xi(:,i));
        end
        ftemp = ftemp.*f;
    end
    f = weight * ftemp;
elseif d == 1
    % For large problems, try more selective looping

    % First sort y and carry along weights
    [ty,idx] = sort(ty);
    weight = weight(idx);

    % Loop over evaluation points
    f = zeros(1,m);

    if isinf(cutoff)
        if reflectionCDF
            fc = compute_CDFreduction(L,U,u,Inf,n,ty,weight,kernel);
        end
        for k=1:m
            % Sum contributions from all
            z = (txi(k)-ty)/u;
            if reflectionPDF
                zleft = (txi(k)+ty-2*L)/u;
                zright = (txi(k)+ty-2*U)/u;
                f(k) = weight * (feval(kernel,z)+feval(kernel,zleft)+feval(kernel,zright));
            elseif reflectionCDF
                zleft = (txi(k)+ty-2*L)/u;
                zright = (txi(k)+ty-2*U)/u;
                fk = weight * (feval(kernel,z)+feval(kernel,zleft)+feval(kernel,zright));
                f(k) = fk - fc;
            else
                f(k) = weight * feval(kernel,z);
            end
        end
            
    else
        % Sort evaluation points and remember their indices
        [stxi,idx] = sort(txi);

        jstart = 1;       % lowest nearby point
        jend = 1;         % highest nearby point
        halfwidth = cutoff*u;
        
        % Calculate reduction for reflectionCDF
        if reflectionCDF
            fc = compute_CDFreduction(L,U,u,halfwidth,n,ty,weight,kernel);
        end
        
        for k=1:m
            % Find nearby data points for current evaluation point
            lo = stxi(k) - halfwidth;
            while(ty(jstart)<lo && jstart<n)
                jstart = jstart+1;
            end
            hi = stxi(k) + halfwidth;
            jend = max(jend,jstart);
            while(ty(jend)<=hi && jend<n)
                jend = jend+1;
            end
            nearby = jstart:jend;

            % Sum contributions from these points
            z = (stxi(k)-ty(nearby))/u;
            if reflectionPDF
                zleft = (stxi(k)+ty(nearby)-2*L)/u;
                zright = (stxi(k)+ty(nearby)-2*U)/u;
                fk = weight(nearby) * (feval(kernel,z)+feval(kernel,zleft)+feval(kernel,zright));
            elseif reflectionCDF
                zleft = (stxi(k)+ty(nearby)-2*L)/u;
                zright = (stxi(k)+ty(nearby)-2*U)/u;
                fk = weight(nearby) * feval(kernel,z);
                fk = fk + sum(weight(1:jstart-1));
                if jstart == 1
                    fk = fk + weight(nearby) * feval(kernel,zleft);
                    fk = fk + sum(weight(jend+1:end));
                else
                    fk = fk + sum(weight);
                end
                if jend == n
                    fk = fk + weight(nearby) * feval(kernel,zright);
                end
                fk = fk - fc;
            elseif ~iscdf
                fk = weight(nearby) * feval(kernel,z);
            elseif iscdf
                fk = weight(nearby) * feval(kernel,z);
                fk = fk + sum(weight(1:jstart-1));
            end
            f(k) = fk;
        end

        % Restore original x order
        f(idx) = f;
    end
    
    if needUntransform
        f = untransform_f(f,L,U,xi);
    end

else % d > 1
    % Calculate reduction for reflectionCDF
    if reflectionCDF
        cutoff = Inf;
        fc = zeros(n,d);
        for j = 1:d
            fc(:,j) = compute_CDFreduction_mv(L(j),U(j),u(j),ty(:,j),kernel);
        end
    end
        
    if isinf(cutoff)
        f = zeros(1,m);
        for i = 1:m
            ftemp = ones(n,1);
            for j = 1:d
                z = (txi(i,j) - ty(:,j))./u(j);
                if reflectionPDF
                    zleft = (txi(i,j) + ty(:,j)-2*L(j))./u(j);
                    zright = (txi(i,j) + ty(:,j)-2*U(j))./u(j);
                    fk = feval(kernel,z)+feval(kernel,zleft)+feval(kernel,zright); 
                elseif reflectionCDF
                    zleft = (txi(i,j) + ty(:,j)-2*L(j))./u(j);
                    zright = (txi(i,j) + ty(:,j)-2*U(j))./u(j);
                    fk = feval(kernel,z)+feval(kernel,zleft)+feval(kernel,zright); 
                    fk = fk - fc(:,j);
                else
                    fk = feval(kernel,z); 
                end
                if needUntransform(j)
                    fk = untransform_f(fk,L(j),U(j),xi(i,j));
                end
                ftemp = ftemp.*fk;
            end
            f(i) = weight * ftemp;
        end
    else
        halfwidth = cutoff*u;
        index = (1:n)';
        f = zeros(1,m);
        for i = 1:m
            Idx = true(n,1);
            cdfIdx = true(n,1);
            cdfIdx_allBelow = true(n,1);
            for j = 1:d
                dist = txi(i,j) - ty(:,j);
                currentIdx = abs(dist) <= halfwidth(j);
                Idx = currentIdx & Idx; % pdf boundary
                if iscdf
                    currentCdfIdx = dist >= -halfwidth(j);
                    cdfIdx = currentCdfIdx & cdfIdx; % cdf boundary1, equal or below the query point in all dimension
                    currentCdfIdx_below = dist - halfwidth(j) > 0;                   
                    cdfIdx_allBelow = currentCdfIdx_below & cdfIdx_allBelow; % cdf boundary2, below the pdf lower boundary in all dimension
                end
            end
            if ~iscdf
                nearby = index(Idx);
            else
                nearby = index((Idx|cdfIdx)&(~cdfIdx_allBelow));
            end
            if ~isempty(nearby)
                ftemp = ones(length(nearby),1);
                for k =1:d
                    z = (txi(i,k) - ty(nearby,k))./u(k);
                    if reflectionPDF
                        zleft = (txi(i,k) + ty(nearby,k)-2*L(k))./u(k);
                        zright = (txi(i,k) + ty(nearby,k)-2*U(k))./u(k);
                        fk = feval(kernel,z)+feval(kernel,zleft)+feval(kernel,zright);
                    else
                        fk = feval(kernel,z);
                    end
                    if needUntransform(k)
                        fk = untransform_f(fk,L(k),U(k),xi(i,k));
                    end
                    ftemp = ftemp.*fk;
                end
                f(i) = weight(nearby) * ftemp;
            end
            if iscdf && any(cdfIdx_allBelow)
                f(i) = f(i) + sum(weight(cdfIdx_allBelow));
            end
        end
    end
end

% -----------------------------
function [x1,p] = compute_icdf(xi,xispecified,m,u,Lin,Uin,...
    weight,cutoff,kernelname,ty,yData,foldpoint,maxp,d,bdycorr)
if isequal(bdycorr, 'log')
    L = -Inf; % log correction has no bounds
    U = Inf;
else
    L = Lin;
    U = Uin;
end

if xispecified
    p = xi;
else
    p = (1:m)/(m+1);
end

[Fi,xi,cutoff,u] = compute_initial_icdf(m,u,L,U,weight,cutoff,...
    kernelname,ty,yData,foldpoint,d,bdycorr);

[kernel_c,iscdf_c] = statkskernelinfo('cdf',kernelname,d);
[kernel_p,iscdf_p] = statkskernelinfo('pdf',kernelname,d);


% Get starting values for ICDF(p) by inverse linear interpolation of
% the gridded CDF, plus some clean-up
x1 = interp1(Fi,xi,p);               % interpolate for p in a good range
x1(isnan(x1) & p<min(Fi)) = min(xi); % use lowest x if p>0 too low
x1(isnan(x1) & p>max(Fi)) = max(xi); % use highest x if p<1 too high
x1(p<0 | p>maxp) = NaN;              % out of range
x1(p==0) = L;                        % use lower bound if p==0
x1(p==maxp) = U;                     % and upper bound if p==1 or other max

% Now refine the ICDF using Newton's method for cases with 0<p<1
notdone = find(p>0 & p<maxp);
maxiter = 100;
min_dF0 = sqrt(eps(class(p)));
for iter = 1:maxiter
    if isempty(notdone), break; end
    x0 = x1(notdone);

    % Compute cdf and derivative (pdf) at this value
    F0 = compute_pdf_cdf(x0,true,m,L,U,weight,kernel_c,cutoff,...
                         iscdf_c,u,ty,foldpoint,d,false,bdycorr);
    dF0 = compute_pdf_cdf(x0,true,m,L,U,weight,kernel_p,cutoff,...
                          iscdf_p,u,ty,foldpoint,d,false,bdycorr);
                      
    % dF0 is always >= 0. Prevent dF0 from becoming too small.
    dF0 = max(dF0,min_dF0);

    % Perform a Newton's step
    dp = p(notdone) - F0;
    dx = dp ./ dF0;
    x1(notdone) = x0 + dx;

    % Continue if the x and function (probability) change are large
    notdone = notdone(  abs(dx) > 1e-6*abs(x0) ...
                      & abs(dp) > 1e-8 ...
                      & x0 < foldpoint);
end
if ~isempty(notdone)
    warning(message('stats:ksdensity:NoConvergence', num2str( p( notdone( 1 ) ) )));
end

% -----------------------------
function [Fi,xi,cutoff,u] = compute_initial_icdf(m,u,L,U,weight,cutoff,...
    kernelname,ty,yData,foldpoint,d,bdycorr)
% To get starting x values for the ICDF evaluated at p, first create a
% grid xi of values spanning the data on which to evaluate the CDF
sy = sort(yData);
xi = linspace(sy(1), sy(end), 100);

% Estimate the CDF on the grid
[kernel_c,iscdf_c,kernelcutoff] = statkskernelinfo('cdf',kernelname,d);

[Fi,xi,u] = compute_pdf_cdf(xi,true,m,L,U,weight,kernel_c,cutoff,...
                            iscdf_c,u,ty,foldpoint,d,false,bdycorr);

if isequal(kernelname,'normal')
    % Truncation for the normal kernel creates small jumps in the CDF.
    % That's not a problem for the CDF, but it causes convergence problems
    % for ICDF calculation, so use a cutoff large enough to make the jumps
    % smaller than the convergence criterion.
    cutoff = max(cutoff,6);
else
    % Other kernels have a fixed finite width.  Ignore any requested
    % truncation for these kernels; it would cause convergence problems if
    % smaller than the kernel width, and would have no effect if larger.
    cutoff = kernelcutoff;
end

% If there are any gaps in the data wide enough to create regions of
% exactly zero density, include points at the edges of those regions
% in the grid, to make sure a linear interpolation smooth of the gridded
% CDF captures them as constant
halfwidth = cutoff*u;
gap = find(diff(sy) > 2*halfwidth);
if ~isempty(gap)
    sy = sy(:)';
    xi = sort([xi, sy(gap)+halfwidth, sy(gap+1)-halfwidth]);
    [Fi,xi,u] = compute_pdf_cdf(xi,true,m,L,U,weight,kernel_c,...
                                cutoff,iscdf_c,u,ty,foldpoint,d,false,bdycorr);
end

% Find any regions where the CDF is constant, these will cause problems
% inverse interpolation for x at p
t = (diff(Fi) == 0);
if any(t)
    % Remove interior points in constant regions, they're unnecessary
    s = ([false t] & [t false]);
    Fi(s) = [];
    xi(s) = [];
    % To make Fi monotonic, nudge up the CDF value at the end of each
    % constant region by the smallest amount possible.
    t = 1 + find(diff(Fi) == 0);
    Fi(t) = Fi(t) + eps(Fi(t));
    % If the CDF at the point following is that same value, just remove
    % the nudge.
    if (t(end) == length(Fi)), t(end) = []; end
    s = t(Fi(t) >= Fi(t+1));
    Fi(s) = [];
    xi(s) = [];
end

% ----------
function x = transform(y,L,U,d,bdycorr)
if isequal(bdycorr,'log')
    idx_unbounded = find(L==-Inf & U==Inf);
    idx_positive = find(L==0 & U==Inf);
    idx_bounded = setdiff(1:d,[idx_unbounded,idx_positive]);
    
    x = zeros(size(y));
    if any(idx_unbounded)
        x(:,idx_unbounded) = y(:,idx_unbounded);
    end
    if any(idx_positive)
        x(:,idx_positive) = log(max(0,y(:,idx_positive)));
    end
    if any(idx_bounded)
        y(:,idx_bounded) = bsxfun(@max,L(:,idx_bounded),bsxfun(@min,y(:,idx_bounded),U(:,idx_bounded)));
        
        x(:,idx_bounded) = log(bsxfun(@minus,y(:,idx_bounded,:),L(:,idx_bounded))) -...
            log(bsxfun(@minus,U(:,idx_bounded),y(:,idx_bounded)));
        [i,j]=find(y(:,idx_bounded)==Inf);
        x(i,idx_bounded(j)) = Inf;
    end
else
    x = y;
end

function y = untransform(x,L,U)
if L==-Inf && U==Inf   % unbounded support
    y = x;
elseif L==0 && U==Inf  % positive support
    y = exp(x);
else                   % finite support [L, U]
    t = x<0;
    y = x;
    y(t) = (U*exp(x(t))+L) ./ (exp(x(t))+1);
    t = ~t;
    y(t) = (U+L*exp(-x(t))) ./ (1+exp(-x(t)));
end

function f = untransform_f(f,L,U,xi)
if L==0 && U==Inf   % positive support
    f = bsxfun(@times,f,1./xi');
elseif U<Inf        % bounded support
    tf = (U-L) ./ ((xi-L) .* (U-xi));
    f = bsxfun(@times,f,tf');
end

function fc = compute_CDFreduction(L,U,u,halfwidth,n,ty,weight,kernel)
jstart = 1;
jend = 1;
hi = L + halfwidth;
while(ty(jend)<=hi && jend<n)
    jend = jend+1;
end
nearby = jstart:jend;
z = (L - ty(nearby))/u;
zleft = (ty(nearby) - L)/u;
zright = (L + ty(nearby) -2*U)/u;
if jend == n
    fc = weight(nearby) * (feval(kernel,z)+feval(kernel,zleft)+feval(kernel,zright));
else
    fc = weight(nearby) * (feval(kernel,z)+feval(kernel,zleft));
    fc = fc + sum(weight(jend+1:end));
end

function fc = compute_CDFreduction_mv(L,U,u,ty,kernel)
z = (L - ty(:))/u;
zleft = (ty(:) - L)/u;
zright = (L + ty(:) -2*U)/u;
fc = feval(kernel,z)+feval(kernel,zleft)+feval(kernel,zright);