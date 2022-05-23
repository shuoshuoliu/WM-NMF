function [fout,xout,u] = statkscompute(ftype,xi,xispecified,m,u,L,U,weight,cutoffin,...
                                     kernelname,ty,yData,foldpoint,maxp,varargin)
%#codegen
%STATKSCOMPUTE Perform computations for kernel smoothing density
%estimation.

%   Copyright 2018 The MathWorks, Inc.

coder.inline('always');

if nargin < 15
    d = 1;
    bdycorr = 'log';
else
    d = varargin{1};    
    bdycorr = varargin{2};
    if isempty(bdycorr)
        bdycorr = 'log';
    end
end

[kernel,iscdf,kernelcutoff,kernelname,ftype] = statkskernelinfo(ftype,kernelname);

if isempty(cutoffin)
    cutoff = kernelcutoff;
else
    cutoff = cutoffin;
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
                                    cutoff,iscdf,u,ty,foldpoint,d,bdycorr);

    % If another function based on the cdf, compute it now
    if isequal(ftype,'survivor')
        fout = 1-fout;
    elseif isequal(ftype,'cumhazard')
        fout = 1-fout;
        t = (fout>0);
        fout(~t) = coder.internal.nan;
        fout(t) = -log(fout(t));
    end
end

end


% -----------------------------
function [fout,xout,u]=compute_pdf_cdf(xiIn,xispecified,m,L,U,weight,...
                          kernel,cutoff,iscdf,u,ty,foldpoint,d,bdycorr)

coder.inline('always');

if coder.target('MATLAB')
    NONFINITES = true;
else
    NONFINITES = eml_option('NonFinitesSupport');
end

foldwidth = min(cutoff,3);
issubdist = isfinite(foldpoint);
if ~xispecified
    xin = compute_default_xi(ty,foldwidth,issubdist,m,u,U,L,d,bdycorr);
else
    xin = cast(xiIn,'like',ty);
end

% Compute transformed values of evaluation points that are in bounds
if d == 1
    fout = zeros(size(xin),'like',ty); % f has the same orientation as xi
    xisize = numel(xin);
else
    xisize = size(xin,1);
    fout = zeros(xisize,1,'like',ty);
end
if iscdf && all(isfinite(U))
	li = bsxfun(@ge,xin,U);
    logicind = all(li, 2);
    fout(logicind) = sum(weight);
end
xout = xin;


if d == 1
    xii = xin(:);
else
    xii = xin;
end

[nx,dx] = size(xii);
if coder.internal.isConst(nx) && coder.internal.isConst(dx)
    coder.varsize('xi',[nx+1,dx],[1,0]);
    coder.varsize('inbounds',[nx+1,1],[1,0]);
    coder.varsize('f',[1,nx+1],[0,1]);
else    
    coder.varsize('xi','inbounds','f');
end

if all(L==-coder.internal.inf) && all(U==coder.internal.inf)   % unbounded support
    inbounds = true(xisize,1);
    xi = xii;
elseif all(L==0) && all(U==coder.internal.inf)  % positive support
    inbounds = all(xii>0,2);
    xi = xii(inbounds,:);
else % finite support [L, U]
    inbounds = all(bsxfun(@gt,xii,L),2) & all(bsxfun(@lt,xii,U),2);
    xi = xii(inbounds,:);
end
txi = transform(xi,L,U,d,bdycorr);
if d == 1
    foldpoint = transform(foldpoint,L,U,1,bdycorr);
end

% If the density is censored at the end, add new points so that we can fold
% them back across the censoring point as a crude adjustment for bias.
if issubdist && NONFINITES
    foldbound = foldpoint - foldwidth*u;
    needfold = (txi >= foldbound);
    txifold = (2*foldpoint) - txi(needfold);
    nfold = sum(needfold);
else
    nfold = 0;
    % Codegen: to make txifold and needfold defined on all execution paths
    txifold = cast(0,'like',ty); 
    needfold = true(1,d);
end

if isempty(xi)
    f = cast(xi(:),'like',ty);
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
                xi(end+1) = xi(end); %#ok<EMGRO>
                f(end+1) = 0; %#ok<EMGRO>
                inbounds(end+1) = true; %#ok<EMGRO>
            end
        end
    end
end

if iscdf
    % Guard against roundoff.  Lower boundary of 0 should be no problem.
    f = min(1,f);
else
    f = f ./ prod(u);
end
fout(inbounds) = f;
if d==1
    xout(inbounds) = xi;
end

end

% -----------------------------
function xi = compute_default_xi(ty,foldwidth,issubdist,m,u,U,L,d,bdycorr)
% Get XI values at which to evaluate the density

coder.inline('always');

if coder.target('MATLAB')
    NONFINITES = true;
else
    NONFINITES = eml_option('NonFinitesSupport');
end

% Compute untransformed values of lower and upper evaluation points
ximin = min(ty,[],1) - foldwidth*u;
if issubdist && NONFINITES
    ximax = max(ty,[],1);
else
    ximax = max(ty,[],1) + foldwidth*u;
end

if isequal(bdycorr,'log')
    for i =1:coder.internal.indexInt(d)
        ximin(i) = untransform(ximin(i),L(i),U(i));
        ximax(i) = untransform(ximax(i),L(i),U(i));
    end
else
    for i =1:coder.internal.indexInt(d)
        ximin(i) = max(ximin(i),L(i));
        ximax(i) = min(ximax(i),U(i));
    end
end

if d == 1
    xi = linspace(ximin(1), ximax(1), m);
else
    xt1 = linspace(ximin(1), ximax(1), m);
    xt2 = linspace(ximin(2), ximax(2), m);
    [x1,x2] = meshgrid(xt1,xt2);
    xi = [x1(:) x2(:)];
end

end


% -----------------------------
function f = dokernel(iscdf,txi,ty,u,weight,kernel,cutoff,d,L,U,xi,bdycorr)
% Now compute density estimate at selected points

coder.inline('always');

if coder.target('MATLAB')
    NONFINITES = true;
else
    NONFINITES = eml_option('NonFinitesSupport');
end

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

if NONFINITES
    lTimes2 = 2*L;
    uTimes2 = 2*U;
else
    lTimes2 = max(-coder.internal.inf,2*L);
    uTimes2 = min(coder.internal.inf,2*U);
end

if n*m<=blocksize && ~iscdf
    % For small problems, compute kernel density estimate in one operation
    ftemp = ones(n,m,'like',ty);
    for i = 1:coder.internal.indexInt(d)
        z = zeros(n,m,'like',ty);        
        for ni = 1:coder.internal.indexInt(n)
            for nj = 1:coder.internal.indexInt(m)
                z(ni,nj) = (txi(nj,i) - ty(ni,i))/u(i);
            end
        end
        if reflectionPDF
            zleft = zeros(n,m,'like',ty);
            zright = zeros(n,m,'like',ty);
            for ni = 1:coder.internal.indexInt(n)
                for nj = 1:coder.internal.indexInt(m)
                    zleft(ni,nj) = (ty(ni) + txi(nj) - lTimes2(i))/u(i);
                    zright(ni,nj) = (ty(ni) + txi(nj) - uTimes2(i))/u(i);
                end
            end
            f =  kernelsum(kernel,z,zleft,zright);
        else
            f =  ksdkernel(kernel,z);
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
    f = zeros(1,m,'like',ty);

    if isinf(cutoff)
        if reflectionCDF
            fc = compute_CDFreduction(L,U,u,coder.internal.inf,n,ty,weight,kernel);
        end
        for k=1:coder.internal.indexInt(m)
            % Sum contributions from all
            z = (txi(k)-ty)/u;
            
            if reflectionPDF
                zleft = (txi(k)+ty-lTimes2)/u;
                zright = (txi(k)+ty-uTimes2)/u;                
                f(k) = weight * kernelsum(kernel,z,zleft,zright);
            elseif reflectionCDF
                zleft = (txi(k)+ty-lTimes2)/u;
                zright = (txi(k)+ty-uTimes2)/u;                
                fk = weight * kernelsum(kernel,z,zleft,zright);
                f(k) = fk - fc;
            else
                f(k) = weight * ksdkernel(kernel,z);
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
        
        for k=1:coder.internal.indexInt(m)
            % Find nearby data points for current evaluation point
            lo = stxi(k) - halfwidth;
            while(ty(jstart)<lo(1) && jstart<n)
                jstart = jstart+1;
            end
            hi = stxi(k) + halfwidth;
            jend = max(jend,jstart);
            while(ty(jend)<=hi(1) && jend<n)
                jend = jend+1;
            end
            nearby = jstart:jend;

            % Sum contributions from these points
            z = (stxi(k)-ty(nearby))/u;
            if reflectionPDF
                zleft = (stxi(k)+ty(nearby)-lTimes2)/u;
                zright = (stxi(k)+ty(nearby)-uTimes2)/u;                
                fk = weight(nearby) * kernelsum(kernel,z,zleft,zright);
            elseif reflectionCDF
                zleft = (stxi(k)+ty(nearby)-lTimes2)/u;
                zright = (stxi(k)+ty(nearby)-uTimes2)/u;
                fk = weight(nearby) * ksdkernel(kernel,z);
                fk = fk + sum(weight(1:jstart-1));
                if jstart == 1
                    fk = fk + weight(nearby) * ksdkernel(kernel,zleft);
                    fk = fk + sum(weight(jend+1:end));
                else
                    fk = fk + sum(weight);
                end
                if jend == n
                    fk = fk + weight(nearby) * ksdkernel(kernel,zright);
                end
                fk = fk - fc;
            elseif ~iscdf
                fk = weight(nearby) * ksdkernel(kernel,z);
            elseif iscdf
                fk = weight(nearby) * ksdkernel(kernel,z);
                fk = fk + sum(weight(1:jstart-1));
            end
            f(k) = fk;
        end

        % Restore original x order
        f(idx) = f;
    end
    
    if needUntransform
        f = untransform_f(f,L(1),U(1),xi(:,1));
    end

else % d > 1
    % Calculate reduction for reflectionCDF
    if reflectionCDF
        cutoff = coder.internal.inf;
        fc = zeros(n,d,'like',ty);
        for j = 1:coder.internal.indexInt(d)
            fc(:,j) = compute_CDFreduction_mv(L(j),U(j),u(j),ty(:,j),kernel);
        end
    end
        
    if isinf(cutoff)
        f = zeros(1,m,'like',ty);
        for i = 1:coder.internal.indexInt(m)
            ftemp = ones(n,1,'like',ty);
            for j = 1:coder.internal.indexInt(d)
                z = (txi(i,j) - ty(:,j))./u(j);
                if reflectionPDF
                    zleft = (txi(i,j) + ty(:,j)-lTimes2(j))./u(j);
                    zright = (txi(i,j) + ty(:,j)-uTimes2(j))./u(j);
                    fk = kernelsum(kernel,z,zleft,zright);
                elseif reflectionCDF
                    zleft = (txi(i,j) + ty(:,j)-lTimes2(j))./u(j);
                    zright = (txi(i,j) + ty(:,j)-uTimes2(j))./u(j);
                    fk = kernelsum(kernel,z,zleft,zright);
                    fk = fk - fc(:,j);
                else
                    fk = ksdkernel(kernel,z); 
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
        f = zeros(1,m,'like',ty);
        for i = 1:coder.internal.indexInt(m)
            Idx = true(n,1);
            cdfIdx = true(n,1);
            cdfIdx_allBelow = true(n,1);
            for j = 1:coder.internal.indexInt(d)
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
                ftemp = ones(length(nearby),1,'like',ty);
                for k =1:coder.internal.indexInt(d)
                    z = (txi(i,k) - ty(nearby,k))./u(k);
                    if reflectionPDF
                        zleft = (txi(i,k) + ty(nearby,k)-lTimes2(k))./u(k);
                        zright = (txi(i,k) + ty(nearby,k)-uTimes2(k))./u(k);
                        fk = kernelsum(kernel,z,zleft,zright);
                    else
                        fk = ksdkernel(kernel,z);
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

end


% -----------------------------
function [x1,p] = compute_icdf(xi,xispecified,m,u,Lin,Uin,...
    weight,cutoff,kernelname,ty,yData,foldpoint,maxp,d,bdycorr)

coder.inline('always');

if coder.target('MATLAB')
    DMA = true;
else
    DMA = ~strcmpi(eml_option('UseMalloc'),'Off');
end


if isequal(bdycorr, 'log')
    L = -coder.internal.inf; % log correction has no bounds
    U = coder.internal.inf;
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

[kernel_c,iscdf_c] = statkskernelinfo('cdf',kernelname);
[kernel_p,iscdf_p] = statkskernelinfo('pdf',kernelname);


% Get starting values for ICDF(p) by inverse linear interpolation of
% the gridded CDF, plus some clean-up
x1 = interp1(Fi,xi,p);               % interpolate for p in a good range
x1(isnan(x1) & p<min(Fi)) = min(xi); % use lowest x if p>0 too low
x1(isnan(x1) & p>max(Fi)) = max(xi); % use highest x if p<1 too high
x1(p<0 | p>maxp) = coder.internal.nan;              % out of range
x1(p==0) = L;                        % use lower bound if p==0
x1(p==maxp) = U;                     % and upper bound if p==1 or other max

% Now refine the ICDF using Newton's method for cases with 0<p<1
notdone = find(p>0 & p<maxp);
maxiter = 100;
min_dF0 = sqrt(eps(class(p)));
for iter = 1:coder.internal.indexInt(maxiter)
    if isempty(notdone), break; end
    x0 = x1(notdone);

    % Compute cdf and derivative (pdf) at this value
    F0 = compute_pdf_cdf(x0,true,m,L,U,weight,kernel_c,cutoff,...
                         iscdf_c,u,ty,foldpoint,d,bdycorr);
    dF0 = compute_pdf_cdf(x0,true,m,L,U,weight,kernel_p,cutoff,...
                          iscdf_p,u,ty,foldpoint,d,bdycorr);
                      
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
    if DMA
        p_str = sprintf('%f',p( notdone(1)));
        coder.internal.warning('stats:ksdensity:NoConvergence', p_str);
    else
        coder.internal.warning('stats:ksdensity:NoConvergenceLess');
    end
end

end


% -----------------------------
function [Fi,xiout,cutoff,u] = compute_initial_icdf(m,u,L,U,weight,cutoff,...
    kernelname,ty,yData,foldpoint,d,bdycorr)
% To get starting x values for the ICDF evaluated at p, first create a
% grid xi of values spanning the data on which to evaluate the CDF

coder.inline('always');
n = numel(yData);
if coder.internal.isConst(n)
    coder.varsize('gap',[n-1,1],[1,0]);
    coder.varsize('xi',[1,100+2*n-2],[0,1]);
end

sy = sort(yData);
xi = linspace(sy(1), sy(end), 100);

% Estimate the CDF on the grid
[kernel_c,iscdf_c,kernelcutoff] = statkskernelinfo('cdf',kernelname);

[Fi,xi,u] = compute_pdf_cdf(xi,true,m,L,U,weight,kernel_c,cutoff,...
                            iscdf_c,u,ty,foldpoint,d,bdycorr);

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
    xiout = sort([xi, sy(gap)+halfwidth, sy(gap+1)-halfwidth]);
    [Fi,xiout,u] = compute_pdf_cdf(xiout,true,m,L,U,weight,kernel_c,...
                                cutoff,iscdf_c,u,ty,foldpoint,d,bdycorr);
else
    xiout = xi;
end

% Find any regions where the CDF is constant, these will cause problems
% inverse interpolation for x at p
t = (diff(Fi) == 0);
if any(t)
    % Remove interior points in constant regions, they're unnecessary
    s = ([false t] & [t false]);
    Fi(s) = [];
    xiout(s) = [];
    % To make Fi monotonic, nudge up the CDF value at the end of each
    % constant region by the smallest amount possible.
    t = 1 + find(diff(Fi) == 0);
    Fi(t) = Fi(t) + eps(Fi(t));
    % If the CDF at the point following is that same value, just remove
    % the nudge.
    if (t(end) == length(Fi))
        t(end) = []; 
    end
    s = t(Fi(t) >= Fi(t+1));
    Fi(s) = [];
    xiout(s) = [];
end

end


% ----------
function x = transform(y,L,U,d,bdycorr)

coder.inline('always');
[n,~] = size(y);
x = y;
if isequal(bdycorr,'log')
    for i = 1:coder.internal.indexInt(d)
        if L(i)==-coder.internal.inf && U(i)==coder.internal.inf % unbounded support
            continue;
        elseif L(i)==0 && U(i)==coder.internal.inf % positive support
            x(:,i) = log(max(0,y(:,i)));
        else % finite support [L, U]
            y(:,i) = max(L(:,i),min(y(:,i),U(:,i)));
            
            x(:,i) = log(y(:,i)-L(:,i)) - log(U(:,i)-y(:,i));
            
            for j=1:coder.internal.indexInt(n)
               if y(j,i)== coder.internal.inf
                   x(j,i) = coder.internal.inf;
               end
            end
            
        end
    end    
end

end


function y = untransform(x,L,U)

coder.inline('always');
if L==-coder.internal.inf && U==coder.internal.inf   % unbounded support
    y = x;
elseif L==0 && U==coder.internal.inf  % positive support
    y = exp(x);
else                   % finite support [L, U]
    t = x<0;
    y = x;
    y(t) = (U*exp(x(t))+L) ./ (exp(x(t))+1);
    t = ~t;
    y(t) = (U+L*exp(-x(t))) ./ (1+exp(-x(t)));
end
end


function f = untransform_f(f,L,U,xi)

coder.inline('always');
if L==0 && U==coder.internal.inf   % positive support
    f = bsxfun(@times,f,1./xi');
elseif U<coder.internal.inf        % bounded support
    tf = (U-L) ./ ((xi-L) .* (U-xi));
    f = bsxfun(@times,f,tf');
end

end


function fc = compute_CDFreduction(L,U,u,halfwidth,n,ty,weight,kernel)

coder.inline('always');
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
    fc = weight(nearby) * kernelsum(kernel,z,zleft,zright);
else
    fc = weight(nearby) * (ksdkernel(kernel,z) + ksdkernel(kernel,zleft));
    fc = fc + sum(weight(jend+1:end));
end

end


function fc = compute_CDFreduction_mv(L,U,u,ty,kernel)

coder.inline('always');
z = (L - ty(:))/u;
zleft = (ty(:) - L)/u;
zright = (L + ty(:) -2*U)/u;
fc = kernelsum(kernel,z,zleft,zright);

end 

function f = kernelsum(kernel,z,zleft,zright)

coder.inline('always');
f = ksdkernel(kernel,z) + ksdkernel(kernel,zleft) + ksdkernel(kernel,zright);

end

function [kernel,iscdf,kernelcutoff,kernelname,ftype] =...
    statkskernelinfo(ftype,kernelname)
%#codegen
%STATKSKERNELINFO Obtain kernel information from name.

coder.inline('always');

% Set a flag indicating we are to compute the cdf; later on
% we may transform to another function that is a transformation
% of the cdf
iscdf = isequal(ftype,'cdf') | isequal(ftype,'survivor') ...
    | isequal(ftype,'cumhazard');

if isempty(kernelname)
    if iscdf
        kernel = 'cdf_nl';
    else
        kernel = 'normal';
    end
    kernelname = 'normal';
    kernelcutoff = 4;
elseif ischar(kernelname)
    % If this is an abbreviation of our own methods, expand the name now.
    % If the string matches the start of both variants of the Epanechnikov
    % spelling, that is not an error so pretend it matches just one.
    
    switch lower(kernelname)
        case 'normal'
            if iscdf
                kernel = 'cdf_nl';
            else
                kernel = 'normal';
            end
            kernelcutoff = 4;
        case 'epanechnikov'
            if iscdf
                kernel = 'cdf_ep';
            else
                kernel = 'epanechnikov';
            end
            kernelcutoff = sqrt(5);
        case 'box'
            if iscdf
                kernel = 'cdf_bx';
            else
                kernel = 'box';
            end
            kernelcutoff = sqrt(3);
        case 'triangle'
            if iscdf
                kernel = 'cdf_tr';
            else
                kernel = 'triangle';
            end
            kernelcutoff = sqrt(6);
            
        otherwise % custom kernel specified by name
            if isequal(ftype,'icdf')
                coder.internal.errorIf(true,'stats:ksdensity:IcdfNotAllowed');
            end
            kernel = coder.const(kernelname);
            kernelcutoff = coder.internal.inf;
    end

elseif isa(kernelname,'function_handle') % custom kernel
    coder.internal.errorIf(isequal(ftype,'icdf'),'stats:ksdensity:IcdfNotAllowed');
    
    kernelStr = func2str(kernelname);
    kernel = coder.const(kernelStr);
    kernelcutoff = coder.internal.inf;
else
    coder.internal.errorIf(true,'stats:ksdensity:BadKernel');
end

end

function f = ksdkernel(kernelname,z)

coder.inline('always');
switch kernelname
    
    case 'normal'
        f = normal(z);
        
    case 'epanechnikov'
        f = epanechnikov(z);
        
    case 'box'
        f = box(z);
        
    case 'triangle'
        f = triangle(z);
        
    case 'cdf_nl'
        f = cdf_nl(z);
        
    case 'cdf_ep'
        f = cdf_ep(z);
        
    case 'cdf_bx'
        f = cdf_bx(z);
        
    case 'cdf_tr'
        f = cdf_tr(z);        
        
    otherwise        
        kernelfh = str2func(kernelname);
        f = kernelfh(z);        
end

end


% -----------------------------
% The following are functions that define smoothing kernels k(z).
% Each function takes a single input Z and returns the value of
% the smoothing kernel.  These sample kernels are designed to
% produce outputs that are somewhat comparable (differences due
% to shape rather than scale), so they are all probability
% density functions with unit variance.
%
% The density estimate has the form
%    f(x;k,h) = mean over i=1:n of k((x-y(i))/h) / h

function f = normal(z)
%NORMAL Normal density kernel.

coder.inline('always');
f = exp(-0.5 * z .^2) ./ sqrt(2*pi);

end


function f = epanechnikov(z)
%EPANECHNIKOV Epanechnikov's asymptotically optimal kernel.

coder.inline('always');
a = sqrt(5);
z = max(-a, min(z,a));
f = max(0,.75 * (1 - .2*z.^2) / a);

end


function f = box(z)
%BOX    Box-shaped kernel

coder.inline('always');
a = cast(sqrt(3),'like',z);
f = (abs(z)<=a) ./ (2 * a);

end


function f = triangle(z)
%TRIANGLE Triangular kernel.

coder.inline('always');
a = sqrt(6);
z = abs(z);
% In case z is Inf
f = zeros(size(z),'like',z);
indomain = z<=a;
f(indomain) = (1 - z(indomain)/a) / a;

end

% -----------------------------
% The following are functions that define cdfs for smoothing kernels.

function f = cdf_nl(z)
%CDF_NL Normal kernel, cdf version

coder.inline('always');
f = normcdf(z);

end


function f = cdf_ep(z)
%CDF_EP Epanechnikov's asymptotically optimal kernel, cdf version

coder.inline('always');
a = sqrt(5);
z = max(-a, min(z,a));
f = ((z+a) - (z.^3+a.^3)/15) * 3 / (4*a);

end


function f = cdf_bx(z)
%CDF_BX Box-shaped kernel, cdf version

coder.inline('always');
a = sqrt(3);
f = max(0, min(1,(z+a)/(2*a)));

end


function f = cdf_tr(z)
%CDF_TR Triangular kernel, cdf version

coder.inline('always');
a = sqrt(6);
denom = 12;  % 2*a^2
f = zeros(size(z),'like',z);            % -Inf < z < -a
t = (z>-a & z<0);
f(t) = (a + z(t)).^2 / denom;           % -a < z < 0
t = (z>=0 & z<a);
f(t) = .5 + z(t).*(2*a-z(t)) / denom;   % 0 < z < a
t = (z>a);
f(t) = 1;                               % a < z < Inf


end