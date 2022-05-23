function [fout,xout,u,plottype,ksinfo] = mvksdensity(yData,varargin)
%MVKSDENSITY Compute multivariate kernel density estimate
%   F = MVKSDENSITY(X,XI,'Bandwidth',BW) computes a probability density
%   estimate of the sample in the N-by-D matrix X, at the values in XI,
%   using the required name-value pair argument value BW for the bandwidth.
%   N is the number of points. D is the number of dimensions. XI is a
%   d-column matrix that specifies the values where the density estimate is
%   to be evaluated. BW is a D-element vector specifying the bandwidth of
%   the  kernel smoothing window in each dimension or a scalar value for
%   all dimensions. F is the vector of density value estimates. The
%   estimation is based on a product Gaussian kernel function.
%
%   F = MVKSDENSITY(...,'PARAM1',val1,'PARAM2',val2,...) specifies
%   parameter name/value pairs to control the density estimation.  Valid
%   parameters are the following:
%
%      Parameter    Value
%      'Kernel'     The type of kernel smoother to use, chosen from among
%                   'normal' (default), 'box', 'triangle', and
%                   'epanechnikov'. The same kernel is applied to each
%                   dimension.
%      'Support'    Either 'unbounded' (default) if the density can extend
%                   over the whole real line, or 'positive' to restrict it
%                   to positive values, or a two-row matrix S giving limits
%                   for the support of the density. S(R,1) is the lower
%                   limit for the Rth dimension, and S(R,2) is the upper
%                   limit for the Rth dimension.
%      'Weights'    Vector of the same length as X, giving the weight to
%                   assign to each X value (default is equal weights).
%      'Function'   The function type to estimate, chosen from among 'pdf',
%                   'cdf', or 'survivor' for the density, cumulative
%                   probability, or survivor, respectively.
%      'BoundaryCorrection'  The method for boundary correction for 'pdf'
%                   kernel density estimation with bounded support. Choices
%                   are 'log' and 'reflection'. The 'log' method first
%                   applies a log transformation to convert bounded data to
%                   be unbounded, then transforms back to the original
%                   scale after performing the density estimation. The
%                   'reflection' method adds the reflection of all points
%                   at the boundaries, then computes the density estimate
%                   inside the bounded range including contributions from
%                   the augmented data. Default is 'log'.
%
%   In place of the kernel functions listed above, you can specify another
%   kernel function by using @ (such as @normpdf) or quotes (such as
%   'normpdf'). MVKSDENSITY calls the function in each dimension with a
%   single argument that is an array containing distances between data
%   values in X and locations in XI where the density is evaluated
%   normalized by the bandwidth in that dimension.
%
%   If the 'support' parameter is 'positive', MVKSDENSITY transforms each
%   dimension of X using a log function, estimates the density of the
%   transformed values, and transforms back to the original scale.  If
%   'support' is a two-row matrix, MVKSDENSITY uses the transformation
%   log((Xi-Li)/(Ui-Xi)) in each dimension, where Xi is the elements of X
%   in the i-th dimension, Li is the lower limit and Ui is the upper limit
%   of the i-th dimension. 'support' can also be a combination of positive,
%   unbounded, and bounded variables specified as [0 -Inf L;Inf Inf U]. The
%   'bandwidth' parameter and U outputs are on the scale of the transformed
%   values.
%
%   Example: generate a mixture of two three-dimensional normal
%   distributions, and compute the estimated density.
%      gridx1 = -0.25:.5:1.25;
%      gridx2 = 0:1:15;
%      gridx3 = 5:1:15;
%      [x1,x2,x3] = ndgrid(gridx1,gridx2,gridx3);
%      x1 = x1(:,:)';
%      x2 = x2(:,:)';
%      x3 = x3(:,:)';
%      xi = [x1(:) x2(:) x3(:)];
%      X = [0+.5*rand(20,1) 5+2.5*rand(20,1) 8+rand(20,1);
%           .75+.25*rand(10,1) 8.75+1.25*rand(10,1) 10+.5*rand(10,1)];
%      f = mvksdensity(X,xi,'bandwidth',[0.17,1,0.54]);
%
%   See also KSDENSITY, MESHGRID, NDGRID.

%   Copyright 2015-2016 The MathWorks, Inc.

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

[yData,~,d,ymin,ymax,xispecified,xi,u,npoints,kernelname,...
    support,weight,cens,cutoff,bdycorr,ftype,plottype,isXChunk] = parse_args(yData,varargin{:});

% Process missing values and do error checking
[bad,~,yData,weight,cens] = statslib.internal.removenanrows(yData,weight,cens);
if isempty(yData)
    error(message('stats:mvksdensity:NanX'));  % all obs missing
elseif (bad==2) || numel(weight)>length(weight)
    error(message('stats:ksdensity:InputSizeMismatchWeight'));
elseif (bad==3) || numel(cens)>length(cens)
    error(message('stats:ksdensity:InputSizeMismatchCensoring'));
end
n = size(yData,1);
weight = standardize_weight(weight,n);
if d==1
    cens = standardize_cens(cens,n);
end

[L,U] = compute_finite_support(support,ymin,ymax,d);
ty = apply_support(yData,L,U,d,bdycorr);

[ty,~,weight,u,foldpoint,maxp] = apply_censoring_get_bandwidth(cens,yData,ty,n,ymax,weight,u,d);

[fout,xout,u] = statkscompute(ftype,xi,xispecified,npoints,u,L,U,weight,cutoff,...
    kernelname,ty,yData,foldpoint,maxp,d,isXChunk,bdycorr);

if nargout>=5
    ksinfo.ty = ty;
    ksinfo.weight = weight;
    ksinfo.foldpoint = foldpoint;
    ksinfo.L = L;
    ksinfo.U = U;
    ksinfo.maxp = maxp;
end

% -----------------------------
function [yData,n,d,ymin,ymax,xispecified,xi,u,npoints,kernelname,...
    support,weight,cens,cutoff,bdycorr,ftype,plottype,isXChunk] = parse_args(yData,varargin)

if ~ismatrix(yData) || isempty(yData)
    error(message('stats:mvksdensity:BadX'));
end

[n,d] = size(yData);

xi = zeros(0,d);
xispecified = false;
if ~isempty(varargin)
    if ~ischar(varargin{1})
        xi = varargin{1};
        varargin(1) = [];
        if ~isempty(xi)
            xispecified = true;
        end
    end
end

ymin = min(yData,[],1);
ymax = max(yData,[],1);

if d==1
    if xispecified && ~isvector(xi)
        error(message('stats:ksdensity:VectorRequired2'));
    end
elseif d==2
    if  xispecified && (~ismatrix(xi)||size(xi,2)~=2)
        error(message('stats:mvksdensity:BadXId',2));
    end
else
    if ~xispecified || ~ismatrix(xi) || size(xi,2)~=d
        error(message('stats:mvksdensity:BadXId',d));
    end
end

% Process additional name/value pair arguments
okargs_n   = {'bandwidth' 'kernel'   'support'   'weights'        'function'  'cutoff' 'boundarycorrection'};
defaults_n = {[]          'normal'   'unbounded'  repmat(1/n,n,1) 'pdf'       []       'log'};
okargs_2   = {'plotfcn'};
defaults_2 = {'surf'};
okargs_1   = { {'npoints','NumPoints'} 'censoring'};
defaults_1 = {          100             false(n,1)};

if d==1
    % Process additional name/value pair arguments
    okargs = [okargs_n okargs_1];
    defaults = [defaults_n defaults_1];
    [u,kernelname,support,weight,ftype,cutoff,bdycorr,npoints,cens,~,extra] = ...
        internal.stats.parseArgs(okargs, defaults, varargin{:});
    plottype = '';
else
    if d==2
        okargs = [okargs_n okargs_2];
        defaults = [defaults_n defaults_2];
        [u,kernelname,support,weight,ftype,cutoff,bdycorr,plottype,~,extra] = ...
            internal.stats.parseArgs(okargs, defaults, varargin{:});
        if ~xispecified
            npoints = 30;
        else
            npoints = [];
        end
        plottype0 = {'contour','plot3','surf','surfc'};
        plottype = internal.stats.getParamVal(plottype,plottype0,'PlotFcn');
    else
        [u,kernelname,support,weight,ftype,cutoff,bdycorr,~,extra] = ...
            internal.stats.parseArgs(okargs_n, defaults_n, varargin{:});
        npoints = [];
        plottype = '';
    end
    cens = false(n,1);
end

% Check for boundary correction method
okbdycorr = {'log','reflection'};
bdycorr = internal.stats.getParamVal(bdycorr,okbdycorr,'BoundaryCorrection');

isXChunk = false;
if ~isempty(extra)
    [u,isXChunk] = internal.stats.parseArgs({'width','isXChunk'},{u,isXChunk},extra{:});
end

% Regularize some vector outputs
if isscalar(weight)
    weight = repmat(weight,n,1);   % scalar -> vector
elseif isrow(weight)
    weight = weight(:);
end
if isrow(cens)
    cens = cens(:);
end

% -----------------------------
function weight = standardize_weight(weight,n)
if isempty(weight)
    weight = ones(1,n);
elseif numel(weight)==1
    weight = repmat(weight,1,n);
elseif numel(weight)~=n || numel(weight)>length(weight)
    error(message('stats:ksdensity:InputSizeMismatchWeight'));
else
    weight = weight(:)';
end
weight = weight / sum(weight);

% -----------------------------
function cens = standardize_cens(cens,n)
if isempty(cens)
    cens = false(1,n);
elseif any(cens(:)~=0 & cens(:)~=1) 
    error(message('stats:ksdensity:BadCensoring'));
elseif numel(cens)~=n || numel(cens)>length(cens)
    error(message('stats:ksdensity:InputSizeMismatchCensoring'));
elseif all(cens)
    error(message('stats:ksdensity:CompleteCensoring'));
end
cens = cens(:);

% -----------------------------
function [L,U] = compute_finite_support(support,ymin,ymax,d)
if isnumeric(support)
    if d==1
        if numel(support)~=2
            error(message('stats:ksdensity:BadSupport1'));
        end
        L = support(1);
        U = support(2);
    else
        if size(support,1)~=2 || size(support,2)~=d
            error(message('stats:mvksdensity:BadSupport',d));
        end
        L = support(1,:);
        U = support(2,:);
    end
    if any(L>=ymin) || any(U<=ymax)
        error(message('stats:ksdensity:BadSupport2'));
    end
elseif ischar(support) && ~isempty(support)
    okvals = {'unbounded' 'positive'};    
    try
        support = validatestring(support,okvals);
    catch
        error(message('stats:ksdensity:BadSupport3'))
    end
    if isequal(support,'unbounded')
        L = -Inf(1,d);
        U = Inf(1,d);
    else
        L = zeros(1,d);
        U = Inf(1,d);
    end
    if isequal(support,'positive') && any(ymin<=0)
        error(message('stats:ksdensity:BadSupport4'))
    end
else
    error(message('stats:ksdensity:BadSupport5'))
end

% -----------------------------
function ty = apply_support(yData,L,U,d,bdycorr)
% Compute transformed values of data
if isequal(bdycorr,'log')
    idx_unbounded = find(L==-Inf & U==Inf);
    idx_positive = find(L==0 & U==Inf);
    idx_bounded = setdiff(1:d,[idx_unbounded,idx_positive]);
    
    ty = zeros(size(yData));
    if any(idx_unbounded)   % unbounded support
        ty(:,idx_unbounded) = yData(:,idx_unbounded);
    end
    if any(idx_positive)    % positive support
        ty(:,idx_positive) = log(yData(:,idx_positive));
    end
    if any(idx_bounded)     % finite support [L, U]
        ty(:,idx_bounded) = log(bsxfun(@minus,yData(:,idx_bounded),L(:,idx_bounded))) -...
            log(bsxfun(@minus,U(:,idx_bounded),yData(:,idx_bounded))); % same as log((x-L)./(U-x))
    end
else
    ty = yData;
end

% -----------------------------
function [ty,ymax,weight,u,foldpoint,maxp] = apply_censoring_get_bandwidth(cens,yData,ty,n,ymax,weight,u,d)

% Deal with censoring
iscensored = any(cens);
if iscensored
    % Compute empirical cdf and create an equivalent weighted sample
    [F,XF] = ecdf(ty, 'censoring',cens, 'frequency',weight);
    weight = diff(F(:))';
    ty = XF(2:end);
    N = sum(~cens);
    ymax = max(yData(~cens),[],1);
    foldpoint = min(yData(cens & yData>=ymax),[],1); % for bias adjustment
    issubdist = ~isempty(foldpoint);  % sub-distribution, integral < 1
    maxp = F(end);
else
    N = n;
    issubdist = false;
    maxp = 1;
end
if ~issubdist
    foldpoint = Inf; % no bias adjustment is needed
end

% Get bandwidth if not already specified
if isempty(u)
    if d>2
        error(message('stats:mvksdensity:BadBandwidth',d));
    else
        if ~iscensored
            % Get a robust estimate of sigma
            sig = mad(ty,1,1) / 0.6745;
        else
            % Estimate sigma using quantiles from the empirical cdf
            Xquant = interp1(F,XF,[.25 .5 .75]);
            if ~any(isnan(Xquant))
                % Use interquartile range to estimate sigma
                sig = (Xquant(3) - Xquant(1)) / (2*0.6745);
            elseif ~isnan(Xquant(2))
                % Use lower half only, if upper half is not available
                sig = (Xquant(2) - Xquant(1)) / 0.6745;
            else
                % Can't easily estimate sigma, just get some indication of spread
                sig = ty(end) - ty(1);
            end
        end
        idx = sig<=0;
        if any(idx)
            sig(idx) = max(ty(:,idx),[],1)-min(ty(:,idx),[],1);
        end
        if sig>0
            % Default window parameter is optimal for normal distribution
            % Scott's rule
            u = sig * (4/((d+2)*N))^(1/(d+4));
        else
            u = ones(1,d);
        end
    end
else
    if isscalar(u)
        u = u*ones(1,d);
    elseif any(size(u(:))~=[d,1])
        error(message('stats:mvksdensity:BadBandwidth',d));
    end
end
