classdef RicianDistribution < prob.ToolboxFittableParametricDistribution
%RicianDistribution Rician probability distribution.
%    An object of the RicianDistribution class represents a Rician
%    probability distribution with specific values of the parameters S and
%    SIGMA. This distribution object can be created directly using the
%    MAKEDIST function or fit to data using the FITDIST function.
%
%    RicianDistribution methods:
%       cdf                   - Cumulative distribution function
%       icdf                  - Inverse cumulative distribution function
%       iqr                   - Interquartile range
%       mean                  - Mean
%       median                - Median
%       negloglik             - Negative log likelihood function
%       paramci               - Confidence intervals for parameters
%       pdf                   - Probability density function
%       proflik               - Profile likelihood function
%       random                - Random number generation
%       std                   - Standard deviation
%       truncate              - Truncation distribution to an interval
%       var                   - Variance
%
%    RicianDistribution properties:    
%       DistributionName      - Name of the distribution
%       s                     - Value of the s parameter
%       sigma                 - Value of the sigma parameter
%       NumParameters         - Number of parameters
%       ParameterNames        - Names of parameters
%       ParameterDescription  - Descriptions of parameters
%       ParameterValues       - Vector of values of parameters
%       Truncation            - Two-element vector indicating truncation limits
%       IsTruncated           - Boolean flag indicating if distribution is truncated
%       ParameterCovariance   - Covariance matrix of estimated parameters
%       ParameterIsFixed      - Two-element boolean vector indicating fixed parameters
%       InputData             - Structure containing data used to fit the distribution
%
%    See also fitdist, makedist.

%    Copyright 2012-2018 The MathWorks, Inc.

    properties
    end
    properties(Dependent=true)
%S  The S property represents the noncentrality parameter S of the Rician
%distribution.
%
%    See also ParameterValues.
        s

%SIGMA  The SIGMA property represents the scale parameter SIGMA of the
%Rician distribution.
%
%    See also ParameterValues.
        sigma
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameRician'));

%NumParameter Number of parameters.
%    NumParameters is the number of parameters in the distribution.
%
%    See also ParameterValues.
        NumParameters = 2;
        
%ParameterNames Parameter names.
%    ParameterNames is a cell array of strings containing the names of the
%    parameters of the probability distribution.
%
%    See also ParameterValues, ParameterDescription.
        ParameterNames = {'s' 'sigma'};

%ParameterDescription Parameter description.
%    ParameterNames is a cell array of strings containing short
%    descriptions of the parameters of the probability distribution.
%
%    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionNoncentrality')) ... 
                                getString(message('stats:probdists:ParameterDescriptionScale'))};
    end
    
    properties(GetAccess='public',SetAccess='protected')
%ParameterValues Parameter values.
%    ParameterVales is a vector containing the values of the parameters of
%    the probability distribution.
%
%    See also S, SIGMA.
        ParameterValues
    end
    
    methods(Hidden)
        function pd = RicianDistribution(s,sigma)
            if nargin==0
                s = 1;
                sigma = 1;
            end
            checkargs(s,sigma)
          
            pd.ParameterValues = [s sigma];
            pd.ParameterIsFixed = [true true];
            pd.ParameterCovariance = zeros(pd.NumParameters);
        end
    end
    methods
        function m = mean(this)
            requireScalar(this)
            if this.IsTruncated
                m = truncatedMoment(this,1);
                return
            end
            t = .5 .* (this.s./this.sigma).^2;
            m = this.sigma.*sqrt(.5.*pi) .* ((1+t).*besseli(0,.5.*t,1) + t.*besseli(1,.5.*t,1));
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            m = mean(this);
            v = 2.*this.sigma.^2 + this.s.^2 - m.^2;
        end
    end
    methods
        function this = set.s(this,s)
            checkargs(s,this.sigma)
            this.ParameterValues(1) = s;
            this = invalidateFit(this);
        end
        function s = get.s(this)
            s = this.ParameterValues(1);
        end
        function this = set.sigma(this,sigma)
            checkargs(this.s,sigma)
            this.ParameterValues(2) = sigma;
            this = invalidateFit(this);
        end
        function sigma = get.sigma(this)
            sigma = this.ParameterValues(2);
        end
    end
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the Rician distribution to data. You
%    should call the FITDIST function instead of calling this method
%    directly.
%
%    See also FITDIST.
            [x,cens,freq,opt] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            p = ricefit(x,0.05,cens,freq,opt);
            [nll,cov] = ricelike(p,x,cens,freq);
            pd = prob.RicianDistribution.makeFitted(p,nll,cov,x,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = ricelike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = ricecdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = ricepdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = riceinv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = ricernd(varargin{:});
        end
        function pd = makeFitted(p,nll,cov,x,cens,freq)
            pd = prob.RicianDistribution(p(1),p(2));
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.RicianDistribution');
            info.name = prob.RicianDistribution.DistributionName;
            info.code = 'rician';
            info.censoring = true;
            info.support = [0 Inf];
            info.islocscale = false;
            info.optimopts = true;
            info.logci = [false true];
            info.plim = [0 realmin;Inf Inf];
        end
    end
end % classdef

function checkargs(s,sigma)
if ~(isscalar(s) && isnumeric(s) && isreal(s) && isfinite(s) && s>=0)
    error(message('stats:probdists:NonnegativeParameter','S'))
end
if ~(isscalar(sigma) && isnumeric(sigma) && isreal(sigma) && isfinite(sigma) && sigma>0)
    error(message('stats:probdists:PositiveParameter','SIGMA'))
end
end
% ==== Rician distribution functions ====

% these distribution functions do not yet handle arrays of parameters

function y = ricepdf(x,s,sigma)
%RICEPDF Rician probability density function (pdf).
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

x(x<0) = 0;
sigsq = sigma.^2;
rsq = (x.^2 + s.^2)./(2.*sigsq);
z = x./sigsq;
expon = rsq - z.*s;
y = z .* exp(-expon) .* besseli(0,z.*s,1);
y(expon > (log(realmax(class(x)))-1)) = 0; % fix up 0*Inf
end

function p = ricecdf(x,s,sigma,uflag)
%RICECDF Rician cumulative distribution function (cdf).
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

x(x<0) = 0;
if nargin > 3
    if ~strcmpi(uflag,'upper')
        error(message('stats:cdf:UpperTailProblem'));
    else
        p = ncx2cdf((x./sigma).^2, 2, (s./sigma).^2,'upper');
    end
else
    p = ncx2cdf((x./sigma).^2, 2, (s./sigma).^2);
end
end

function x = riceinv(p,s,sigma)
%RICEINV Inverse of the Rician cumulative distribution function (cdf).
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

x = sigma .* sqrt(ncx2inv(p, 2, (s./sigma).^2));
end

function r = ricernd(s,sigma,varargin)
%RICERND Random arrays from the Rician distribution.
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

[err, sizeOut] = internal.stats.statsizechk(2,s,sigma,varargin{:});
if err > 0
    error(message('stats:ricernd:InconsistentSizes'));
end

r = sigma .* sqrt(ncx2rnd(2, (s./sigma).^2, sizeOut));
end

function [nlogL,acov] = ricelike(params,data,cens,freq)
%RICELIKE Negative log-likelihood for the Rician distribution.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end

nlogL = rice_nloglf(params, data, cens, freq);
if nargout > 1
    acov = mlecov(params, data, 'nloglf',@rice_nloglf, 'cens',cens, 'freq',freq);
end
end

% ==== Rician fitting functions ====

function [phat,pci] = ricefit(x,alpha,cens,freq,opts)
%RICEFIT Parameter estimates and confidence intervals for Rician data.

if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = zeros(size(x)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); end
if nargin < 5, opts = []; end

if any(x <= 0)
    error(message('stats:ricefit:BadData'));
end

% Moment estimators of the uncensored data as starting point
% E[x.^2] = s.^2 + 2.*sigma.^2
% E[x.^4] = s.^4 + 8.*s.^2.*sigma.^2 + 8.*sigma.^4
xsqunc = x(cens == 0).^2;
meanxsq = mean(xsqunc); meanx4th = mean(xsqunc.^2);
if meanxsq^2 < meanx4th && meanx4th < 2*meanxsq^2
    s4th = 2*meanxsq^2 - meanx4th;
    ssq = sqrt(s4th);
    sigsq = .5*(meanxsq - ssq);
    start = [sqrt(ssq) sqrt(sigsq)];
else
    start = cast([1 1],class(x));
end

% The default options include turning fminsearch's display off.  This
% function gives its own warning/error messages, and the caller can turn
% display on to get the text output from fminsearch if desired.
options = statset(statset('ricefit'), opts);
tolBnd = options.TolBnd;
options = optimset(options);

% Maximize the log-likelihood with respect to mu and sigma.
dfltOptions = struct('DerivativeCheck','off', 'HessMult',[], ...
    'HessPattern',ones(2,2), 'PrecondBandWidth',Inf, ...
    'TypicalX',ones(2,1,class(x)), 'MaxPCGIter',1, 'TolPCG',0.1);
funfcn = {'fungrad' 'ricefit' @rice_nloglf [] []};
[phat, ~, ~, err, output] = dfswitchyard('statsfminbx', ...
                     funfcn, start, [-Inf; tolBnd], [Inf; Inf], ...
                     options, dfltOptions, 1, x, cens, freq);
if (err == 0)
    % fminsearch may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if output.funcCount >= options.MaxFunEvals
        warning(message('stats:ricefit:EvalLimit'));
    else
        warning(message('stats:ricefit:IterLimit'));
    end
elseif (err < 0)
    error(message('stats:ricefit:NoSolution'));
end

% Compute CIs using a normal approximation for phat.
if nargout > 1
    acov = mlecov(phat, x, 'nloglf',@rice_nloglf, 'cens',cens, 'freq',freq);
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';
    pci = norminv([probs probs], [phat; phat], [se; se]);
end
end

function [nll,ngrad] = rice_nloglf(parms, x, cens, freq)
%RICE_NLOGLF Objective function for Rician maximum likelihood.
s = parms(1);
sigma = parms(2);

theta = s/sigma;
z = x./sigma;
ztheta = z.*theta;
bess0 = besseli(0, ztheta, 1); % I0(z.*theta)*exp(-z.*theta)
rsq = (z.^2 + theta.^2)./2;
L = -rsq + log(bess0) + log(z./sigma) + ztheta;
ncen = sum(freq.*cens);
if ncen > 0
    cen = (cens == 1);
    zcen = z(cen);
    Q = marcumq(theta,zcen);
    L(cen) = log(Q);
end
nll = -sum(freq .* L);

if nargout > 1
    bess1 = besseli(1, ztheta, 1);
    dlogbess0 = bess1 ./ bess0;
    dL1 = (-theta + dlogbess0.*z) ./ sigma;
    dL2 = -2 * (1 - rsq + dlogbess0.*ztheta) ./ sigma;
    if ncen > 0
        t = exp(-rsq(cen) + ztheta(cen));
        dQdtheta = zcen.*bess1(cen).*t;
        dQdz = -zcen.*bess0(cen).*t;
        dtheta1 = 1./sigma;
        dtheta2 = -theta./sigma;
        % dz1 = 0;
        dz2 = -zcen./sigma;
        dL1(cen) = dQdtheta.*dtheta1 ./ Q;
        dL2(cen) = (dQdtheta.*dtheta2 + dQdz.*dz2) ./ Q;
    end
    ngrad = -[sum(freq .* dL1) sum(freq .* dL2)];
end
end

function Q = marcumq(a,b)
% Q = MARCUMQ(A,B) returns Marcum's "Q" function.

if isa(a,'single') || isa(b,'single')
   Q = NaN(size(b),'single');
else
   Q = NaN(size(b));
end
Q(a~=Inf & b==0) = 1;
Q(a~=Inf & b==Inf) = 0;
Q(a==Inf & b~=Inf) = 1;
z = (isnan(Q) & a==0 & b~=Inf);
if (any(z))
   Q(z) = exp((-b(z).^2)./2);
end

z = isnan(Q) & ~isnan(a) & ~isnan(b);
if (any(z(:)))
%    aa = (a(z).^2)./2;
   aa = (a.^2)./2;
   bb = (b(z).^2)./2;

   d = exp(-aa);
   h = d;
   f = bb.*exp(-bb);
   k = 1;
   delta = f .* h;
   sum = delta;
   j = (delta > sum.*eps(class(delta)));
   while any(j)
      d = aa.*d./k;
      h = h + d;
      f = bb.*f./(k+1);
      delta = f .* h;
      sum(j) = sum(j) + delta(j);
      j = (delta > sum.*eps(class(delta)));
      k = k + 1;
   end
   Q(z) = 1 - sum;
end
end
