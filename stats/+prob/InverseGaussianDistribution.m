classdef InverseGaussianDistribution < prob.ToolboxFittableParametricDistribution
%InverseGaussianDistribution Inverse Gaussian probability distribution.
%    An object of the InverseGaussianDistribution class represents an inverse
%    Gaussian probability distribution with specific values of the MU and
%    LAMBDA parameters. This distribution object can be created directly
%    using the MAKEDIST function or fit to data using the FITDIST function.
%
%    InverseGaussianDistribution methods:
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
%    InverseGaussianDistribution properties:    
%       DistributionName      - Name of the distribution
%       mu                    - Value of the mu parameter (scale)
%       lambda                - Value of the lambda parameter (shape)
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

%    Copyright 2012-2019 The MathWorks, Inc.

    properties
    end
    properties(Dependent=true)
%MU Scale parameter
%    The MU property represents the scale parameter of the inverse
%    Gaussian distribution.
%
%    See also LAMBDA.
        mu

%LAMBDA Shape parameter
%    The LAMBDA property represents the shape parameter of the inverse
%    Gaussian distribution.
%
%    See also MU.
        lambda
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameInverseGaussian'));

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
        ParameterNames = {'mu' 'lambda'};

%ParameterDescription Parameter description.
%    ParameterNames is a cell array of strings containing short
%    descriptions of the parameters of the probability distribution.
%
%    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionScale')) ...
                                getString(message('stats:probdists:ParameterDescriptionShape'))};
    end
    properties(GetAccess='public',SetAccess='protected')
%ParameterValues Parameter values.
%    ParameterVales is a vector containing the values of the parameters of
%    the probability distribution.
%
%    See also MU, SIGMA.
        ParameterValues
    end
    methods(Hidden)
        function pd = InverseGaussianDistribution(mu,lambda)
            if nargin==0
                mu = 1;
                lambda = 1;
            end
            checkargs(mu,lambda)
            
            pd.ParameterValues = [mu lambda];
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
            m = this.mu;
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            v = this.mu.^3 ./ this.lambda;
        end
        function this = set.mu(this,mu)
            checkargs(mu,this.lambda)
            this.ParameterValues(1) = mu;
            this = invalidateFit(this);
        end
        function this = set.lambda(this,lambda)
            checkargs(this.mu,lambda)
            this.ParameterValues(2) = lambda;
            this = invalidateFit(this);
        end
        function mu = get.mu(this)
            mu = this.ParameterValues(1);
        end
        function lambda = get.lambda(this)
            lambda = this.ParameterValues(2);
        end
    end
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the inverse Gaussian distribution to data.
%    You should call the FITDIST function instead of calling this method directly.
%
%    See also FITDIST.
            [x,cens,freq,opt] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            p = invgfit(x,0.05,cens,freq,opt);
            [nll,cov] = invglike(p,x,cens,freq);
            pd = prob.InverseGaussianDistribution.makeFitted(p,nll,cov,x,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = invglike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = invgcdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = invgpdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = invginv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = invgrnd(varargin{:});
        end
        function pd = makeFitted(p,nll,cov,x,cens,freq)
            pd = prob.InverseGaussianDistribution(p(1),p(2));
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.InverseGaussianDistribution');
            info.name = prob.InverseGaussianDistribution.DistributionName;
            info.code = 'inverse gaussian';
            info.censoring = true;
            info.support = [0 Inf];
        end    end
end % classdef

function checkargs(mu,lambda)
if ~(isscalar(mu) && isnumeric(mu) && isreal(mu) && isfinite(mu) && mu>0)
    error(message('stats:probdists:PositiveParameter','MU'))
end
if ~(isscalar(lambda) && isnumeric(lambda) && isreal(lambda) && isfinite(lambda) && lambda>0)
    error(message('stats:probdists:PositiveParameter','LAMBDA'))
end
end



function y = invgpdf(x, mu, lambda)
%INVGPDF inverse Gaussian probability density function (pdf).
% Inputs need to match in size. This is guaranteed for the common case of
% scalar parameters.
[x,mu,lambda] = internal.stats.checkCompatibleSizes(x,mu,lambda);

mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;

nonpos = (x <= 0);
x(nonpos)= realmin;
y = sqrt(lambda./(2.*pi.*x.^3)) .* exp(-0.5.*lambda.*(x./mu - 2 + mu./x)./mu);
% this would happen automatically for x==0, but generates DivideByZero warnings
y(nonpos) = 0;
end

function p = invgcdf(x, mu, lambda,uflag)
%INVGCDF inverse Gaussian cumulative distribution function (cdf).
[x,mu,lambda] = internal.stats.checkCompatibleSizes(x,mu,lambda);
mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;

nonpos = (x <= 0);
posinf = (x == Inf);
x(nonpos)= realmin;
z1 = (x./mu - 1).*sqrt(lambda./x);
z2 = -(x./mu + 1).*sqrt(lambda./x);
if nargin > 3 
    if ~strcmpi(uflag,'upper')
        error(message('stats:cdf:UpperTailProblem'));
    else
        p = 0.5.*erfc(z1./sqrt(2)) - exp(2.*lambda./mu) .* 0.5.*erfc(-z2./sqrt(2));
    end
else
    p = 0.5.*erfc(-z1./sqrt(2)) + exp(2.*lambda./mu) .* 0.5.*erfc(-z2./sqrt(2));
end
% this would happen automatically for x==0, but generates DivideByZero warnings
if nargin > 3
    p(nonpos) = 1;
    p(posinf) = 0;
else
    p(nonpos) = 0;
    p(posinf) = 1;
end
end

function x = invginv(p, mu, lambda)
%INVGINV Inverse of the inverse Gaussian cumulative distribution function (cdf).
[p,mu,lambda] = internal.stats.checkCompatibleSizes(p,mu,lambda);
mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;

try
    okParams = (0 < mu) & (0 < lambda & lambda < Inf);
    k = (okParams & (0 < p & p < 1));
catch me
    error(message('stats:invginv:InputSizeMismatch'));
end
allOK = all(k(:));

% Fill in NaNs for out of range cases, fill in edges cases when P is 0 or 1.
if ~allOK
    if isa(p,'single') || isa(mu,'single') || isa(lambda,'single')
       x = NaN(size(k),'single');
    else
       x = NaN(size(k));
    end
    x(p == 0 & okParams) = 0;
    x(p == 1 & okParams) = Inf;

    % Remove the bad/edge cases, leaving the easy cases.  If there's
    % nothing remaining, return.
    if any(k(:))
        if numel(p) > 1, p = p(k); end
        if numel(mu) > 1, mu = mu(k); end
        if numel(lambda) > 1, lambda = lambda(k); end
    else
        return;
    end
end

% Newton's Method to find a root of invgcdf(x,mu,lambda) = p
%
% Choose a starting guess for q.  Use quantiles from a lognormal
% distribution with the same mean (==1) and variance (==lambda0) as
% IG(1,lambda0).
lambda0 = lambda./mu;
sigsqLN = log(1./lambda0 + 1);
muLN = -0.5 .* sigsqLN;
q = exp(muLN - sqrt(2.*sigsqLN).*erfcinv(2*p));

% Limit the number of iterations.
maxiter = 500;
reltol = eps(class(q)).^(3/4);

F = invgcdf(q,1,lambda0);
dF = F - p;
for iter = 1:maxiter
    % Compute the Newton step, but limit its size to prevent stepping to
    % negative or infinite values.
    f = invgpdf(q,1,lambda0);
    h = dF ./ f;
    qNew = max(q/10, min(10*q, q - h)); % Avoid taking too large of a step
    
    % Break out of the iteration loop when the relative size of the last step
    % is small for all elements of q.
    done = (abs(h) <= reltol*q);
    if all(done(:))
        q = qNew;
        break
    end
    
    % Check the error, and backstep for those elements whose error did not
    % decrease.  The direction of h is always correct, because f > 0
    dFold = dF;
    for j = 1:25 % If this fails, the outer loop may too
        F = invgcdf(qNew,1,lambda0);
        dF = F - p;
        worse = (abs(dF) > abs(dFold)) & ~done;
        if ~any(worse(:))
            break
        end
        qNew(worse) = (q(worse) + qNew(worse))/2;
    end
    q = qNew;
end

badcdf = (abs(dF./F) > reltol.^(2/3));
if iter>maxiter || any(badcdf(:))   % too many iterations or cdf is too far off
    didnt = find(~done | badcdf, 1, 'first');
    if numel(mu) == 1, mubad = mu; else mubad = mu(didnt); end
    if numel(lambda) == 1, lambdabad = lambda; else lambdabad = lambda(didnt); end
    if numel(p) == 1, pbad = p; else pbad = p(didnt); end
    warning(message('stats:invginv:NoConvergence', num2str( mubad ), num2str( lambdabad ), num2str( pbad )));
end

% Add in the scale factor, and broadcast the values to the correct place if
% need be.
if allOK
    x = q .* mu;
else
    x(k) = q .* mu;
end
end

function r = invgrnd(mu, lambda, varargin)
%INVGRND Random arrays from the inverse Gaussian distribution.
mu(mu <= 0) = NaN;
lambda(lambda <= 0) = NaN;

[err, sizeOut] = internal.stats.statsizechk(2,mu,lambda,varargin{:});
if err > 0
    error(message('stats:invgrnd:InconsistentSizes'));
end

if isscalar(mu), mu = repmat(mu,sizeOut); end

c = mu.*chi2rnd(1,sizeOut);
r = (mu./(2.*lambda)) .* (2.*lambda + c - sqrt(4.*lambda.*c + c.^2));
invert = (rand(sizeOut).*(mu+r) > mu);
r(invert) = mu(invert).^2 ./ r(invert);
end

function [nlogL,acov] = invglike(params,data,cens,freq)
%INVGLIKE Negative log-likelihood for the inverse Gaussian distribution.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end

nlogL = invg_nloglf(params, data, cens, freq);
if nargout > 1
    acov = mlecov(params, data, 'nloglf',@invg_nloglf, 'cens',cens, 'freq',freq);
end
end

% ==== inverse Gaussian fitting functions ====

function [phat,pci] = invgfit(x,alpha,cens,freq,opts)
%INVGFIT Parameter estimates and confidence intervals for inverse Gaussian data.

if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = zeros(size(x)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); end
if nargin < 5, opts = []; end

if any(x <= 0)
    error(message('stats:invgfit:BadData'));
end

ncen = sum(freq.*cens);
if ncen == 0
    xbar = mean(x);
    phat = [xbar 1./mean(1./x - 1./xbar)];

else
    % MLEs of the uncensored data as starting point
    xunc = x(cens == 0); xbarunc = mean(xunc);
    start = [xbarunc 1./mean(1./xunc - 1./xbarunc)];

    % The default options include turning statsfminbx's display off.  This
    % function gives its own warning/error messages, and the caller can
    % turn display on to get the text output from statsfminbx if desired.
    options = statset(statset('invgfit'), opts);
    tolBnd = options.TolBnd;
    options = optimset(options);
    dfltOptions = struct('DerivativeCheck','off', 'HessMult',[], ...
        'HessPattern',ones(2,2), 'PrecondBandWidth',Inf, ...
        'TypicalX',ones(2,1), 'MaxPCGIter',1, 'TolPCG',0.1);

    % Maximize the log-likelihood with respect to mu and lambda.
    funfcn = {'fungrad' 'invgfit' @invg_nloglf [] []};
    [phat, ~, ~, err, output] = dfswitchyard('statsfminbx',...
                    funfcn, start, [tolBnd; tolBnd], [Inf; Inf], ...
                    options, dfltOptions, 1, x, cens, freq);
    if (err == 0)
        % statsfminbx may print its own output text; in any case give something
        % more statistical here, controllable via warning IDs.
        if output.funcCount >= options.MaxFunEvals
            warning(message('stats:invgfit:EvalLimit'));
        else
            warning(message('stats:invgfit:IterLimit'));
        end
    elseif (err < 0)
        error(message('stats:invgfit:NoSolution'));
    end
end

% Compute CIs using a normal approximation for phat.
if nargout > 1
    acov = mlecov(phat, x, 'nloglf',@invg_nloglf, 'cens',cens, 'freq',freq);
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';
    pci = norminv([probs probs], [phat; phat], [se; se]);
end
end

function [nll,ngrad] = invg_nloglf(params, x, cens, freq)
%INVG_NLOGLF Objective function for inverse Gaussian maximum likelihood.

mu = params(1);
lambda = params(2);

L = .5.*log(lambda./(2*pi)) - 1.5.*log(x) - lambda.*(x./mu-1).^2 ./ (2.*x);
ncen = sum(freq.*cens);
if ncen > 0
    cen = (cens == 1);
    xcen = x(cen);
    tmpsqrt = sqrt(lambda./xcen);
    tmpexp = exp(2.*lambda./mu);
    zcen = -(xcen./mu-1) .* tmpsqrt;
    wcen = -(xcen./mu+1) .* tmpsqrt;
    Phizcen = 0.5.*erfc(-zcen./sqrt(2));
    Phiwcen = 0.5.*erfc(-wcen./sqrt(2));
    Scen = Phizcen - tmpexp .* Phiwcen;
    L(cen) = log(Scen);
end
nll = -sum(freq .* L);

if nargout > 1
    dL1 = lambda.*(x-mu)./mu.^3;
    dL2 = 1./(2.*lambda) - (x./mu-1).^2 ./ (2.*x);
    if ncen > 0
        phizcen = exp(-0.5.*zcen.^2)./sqrt(2.*pi);
        phiwcen = exp(-0.5.*wcen.^2)./sqrt(2.*pi);
        dS1cen = (phizcen - tmpexp.*phiwcen).*(xcen./mu.^2).*tmpsqrt ...
                  + 2.*Phiwcen.*tmpexp.*lambda./mu.^2;
        dS2cen = 0.5.*(phizcen.*zcen - tmpexp.*phiwcen.*wcen)./lambda ...
                  - 2.*Phiwcen.*tmpexp./mu;
        dL1(cen) = dS1cen ./ Scen;
        dL2(cen) = dS2cen ./ Scen;
    end
    ngrad = -[sum(freq .* dL1) sum(freq .* dL2)];
end
end
