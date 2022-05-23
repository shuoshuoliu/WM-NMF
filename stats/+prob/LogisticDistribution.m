classdef LogisticDistribution < prob.ToolboxFittableParametricDistribution
%LogisticDistribution Logistic probability distribution.
%    An object of the LogisticDistribution class represents a logistic
%    probability distribution with a specific mean MU and scale
%    parameter SIGMA. This distribution object can be created directly
%    using the MAKEDIST function or fit to data using the FITDIST function.
%
%    LogisticDistribution methods:
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
%    LogisticDistribution properties:    
%       DistributionName      - Name of the distribution
%       mu                    - Value of the mu parameter (mean)
%       sigma                 - Value of the sigma parameter (scale)
%       NumParameters         - Number of parameters (2)
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

    properties(Dependent=true)
%MU Mean of logistic distribution
%    The MU property represents the parameter that is the mean of the
%    logistic distribution.
%
%    See also SIGMA.
        mu
        
%SIGMA Scale parameter of logistic distribution
%    The SIGMA property represents the parameter that is the scale
%    parameter of the logistic distribution.
%
%    See also MU.
        sigma
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameLogistic'));

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
        ParameterNames = {'mu' 'sigma'};

%ParameterDescription Parameter description.
%    ParameterNames is a cell array of strings containing short
%    descriptions of the parameters of the probability distribution.
%
%    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionLocation')) ...
                                getString(message('stats:probdists:ParameterDescriptionScale'))};
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
        function pd = LogisticDistribution(mu,sigma)
            if nargin==0
                mu = 0;
                sigma = 1;
            end
            checkargs(mu,sigma)

            pd.ParameterValues = [mu sigma];
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
            v = this.sigma.^2 .* pi.^2 ./ 3;
        end
        function this = set.mu(this,mu)
            checkargs(mu,this.sigma);
            this.ParameterValues(1) = mu;
            this = invalidateFit(this);
        end
        function this = set.sigma(this,sigma)
            checkargs(this.mu,sigma);
            this.ParameterValues(2) = sigma;
            this = invalidateFit(this);
        end
        function mu = get.mu(this)
            mu = this.ParameterValues(1);
        end
        function sigma = get.sigma(this)
            sigma = this.ParameterValues(2);
        end

    end
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the logistic distribution to data.
%    Fitting requires the Statistics and Machine Learning Toolbox. You should call the FITDIST
%    function instead of calling this method directly.
%
%    See also FITDIST.

            [x,cens,freq,opt] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            p = logifit(x,0.05,cens,freq,opt);
            [nll,cov] = logilike(p,x,cens,freq);
            pd = prob.LogisticDistribution.makeFitted(p,nll,cov,x,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = logilike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = logicdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = logipdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = logiinv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = logirnd(varargin{:});
        end
        function pd = makeFitted(p,nll,cov,x,cens,freq)
            pd = prob.LogisticDistribution(p(1),p(2));
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.LogisticDistribution');
            info.name = getString(message('stats:dfittool:NameLogistic'));
            info.code = 'logistic';
            info.hasconfbounds = false;
            info.censoring = true;
            info.islocscale = true;
            info.logci = [false true];
        end
    end
end % classdef

function checkargs(mu,sigma)
if ~(isscalar(mu) && isnumeric(mu) && isreal(mu) && isfinite(mu))
    error(message('stats:probdists:ScalarParameter','MU'))
end
if ~(isscalar(sigma) && isnumeric(sigma) && isreal(sigma) && sigma>=0 && isfinite(sigma))
    error(message('stats:probdists:NonnegativeParameter','SIGMA'))
end
end


function y = logipdf(x, mu, sigma)
%LOGIPDF Logistic probability density function (pdf).
if nargin<2, mu = 0; end
if nargin<3, sigma = 1; end
sigma(sigma <= 0) = NaN;

try
    z = (x-mu)./sigma;
catch me
    error(message('stats:cdf:InputSizeMismatch'));
end

k = (z>350); if any(k), z(k) = -z(k); end % prevent Inf/Inf
y = exp(z) ./ ((1 + exp(z)).^2 .* sigma);
end

function p = logicdf(x,varargin)
%LOGICDF Logistic cumulative distribution function (cdf).
if nargin>1 && strcmpi(varargin{end},'upper')
    uflag = true;
    varargin(end) = [];
elseif nargin>1 && ischar(varargin{end})&& ~strcmpi(varargin{end},'upper')
    error(message('stats:cdf:UpperTailProblem'));
else
    uflag = false;
end
p = locallogicdf(uflag,x,varargin{:});

    function p = locallogicdf(uflag,x,mu,sigma)
        if nargin<3, mu = 0; end
        if nargin<4, sigma = 1; end
        sigma(sigma <= 0) = NaN;
        try
            z = (x-mu)./sigma;
        catch me
            error(message('stats:cdf:InputSizeMismatch'));
        end
        if uflag == true
            p = 1 ./ (1 + exp(z));
        else
            p = 1 ./ (1 + exp(-z));
        end
    end
end

function x = logiinv(p, mu, sigma)
%LOGIINV Inverse of the logistic cumulative distribution function (cdf).
if nargin<2, mu = 0; end
if nargin<3, sigma = 1; end
sigma(sigma <= 0) = NaN;

try
    x = logit(p).*sigma + mu;
catch me
    error(message('stats:cdf:InputSizeMismatch'));
end
end


function r = logirnd(mu, sigma, varargin)
%LOGIRND Random arrays from the logistic distribution.
if nargin<1, mu = 0; end
if nargin<2, sigma = 1; end
sigma(sigma <= 0) = NaN;

[err, sizeOut] = internal.stats.statsizechk(2,mu,sigma,varargin{:});
if err > 0
    error(message('stats:logirnd:InconsistentSizes'));
end

p = rand(sizeOut);
try
    r = logit(p).*sigma + mu;
catch me
    error(message('stats:cdf:InputSizeMismatch'));
end
end


function [nlogL,acov] = logilike(params,data,cens,freq)
%LOGILIKE Negative log-likelihood for the logistic distribution.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end

nlogL = logi_nloglf(params, data, cens, freq);
if nargout > 1
    acov = mlecov(params, data, 'nloglf',@logi_nloglf, 'cens',cens, 'freq',freq);
end
end

% ==== Logistic fitting functions ====

function [phat,pci] = logifit(x,alpha,cens,freq,opts)
%LOGIFIT Parameter estimates and confidence intervals for logistic data.

if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = zeros(size(x)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); end
if nargin < 5, opts = []; end

% Moment estimators as starting point
xunc = x(cens == 0);
start = [mean(xunc) std(xunc).*sqrt(3)./pi];

% The default options include turning statsfminbx's display off.  This
% function gives its own warning/error messages, and the caller can turn
% display on to get the text output from statsfminbx if desired.
options = statset(statset('logifit'), opts);
tolBnd = options.TolBnd;
options = optimset(options);
dfltOptions = struct('DerivativeCheck','off', 'HessMult',[], ...
    'HessPattern',ones(2,2), 'PrecondBandWidth',Inf, ...
    'TypicalX',ones(2,1), 'MaxPCGIter',1, 'TolPCG',0.1);

% Maximize the log-likelihood with respect to mu and sigma.
funfcn = {'fungrad' 'logifit' @logi_nloglf [] []};
[phat, ~, ~, err, output] = dfswitchyard( ...
         'statsfminbx',funfcn, start, [-Inf; tolBnd], [Inf; Inf], ...
                     options, dfltOptions, 1, x, cens, freq);
if (err == 0)
    % statsfminbx may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if output.funcCount >= options.MaxFunEvals
        warning(message('stats:logifit:EvalLimit'));
    else
        warning(message('stats:logifit:IterLimit'));
    end
elseif (err < 0)
    error(message('stats:logifit:NoSolution'));
end

if nargout > 1
    acov = mlecov(phat, x, 'nloglf',@logi_nloglf, 'cens',cens, 'freq',freq);
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';

    % Compute the CI for mu using a normal approximation for muhat.
    pci(:,1) = norminv(probs, phat(1), se(1));

    % Compute the CI for sigma using a normal approximation for
    % log(sigmahat), and transform back to the original scale.
    % se(log(sigmahat)) is se(sigmahat) / sigmahat.
    logsigci = norminv(probs, log(phat(2)), se(2)./phat(2));
    pci(:,2) = exp(logsigci);
end
end

function [nll,ngrad] = logi_nloglf(parms, x, cens, freq)
%LOGI_NLOGLF Objective function for logistic maximum likelihood.
mu = parms(1);
sigma = parms(2);
z = (x - mu) ./ sigma;
logitz = 1 ./ (1 + exp(-z));
clogitz = 1 ./ (1 + exp(z));
logclogitz = log(clogitz);
k = (z > 700); if any(k), logclogitz(k) = z(k); end % fix intermediate overflow

L = z + 2.*logclogitz - log(sigma);
ncen = sum(freq.*cens);
if ncen > 0
    cen = (cens == 1);
    L(cen) = logclogitz(cen);
end
if sigma<0
   nll = Inf;
else
   nll = -sum(freq .* L);
end

if nargout > 1
    t = (2.*logitz - 1) ./ sigma;
    dL1 = t;
    dL2 = z.*t - 1./sigma;
    if ncen > 0
        t = logitz(cen) ./ sigma;
        dL1(cen) = t;
        dL2(cen) = z(cen) .* t;
    end
    ngrad = -[sum(freq .* dL1) sum(freq .* dL2)];
end
end

function logitp = logit(p)
%LOGIT Logistic transformation, handling edge and out of range.
logitp = NaN(size(p));
logitp(p==0) = -Inf;
logitp(p==1) = Inf;
ok = (0<p & p<1);
logitp(ok) = log(p(ok)./(1-p(ok)));
end
