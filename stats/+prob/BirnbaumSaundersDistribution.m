classdef BirnbaumSaundersDistribution < prob.ToolboxFittableParametricDistribution
%BirnbaumSaundersDistribution BirnbaumSaunders probability distribution.
%    An object of the BirnbaumSaundersDistribution class represents a
%    Birnbaum-Saunders probability distribution with specific values for
%    the BETA and GAMMA parameters. This distribution object can be created
%    directly using the MAKEDIST function or fit to data using the FITDIST
%    function.
%
%    BirnbaumSaundersDistribution methods:
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
%    BirnbaumSaundersDistribution properties:    
%       DistributionName      - Name of the distribution
%       beta                  - Value of the beta parameter (scale)
%       gamma                 - Value of the gamma parameter (shape)
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

    properties(Dependent=true)
%BETA BETA (scale) parameter
%    The BETA property represents the scale parameter of the
%    Birnbaum-Saunders distribution.
%
%    See also GAMMA.
        beta

%GAMMA GAMMA (shape) parameter
%    The GAMMA property represents the shape parameter of the
%    Birnbaum-Saunders distribution.
%
%    See also GETA.
        gamma
    end
   
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameBirnbaumSaunders'));

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
        ParameterNames = {'beta' 'gamma'};

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
%    See also BETA, GAMMA.
        ParameterValues
    end
    
    methods(Hidden)
        function pd = BirnbaumSaundersDistribution(beta,gamma)
            if nargin==0
                beta = 1;
                gamma = 1;
            end
            checkargs(beta,gamma)
            
            pd.ParameterValues = [beta gamma];
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
            m = this.beta .* (0.5 .* this.gamma.^2 + 1);
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            v = (this.beta.*this.gamma).^2 .* (1.25 .* this.gamma.^2 + 1);
        end
        function this = set.beta(this,beta)
            checkargs(beta,this.gamma)
            this.ParameterValues(1) = beta;
            this = invalidateFit(this);
        end
        function this = set.gamma(this,gamma)
            checkargs(this.beta,gamma)
            this.ParameterValues(2) = gamma;
            this = invalidateFit(this);
        end
        function beta = get.beta(this)
            beta = this.ParameterValues(1);
        end
        function gamma = get.gamma(this)
            gamma = this.ParameterValues(2);
        end
    end
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the Birnbaum-Saudners distribution
%    to data. You should call the FITDIST function instead of calling this
%    method directly.
%
%    See also FITDIST.
            [x,cens,freq,opt] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            p = bisafit(x,0.05,cens,freq,opt);
            [nll,cov] = bisalike(p,x,cens,freq);
            pd = prob.BirnbaumSaundersDistribution.makeFitted(p,nll,cov,x,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = bisalike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = bisacdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = bisapdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = bisainv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = bisarnd(varargin{:});
        end
        function pd = makeFitted(p,nll,cov,x,cens,freq)
            pd = prob.BirnbaumSaundersDistribution(p(1),p(2));
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.BirnbaumSaundersDistribution');
            info.name = prob.BirnbaumSaundersDistribution.DistributionName;
            info.code = 'birnbaumsaunders';
            info.censoring = true;
            info.support = [0 Inf];
            info.optimopts = true;
        end
    end
end % classdef

function checkargs(beta,gamma)
if ~(isscalar(beta) && isnumeric(beta) && isreal(beta) && isfinite(beta) && beta>0)
    error(message('stats:probdists:PositiveParameter','BETA'))
end
if ~(isscalar(gamma) && isnumeric(gamma) && isreal(gamma) && gamma>=0 && isfinite(gamma))
    error(message('stats:probdists:PositiveParameter','GAMMA'))
end
end


% ==== Birnbaum-Saunders distribution functions ====

function y = bisapdf(x, beta, gamma)
%BISAPDF Birnbaum-Saunders probability density function (pdf).
[x,beta,gamma] = internal.stats.checkCompatibleSizes(x,beta,gamma);
beta(beta <= 0) = NaN;
gamma(gamma <= 0) = NaN;
x(x<0) = 0;

z = (sqrt(x./beta) - sqrt(beta./x)) ./ gamma;
w = (sqrt(x./beta) + sqrt(beta./x)) ./ gamma;
ynorm = exp(-0.5 .* z.^2) ./ sqrt(2.*pi);
y = ynorm .* w ./ (2.*x);
y(x==0) = 0;
end

function p = bisacdf(x, beta, gamma,uflag)
%BISACDF Birnbaum-Saunders cumulative distribution function (cdf).
beta(beta <= 0) = NaN;
gamma(gamma <= 0) = NaN;
x(x<0) = 0;
if nargin>3 
    if ~strcmpi(uflag,'upper')
        error(message('stats:cdf:UpperTailProblem'));
    else
        z = (-sqrt(x./beta) + sqrt(beta./x)) ./ gamma;
    end
else
   z = (sqrt(x./beta) - sqrt(beta./x)) ./ gamma;
end
p = 0.5 * erfc(-z ./ sqrt(2));
end

function x = bisainv(p, beta, gamma)
%BISAINV Inverse of the Birnbaum-Saunders cumulative distribution function (cdf).
beta(beta <= 0) = NaN;
gamma(gamma <= 0) = NaN;

p(p < 0 | 1 < p) = NaN;
gamz = -sqrt(2).*erfcinv(2*p) .* gamma;
x = 0.25 .* beta .* (gamz + sqrt(4+gamz.^2)).^2;
x(p == 0 & ~isnan(beta) & ~isnan(gamma)) = 0;
% x(p == 1 & ~isnan(beta) & ~isnan(gamma)) = Inf; % get this automatically
end

function r = bisarnd(beta, gamma, varargin)
%BISARND Random arrays from the Birnbaum-Saunders distribution.
beta(beta <= 0) = NaN;
gamma(gamma <= 0) = NaN;

[err, sizeOut] = internal.stats.statsizechk(2,beta,gamma,varargin{:});
if err > 0
    error(message('stats:bisarnd:InconsistentSizes'));
end

plusminus = 2.*(rand(sizeOut)>.5) - 1; % plus or minus one, w.p. 1/2
gamz = gamma.*randn(sizeOut);
r = 0.5.*beta .* (2 + gamz.^2 + plusminus.*gamz.*sqrt(4 + gamz.^2));
end

function [nlogL,acov] = bisalike(params,data,cens,freq)
%BISALIKE Negative log-likelihood for the Birnbaum-Saunders distribution.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end

nlogL = bisa_nloglf(params, data, cens, freq);
if nargout > 1
    acov = mlecov(params, data, 'nloglf',@bisa_nloglf, 'cens',cens, 'freq',freq);
end
end

% ==== Birnbaum-Saunders fitting functions ====

function [phat,pci] = bisafit(x,alpha,cens,freq,opts)
%BISAFIT Parameter estimates and confidence intervals for Birnbaum-Saunders data.

if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = zeros(size(x)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); end
if nargin < 5, opts = []; end

if any(x <= 0)
    error(message('stats:bisafit:BadData'));
end

% Starting points as suggested by Birnbaum and Saunders
xunc = x(cens==0); xbarunc = mean(xunc); xinvbarunc = mean(1./xunc);
start = [sqrt(xbarunc./xinvbarunc) 2.*sqrt(sqrt(xbarunc.*xinvbarunc) - 1)];

% The default options include turning statsfminbx's display off.  This
% function gives its own warning/error messages, and the caller can turn
% display on to get the text output from statsfminbx if desired.
options = statset(statset('bisafit'), opts);
tolBnd = options.TolBnd;
options = optimset(options);
dfltOptions = struct('DerivativeCheck','off', 'HessMult',[], ...
    'HessPattern',ones(2,2), 'PrecondBandWidth',Inf, ...
    'TypicalX',ones(2,1), 'MaxPCGIter',1, 'TolPCG',0.1);

% Maximize the log-likelihood with respect to mu and sigma.
funfcn = {'fungrad' 'bisafit' @bisa_nloglf [] []};
[phat, ~, ~, err, output] = dfswitchyard('statsfminbx',...
                     funfcn, start, [tolBnd; tolBnd], [Inf; Inf], ...
                     options, dfltOptions, 1, x, cens, freq);
if (err == 0)
    % statsfminbx may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if output.funcCount >= options.MaxFunEvals
        warning(message('stats:bisafit:EvalLimit'));
    else
        warning(message('stats:bisafit:IterLimit'));
    end
elseif (err < 0)
    error(message('stats:bisafit:NoSolution'));
end

% Compute CIs using a normal approximation for phat.
if nargout > 1
    acov = mlecov(phat, x, 'nloglf',@bisa_nloglf, 'cens',cens, 'freq',freq);
    pci = dfswitchyard('statparamci',phat,acov,alpha,true(1,3));
end
end

function [nll,ngrad] = bisa_nloglf(params, data, cens, freq)
%BISA_NLOGLF Objective function for Birnbaum-Saunders maximum likelihood.

beta = params(1);
gamma = params(2);
z = (sqrt(data./beta) - sqrt(beta./data)) ./ gamma;
w = (sqrt(data./beta) + sqrt(beta./data)) ./ gamma;

logphi = -0.5 .* (z.^2 + log(2.*pi));
L = logphi + log(w) - log(2.*data);
ncen = sum(freq.*cens);
if ncen > 0
    cen = (cens == 1);
    zcen = z(cen);
    Scen = 0.5 * erfc(zcen ./ sqrt(2));
    L(cen) = log(Scen);
end
nll = -sum(freq .* L);

if nargout > 1
    dL1 = (w.^2 - 1) .* 0.5.*z./(w.*beta);
    dL2 = (z.^2 - 1) ./ gamma;
    if ncen > 0
        phicen = exp(logphi(cen));
        wcen = w(cen);
        d1Scen = phicen .* 0.5.*wcen./beta;
        d2Scen = phicen .* zcen./gamma;
        dL1(cen) = d1Scen ./ Scen;
        dL2(cen) = d2Scen ./ Scen;
    end
    ngrad = -[sum(freq .* dL1) sum(freq .* dL2)];
end
end
