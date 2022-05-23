classdef NakagamiDistribution < prob.ToolboxFittableParametricDistribution
%NakagamiDistribution Nakagami probability distribution.
%    An object of the NakagamiDistribution class represents an Nakagami
%    probability distribution with specific values of the parameters MU and
%    OMEGA. This distribution object can be created directly using the
%    MAKEDIST function or fit to data using the FITDIST function.
%
%    NakagamiDistribution methods:
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
%    NakagamiDistribution properties:    
%       DistributionName      - Name of the distribution
%       mu                    - Value of the mu parameter
%       omega                 - Value of the omega parameter
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

    properties(Dependent=true)
%MU  The MU property represents the shape parameter MU of the Nakagami
%distribution.
%
%    See also ParameterValues.
        mu
%OMEGA  The OMEGA property represents the scale parameter OMEGA of the
%Nakagami distribution.
%
%    See also ParameterValues.
        omega
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameNakagami'));

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
        ParameterNames = {'mu' 'omega'};

%ParameterDescription Parameter description.
%    ParameterNames is a cell array of strings containing short
%    descriptions of the parameters of the probability distribution.
%
%    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionShape')) ...
                                getString(message('stats:probdists:ParameterDescriptionScale'))};
    end
    
    properties(GetAccess='public',SetAccess='protected')
%ParameterValues Parameter values.
%    ParameterVales is a vector containing the values of the parameters of
%    the probability distribution.
%
%    See also MU, OMEGA.
        ParameterValues
    end
    
    methods(Hidden)
        function pd = NakagamiDistribution(mu,omega)
            if nargin==0
                mu = 1;
                omega = 1;
            end
            checkargs(mu,omega)
          
            pd.ParameterValues = [mu omega];
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
            gamratio = exp(gammaln(this.mu+.5) - gammaln(this.mu));
            m = gamratio .* sqrt(this.omega./this.mu);
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            gamratio = exp(gammaln(this.mu+.5) - gammaln(this.mu));
            v = this.omega .* (1 - gamratio.^2 ./ this.mu);
        end
    end
    methods
        function this = set.mu(this,mu)
            checkargs(mu,this.omega)
            this.ParameterValues(1) = mu;
            this = invalidateFit(this);
        end
        function mu = get.mu(this)
            mu = this.ParameterValues(1);
        end
        function this = set.omega(this,omega)
            checkargs(this.mu,omega)
            this.ParameterValues(2) = omega;
            this = invalidateFit(this);
        end
        function omega = get.omega(this)
            omega = this.ParameterValues(2);
        end
    end
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the Nakagami distribution to data. You
%    should call the FITDIST function instead of calling this method
%    directly.
%
%    See also FITDIST.
            [x,cens,freq,opt] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            p = nakafit(x,0.05,cens,freq,opt);
            [nll,cov] = nakalike(p,x,cens,freq);
            pd = prob.NakagamiDistribution.makeFitted(p,nll,cov,x,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = nakalike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = nakacdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = nakapdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = nakainv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = nakarnd(varargin{:});
        end
        function pd = makeFitted(p,nll,cov,x,cens,freq)
            pd = prob.NakagamiDistribution(p(1),p(2));
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.NakagamiDistribution');
            info.name = prob.NakagamiDistribution.DistributionName;
            info.code = 'nakagami';
            info.censoring = true;
            info.support = [0 Inf];
            info.islocscale = false;
            info.optimopts = true;
            info.logci = [true true];
        end
    end
end % classdef

function checkargs(mu,omega)
if ~(isscalar(mu) && isnumeric(mu) && isreal(mu) && isfinite(mu) && mu>0)
    error(message('stats:probdists:PositiveParameter','MU'))
end
if ~(isscalar(omega) && isnumeric(omega) && isreal(omega) && isfinite(omega) && omega>0)
    error(message('stats:probdists:PositiveParameter','OMEGA'))
end
end

function y = nakapdf(x, mu, omega)
%NAKAPDF Nakagami probability density function (pdf).
mu(mu <= 0) = NaN;
omega(omega <= 0) = NaN;

x(x<0) = 0;
% equivalent to y = 2.*x .* gampdf(x.^2, mu, omega./mu), but the version here
% puts all the x terms into gampdf, so Inf*0, etc. is handled there.
y = 2.*sqrt(omega./mu).*exp(gammaln(mu+.5) - gammaln(mu)) .* gampdf(x.^2, mu+.5, omega./mu);
end

function p = nakacdf(x, mu, omega, uflag)
%NAKACDF Nakagami cumulative distribution function (cdf).
mu(mu <= 0) = NaN;
omega(omega <= 0) = NaN;

x(x<0) = 0;
if nargin > 3
    if ~strcmpi(uflag,'upper')
        error(message('stats:cdf:UpperTailProblem'));
    else
        p = gamcdf(x.^2,mu,omega./mu,'upper');
    end
else
    p = gamcdf(x.^2, mu, omega./mu);
end
end

function x = nakainv(p, mu, omega)
%NAKAINV Inverse of the Nakagami cumulative distribution function (cdf).
mu(mu <= 0) = NaN;
omega(omega <= 0) = NaN;

x = sqrt(gaminv(p,mu,omega./mu));
end

function r = nakarnd(mu, omega, varargin)
%NAKARND Random arrays from the Nakagami distribution.
mu(mu <= 0) = NaN;
omega(omega <= 0) = NaN;

[err, sizeOut] = internal.stats.statsizechk(2,mu,omega,varargin{:});
if err > 0
    error(message('stats:nakarnd:InconsistentSizes'));
end

r = sqrt(gamrnd(mu,omega./mu,sizeOut));
end

function [nlogL,acov] = nakalike(params,data,cens,freq)
%NAKALIKE Negative log-likelihood for the Nakagami distribution.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end

nlogL = naka_nloglf(params, data, cens, freq);
if nargout > 1
    acov = mlecov(params, data, 'nloglf',@naka_nloglf, 'cens',cens, 'freq',freq);
end
end

% ==== Nakagami fitting functions ====

function [phat,pci] = nakafit(x,alpha,cens,freq,opts)
%NAKAFIT Parameter estimates and confidence intervals for Nakagami data.

if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = zeros(size(x)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); end
if nargin < 5, opts = []; end

if any(x <= 0)
    error(message('stats:nakafit:BadData'));
end

phat = gamfit(x.^2,alpha,cens,freq,opts);
phat(2) = phat(1).*phat(2); % (a,b) -> (mu,omega)
if nargout > 1
    acov = mlecov(phat, x, 'nloglf',@naka_nloglf, 'cens',cens, 'freq',freq);
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';
    pci = norminv(repmat(probs,1,numel(phat)), [phat; phat], [se; se]);
    % CI on the log scale for omega?
end
end

function [nll,ngrad] = naka_nloglf(parms, x, cens, freq)
%NAKA_NLOGLF Objective function for Nakagami maximum likelihood.

% do all the calculations in terms of the gamma dist'n
a = parms(1);
b = parms(2)./parms(1); % (mu,omega) -> (a,b)
loggama = gammaln(a);
logb = log(b);

xsq = x.^2;
z = xsq ./ b;
logz = log(z);
L = (a-1).*logz - z - loggama - logb + log(2.*x);
ncen = sum(freq.*cens);
if ncen > 0
    cen = (cens == 1);
    zcen = z(cen);
    if nargout == 1
        Scen = gammainc(zcen,a,'upper');
    else
        [dScen,Scen] = dgammainc(zcen,a,'upper');
    end
    L(cen) = log(Scen);
end
nll = -sum(freq .* L);

if nargout > 1
    dL1 = logz - psi(a);
    dL2 = (z - a)./b;
    if ncen > 0
        dL1(cen) = dScen ./ Scen;
        dL2(cen) = exp(a.*logz(cen) - logb - zcen - loggama) ./ Scen;
    end
    ngrad = -[sum(freq .* dL1) sum(freq .* dL2)];

    % transform back to Nakagami parameters
    ngrad = ngrad * [1 0; -b./a 1./a]; % (a,b) -> (mu,omega)
end
end
