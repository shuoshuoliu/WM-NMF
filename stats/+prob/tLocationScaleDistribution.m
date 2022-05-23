classdef tLocationScaleDistribution < prob.ToolboxFittableParametricDistribution
%tLocationScaleDistribution t location-scale probability distribution.
%    An object of the tLocationScaleDistribution class represents a
%    t location-scale probability distribution with specific values of
%    the parameters MU, SIGMA, and NU. This distribution object can be
%    created directly using the MAKEDIST function or fit to data using
%    the FITDIST function.
%
%    tLocationScaleDistribution methods:
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
%    tLocationScaleDistribution properties:    
%       DistributionName      - Name of the distribution
%       mu                    - Value of the mu parameter (location)
%       sigma                 - Value of the sigma parameter (scale)
%       nu                    - Value of the nu parameter (degrees of freedom)
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
%MU  The MU property represents the location parameter MU of the
%    t location-scale distribution.
%
%    See also ParameterValues.
        mu

%SIGMA  The SIGMA property represents the scale parameter SIGMA of the
%t location-scale distribution.
%
%    See also ParameterValues.
        sigma

%NU  The NU property represents the degrees of freedom of the
%t location-scale distribution.
%
%    See also ParameterValues.
        nu
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NametLocationScale'));
        
%NumParameter Number of parameters.
%    NumParameters is the number of parameters in the distribution.
%
%    See also ParameterValues.
        NumParameters = 3;
        
%ParameterNames Parameter names.
%    ParameterNames is a cell array of strings containing the names of the
%    parameters of the probability distribution.
%
%    See also ParameterValues, ParameterDescription.
        ParameterNames = {'mu' 'sigma' 'nu'};

%ParameterDescription Parameter description.
%    ParameterNames is a cell array of strings containing short
%    descriptions of the parameters of the probability distribution.
%
%    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionLocation')) ...
                                getString(message('stats:probdists:ParameterDescriptionScale')) ...
                                getString(message('stats:probdists:ParameterDescriptionDOF'))};
    end
    
    properties(GetAccess='public',SetAccess='protected')
%ParameterValues Parameter values.
%    ParameterVales is a vector containing the values of the parameters of
%    the probability distribution.
%
%    See also MU, SIGMA, NU.
        ParameterValues
    end
    
    methods(Hidden)
        function pd = tLocationScaleDistribution(mu,sigma,nu)
            if nargin==0
                mu = 0;
                sigma = 1;
                nu = 5;
            end
            checkargs(mu,sigma,nu)
          
            pd.ParameterValues = [mu sigma nu];
            pd.ParameterIsFixed = [true true true];
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
            if this.nu <= 1
                m = NaN;
            else
                m = this.mu;
            end
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            if this.nu <= 2
                v = Inf;
            else
                v = this.sigma.^2 .* this.nu ./ (this.nu - 2);
            end
        end
    end
    methods
        function this = set.mu(this,mu)
            checkargs(mu,this.sigma,this.nu)
            this.ParameterValues(1) = mu;
            this = invalidateFit(this);
        end
        function mu = get.mu(this)
            mu = this.ParameterValues(1);
        end
        function this = set.sigma(this,sigma)
            checkargs(this.mu,sigma,this.nu)
            this.ParameterValues(2) = sigma;
            this = invalidateFit(this);
        end
        function sigma = get.sigma(this)
            sigma = this.ParameterValues(2);
        end
        function this = set.nu(this,nu)
            checkargs(this.mu,this.sigma,nu)
            this.ParameterValues(3) = nu;
            this = invalidateFit(this);
        end
        function nu = get.nu(this)
            nu = this.ParameterValues(3);
        end
    end
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the t location-scale distribution to
%    data. You should call the FITDIST function instead of calling this method
%    directly.
%
%    See also FITDIST.
            [x,cens,freq,opt] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            p = tlsfit(x,0.05,cens,freq,opt);
            [nll,cov] = tlslike(p,x,cens,freq);
            pd = prob.tLocationScaleDistribution.makeFitted(p,nll,cov,x,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = tlslike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = tlscdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = tlspdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = tlsinv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = tlsrnd(varargin{:});
        end
        function pd = makeFitted(p,nll,cov,x,cens,freq)
            pd = prob.tLocationScaleDistribution(p(1),p(2),p(3));
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.tLocationScaleDistribution');
            info.name = prob.tLocationScaleDistribution.DistributionName;
            info.code = 'tlocationscale';
            info.censoring = true;
            info.islocscale = false;
            info.optimopts = true;
            info.logci = [false true true];
        end
    end
end % classdef

function checkargs(mu,sigma,nu)
if ~(isscalar(mu) && isnumeric(mu) && isreal(mu) && isfinite(mu))
    error(message('stats:probdists:ScalarParameter','MU'))
end
if ~(isscalar(sigma) && isnumeric(sigma) && isreal(sigma) && isfinite(sigma) && sigma>0)
    error(message('stats:probdists:PositiveParameter','SIGMA'))
end
if ~(isscalar(nu) && isnumeric(nu) && isreal(nu) && isfinite(nu) && nu>0)
    error(message('stats:probdists:PositiveParameter','NU'))
end
end

% ==== t Location-Scale distribution functions ====

% these distribution functions do not yet handle arrays of parameters

function y = tlspdf(x, mu, sigma, nu)
%TLSPDF T location-scale probability density function (pdf).
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;

y = tpdf((x - mu)./sigma,nu)./sigma;
end

function p = tlscdf(x ,mu, sigma, nu, uflag)
%TLSCDF T location-scale cumulative distribution function (cdf).
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;

if nargin > 4
    if ~strcmpi(uflag,'upper')
        error(message('stats:cdf:UpperTailProblem'));
    else
        p = tcdf((x - mu)./sigma,nu,'upper');
    end
else
    p = tcdf((x - mu)./sigma,nu);
end
end

function x = tlsinv(p, mu, sigma, nu)
%TLSINV Inverse of the t location-scale cumulative distribution function (cdf).
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;

x = tinv(p,nu).*sigma + mu;
end

function r = tlsrnd(mu, sigma, nu, varargin)
%TLSRND Random arrays from the t location-scale distribution.
sigma(sigma <= 0) = NaN;
nu(nu <= 0) = NaN;

[err, sizeOut] = internal.stats.statsizechk(3,mu,sigma,nu,varargin{:});
if err > 0
    error(message('stats:tlsrnd:InconsistentSizes'));
end

r = mu + sigma.*trnd(nu,sizeOut);
end

function [nlogL,acov] = tlslike(params,data,cens,freq)
%TLSLIKE Negative log-likelihood for the t location-scale distribution.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end

nlogL = tls_nloglf(params, data, cens, freq);
if nargout > 1
    acov = mlecov(params, data, 'nloglf',@tls_nloglf, 'cens',cens, 'freq',freq);
end
end

% ==== t location-scale fitting functions ====

function [phat,pci] = tlsfit(x,alpha,cens,freq,opts)

if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = zeros(size(x)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); end
if nargin < 5, opts = []; end

% Robust estimators for the mean and std dev of a normal, and method
% of moments on t-kurtosis for nu
xunc = x(cens == 0);
k = max(kurtosis(xunc), 4);
start = [median(xunc), 1.253.*mad(xunc), 2.*(2.*k-3)./(k-3)];

% The default options include turning fminsearch's display off.  This
% function gives its own warning/error messages, and the caller can turn
% display on to get the text output from fminsearch if desired.
options = statset(statset('tlsfit'), opts);
tolBnd = options.TolBnd;
options = optimset(options);

% Maximize the log-likelihood with respect to mu, sigma, and nu.
[phat,~,err,output] = ...
    fminsearch(@tls_nloglf, start, options, x, cens, freq, tolBnd);
if (err == 0)
    % fminsearch may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if phat(3) > 100 % degrees of freedom became very large
        wmsg = getString(message('stats:tlsfit:SuggestNormal'));
    else
        wmsg = '';
    end
    if output.funcCount >= options.MaxFunEvals
        warning(message('stats:tlsfit:EvalLimit', wmsg));
    else
        warning(message('stats:tlsfit:IterLimit', wmsg));
    end
elseif (err < 0)
    error(message('stats:tlsfit:NoSolution'));
end

if nargout > 1
    acov = mlecov(phat, x, 'nloglf',@tls_nloglf, 'cens',cens, 'freq',freq);
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';

    % Compute the CI for mu using a normal approximation for muhat.
    pci(:,1) = norminv(probs, phat(1), se(1));

    % Compute the CI for sigma using a normal approximation for
    % log(sigmahat), and transform back to the original scale.
    % se(log(sigmahat)) is se(sigmahat) / sigmahat.
    logsigci = norminv(probs, log(phat(2)), se(2)./phat(2));
    pci(:,2) = exp(logsigci);

    % Compute the CI for nu using a normal distribution for nuhat.
    pci(:,3) = norminv(probs, phat(3), se(3));
end
end

function nll = tls_nloglf(parms, x, cens, freq, tolBnd)
%TLS_NLOGLF Objective function for t location-scale maximum likelihood.
mu = parms(1);
sigma = parms(2);
nu = parms(3);

% Restrict sigma and nu to the open interval (0, Inf).
if nargin > 4
    if sigma < tolBnd || nu < tolBnd
        nll = Inf;
        return
    end
end

t = (x - mu) ./ sigma;
w = nu + (t.^2);
logw = log(w);

L = -.5.*(nu+1).*logw + gammaln(.5.*(nu+1)) - gammaln(.5.*nu) + 0.5.*nu.*log(nu) - log(sigma) - .5.*log(pi);
ncen = sum(freq.*cens);
if ncen > 0
    cen = (cens == 1);
    if nu < 1e7  % Use the standard formula for F(-|t|) = 1-F(|t|) < .5
        Scen = betainc(nu ./ w(cen), .5.*nu, 0.5) ./ 2;

        % Reflect for negative t.
        reflect = (t(cen) < 0);
        Scen(reflect) = 1 - Scen(reflect); % Scen < .5, cancellation not a problem

    else  % Use a normal approximation.
        Scen = log(0.5 * erfc(t(cen) ./ sqrt(2)));
    end
    L(cen) = log(Scen);
end
nll = -sum(freq .* L);

% Don't yet have dbetainc, so can't compute an analytic gradient with censoring.
%
% if nargout > 1
%     dL1 = (nu+1).*t./(w.*sigma);
%     dL2 = t.*dL1 - 1./sigma;
%     dL3 = .5.*(-logw - (nu+1)./w + psi(.5.*(nu+1)) - psi(.5.*nu) + log(nu) + 1);
%     if ncen > 0
% %         dL1(cen) = ;
% %         dL2(cen) = ;
% %         dL3(cen) = ;
%     end
%     ngrad = -[sum(freq .* dL1) sum(freq .* dL2) sum(freq .* dL3)];
% end
end
