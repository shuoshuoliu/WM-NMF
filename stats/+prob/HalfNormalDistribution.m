classdef HalfNormalDistribution < prob.ToolboxFittableParametricDistribution
%HalfNormalDistribution Half-normal probability distribution.
%    An object of the HalfNormalDistribution class represents a half-normal
%    probability distribution with specific values of the parameters
%    MU and SIGMA. This distribution object can be created directly 
%    using the MAKEDIST function or fit to data using the FITDIST function.
%
%    HalfNormalDistribution methods:
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
%    HalfNormalDistribution properties:    
%       DistributionName      - Name of the distribution
%       mu                    - Value of the mu parameter
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

%    Copyright 2015-2018 The MathWorks, Inc.

    properties(Dependent=true)
%MU Value of MU parameter
%    The MU property represents the location parameter of the half-normal 
%    distribution. It is also the lower limit of the distribution. If Z has
%    a normal distribution with mean MU and standard deviation SIGMA, then 
%    X = MU+abs(Z-MU) has a half-normal distribution with parameters MU and SIGMA. 
%
%    See also SIGMA.
        mu
        
%SIGMA Value of SIGMA parameter
%    The SIGMA property represents the scale parameter of the half-normal 
%    distribution. If Z has a normal distribution with mean MU and standard
%    deviation SIGMA, then X = MU+abs(Z-MU) has a half-normal distribution 
%    with parameters MU and SIGMA.
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
        DistributionName = getString(message('stats:dfittool:NameHalfNormal'));

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
                                getString(message('stats:probdists:ParameterDescriptionScale')) ...
};
    end
    properties(GetAccess='public',SetAccess='protected')
%ParameterValues Parameter values.
%    ParameterValues is a vector containing the values of the parameters of
%    the probability distribution.
%
%    See also MU, SIGMA.
        ParameterValues
    end
    methods(Hidden)
        function pd = HalfNormalDistribution(mu,sigma)
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
            m = this.mu + sqrt(2/pi)*this.sigma;
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            s2 = this.sigma .^ 2;
            v = (1-2/pi)*s2;
        end
        function ci = paramci(this,varargin)
            [varargin{:}] = convertStringsToChars(varargin{:});
            
            requireScalar(this)
            if isscalar(varargin) && isnumeric(varargin{1})
                % Support syntax of older ProbDistUnivParam/paramci method
                varargin = {'alpha' varargin{1}};
            end
            okargs =   {'alpha' 'parameter' 'type' 'logflag'};
            defaults = {0.05, 1:this.NumParameters 'exact' []};
            [alpha,~,citype] = internal.stats.parseArgs(okargs,defaults,varargin{:});
            if ~(isscalar(alpha) && isnumeric(alpha) && alpha>0 && alpha<1)
                error(message('stats:ecdf:BadAlpha'))
            end
            if isequal(citype,'exact')  && ~all(this.ParameterIsFixed)
                ci = stathnci([this.mu this.sigma], alpha, this.InputData.data);
            else
                ci = paramci@prob.ToolboxFittableParametricDistribution(this,varargin{:});
            end
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
%    FIT is a static method that fits the half-normal distribution to data.
%    You should call the FITDIST function instead of calling this method directly.
%
%    See also FITDIST.

%    PD = PROB.HALFNORMALDISTRIBUTION.FIT(X) creates a 
%    HalfNormalDistribution object with MU=0 and with parameter
%    SIGMA fit by maximum likelihood to the data in X.
%
%    PD = PROB.HALFNORMALDISTRIBUTION.FIT(X,NAME,VALUE) specifies
%    the following optional parameter:
%
%       'MU'           A scalar numeric value giving the location
%                      parameter of the half-normal distribution.
%                      MU is not estimated from the data.

            [mu,~,fitargs] = internal.stats.parseArgs({'mu'},{0},varargin{2:end});
            [xOriginal,cens,freq,~] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{1}, fitargs{:});
            x = prob.ToolboxFittableParametricDistribution.removeCensoring(xOriginal,cens,freq,'half normal');
            params = hnfit(x,mu,0.05);
            [nll,cov] = hnlike(params,x);
            pd = prob.HalfNormalDistribution.makeFitted(params,nll,cov,xOriginal,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = hnlike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = hncdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = hnpdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = hninv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = hnrnd(varargin{:});
        end
        function pd = makeFitted(params,nll,cov,x,cens,freq)
            c = num2cell(params);
            pd = prob.HalfNormalDistribution(c{:});
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [true false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.HalfNormalDistribution');
            info.name = prob.HalfNormalDistribution.DistributionName;
            info.code = 'half normal';
            info.prequired = [true false];
            info.islocscale = true;
            info.optimopts = false;
            info.logci = [false false];
            info.plim = [-Inf 0 ;Inf Inf];
        end
    end
end % classdef

% The following utilities check for valid parameter values
function checkargs(mu,sigma)
if ~(isscalar(mu) && isnumeric(mu) && isreal(mu) && isfinite(mu))
    error(message('stats:probdists:ScalarParameter','MU'))
end
if ~(isscalar(sigma) && isnumeric(sigma) && isreal(sigma) && sigma>=0 && isfinite(sigma))
    error(message('stats:probdists:PositiveParameter','SIGMA'));
end
end

function [phat,pci] = hnfit(x,mu,alpha)
%HNFIT Parameter estimates and confidence intervals for half-normal value data.
%   HNFIT does not estimate MU, and it must be assumed known.
%
%   References:
%      [1] M. Ahsanullah, B.M. G. Kibria, M. Shakil (2014) Normal and 
%          Student's t Distributions and Their Applications, Atlantis 
%          Press.
%      [2] A. Pewsey (2002) Large-Sample Inference for the General
%          Half-Normal Distribution.

% Illegal data return an error.
if any(x<mu)
    error(message('stats:addhn:BadDataLocation'));
end

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
end

n = size(x,1);

x = x-mu;
sigmahat = sqrt(sum(x.*x)./n);

if nargout < 2
    phat = [mu sigmahat];
else
    phat = [mu sigmahat];
    pci = stathnci(phat,alpha,x);
end
end


function parmci=stathnci(parmhat,alpha,x)
%STATSHNCI Half-normal parameter confidence interval.

% Number of observations
n = size(x,1);

mu = parmhat(1);
sigma = parmhat(2);

% Confidence interval
if n > 0
    chi2crit = chi2inv([alpha/2 1-alpha/2],n);
    parmci = [sigma*sqrt(n./chi2crit(2)); ...
            sigma*sqrt(n./chi2crit(1))];
    parmci = [[mu;mu] parmci];
else
    parmci = NaN(2,numel(parmhat),'like',x);
end
end


function [nlogL,acov] = hnlike(params,data)
%HNLIKE Negative log-likelihood for the half-normal distribution.

if nargin < 2
    error(message('stats:addhn:TooFewInputs'));
elseif numel(data) > length(data)
    error(message('stats:addhn:VectorRequired'));
end

mu = params(1);
sigma = params(2);

% Return NaN for out of range parameter or data.
sigma(sigma <= 0) = NaN;
data(data<mu) = NaN;
z = (data-mu) ./ sigma;

% Sum up the individual log-likelihood terms, and return the negative
% log-likelihood.
logL = -.5.*z.*z - log(sqrt(pi./2).*sigma);
nlogL = -sum(logL);

if nargout == 2
    nH = -sum(1 - 3.*z.*z);
    avar =  (sigma.^2) ./ nH;
    acov = [[0 0]; [0 avar]];
end
end


function [varargout] = hncdf(x,varargin)
%HNCDF Half-normal cumulative distribution function (cdf).

if nargin<1
   error(message('stats:addhn:TooFewInputsX'));
end
if nargin>1 && strcmpi(varargin{end},'upper')
    %Compute upper tail
    uflag=true;
    varargin(end)= [];
elseif nargin>1 && ischar(varargin{end})&& ~strcmpi(varargin{end},'upper')
    error(message('stats:cdf:UpperTailProblem'));
else
    uflag=false;
end
[varargout{1:max(1,nargout)}] = localhncdf(uflag,x,varargin{:});
end

function p = localhncdf(uflag,x,mu,sigma)
if nargin < 3
    mu = 0;
end
if nargin < 4
    sigma = 1;
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    z = (x-mu) ./ sigma;
catch ME
    error(message('stats:addhn:InputSizeMismatch'));
end

% The support is 0 <= (x-mu)/sigma, force zero below that.
z(z<0) = 0;

if uflag == true
    p = erfc(z./sqrt(2));
else
    p = erf(z./sqrt(2));
end
end

function y = hnpdf(x,mu,sigma)
%HNPDF Half-normal probability density function (pdf).

if nargin<1
    error(message('stats:addhn:TooFewInputsX'));
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    z = (x-mu) ./ sigma;
catch ME
    error(message('stats:addhn:InputSizeMismatch'));
end

y = sqrt(2/pi)./sigma.*exp(-0.5 * z.^2);

% The support is 0 <= (x-mu)/sigma, force zero below that.
y(z<0) = 0;
end

function x = hninv(p,mu,sigma)
%HNINV Inverse of the half-normal cumulative distribution function (cdf).

if nargin<1
    error(message('stats:addhn:TooFewInputsP'));
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

% Return NaN for out of range parameters or probabilities.
sigma(sigma <= 0) = NaN;
p(p < 0 | 1 < p) = NaN;

x0 = erfinv(p);
try
    x = sqrt(2)*sigma.*x0 + mu;
catch
    error(message('stats:addhn:InputSizeMismatch'));
end
end

function r = hnrnd(mu,sigma,varargin)
%HNRND Random arrays from the half-normal distribution.

if nargin < 2
    error(message('stats:addhn:TooFewInputs'));
end

[err, sizeOut] = internal.stats.statsizechk(2, mu, sigma, varargin{:});
if err > 0
    error(message('stats:addhn:InputSizeMismatch'));
end

% Return NaN for elements corresponding to illegal parameter values.
sigma(sigma < 0) = NaN;

ty = internal.stats.dominantType(mu, sigma);
r = abs(randn(sizeOut,'like',ty)) .* sigma + mu;
end

