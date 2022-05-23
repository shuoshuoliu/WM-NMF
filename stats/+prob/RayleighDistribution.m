classdef RayleighDistribution < prob.ToolboxFittableParametricDistribution
%RayleighDistribution Rayleigh probability distribution.
%    An object of the RayleighDistribution class represents a Rayleigh
%    probability distribution with specific parameter B.
%    This distribution object can be created directly using the MAKEDIST
%    function or fit to data using the FITDIST function.
%
%    RayleighDistribution methods:
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
%    RayleighDistribution properties:    
%       DistributionName      - Name of the distribution
%       B                     - Value of the B parameter
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
%B Value of B parameter
%    The B property represents the defining parameter of the Rayleigh 
%    distribution
%
        B 
    end
    properties(Hidden,Dependent=true)
       b
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameRayleigh'));

%NumParameter Number of parameters.
%    NumParameters is the number of parameters in the distribution.
%
%    See also ParameterValues.
        NumParameters = 1;

%ParameterNames Parameter names.
%    ParameterNames is a cell array of strings containing the names of the
%    parameters of the probability distribution.
%
%    See also ParameterValues, ParameterDescription.
        ParameterNames = {'B'};

%ParameterDescription Parameter description.
%    ParameterNames is a cell array of strings containing short
%    descriptions of the parameters of the probability distribution.
%
%    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionScale'))};
    end
    properties(GetAccess='public',SetAccess='protected')
%ParameterValues Parameter values.
%    ParameterVales is a vector containing the values of the parameters of
%    the probability distribution.
%
%    See also B.
        ParameterValues
    end
    methods(Hidden)
        function pd = RayleighDistribution(b)
            if nargin==0
                b = 1;
            end
            checkargs(b)

            pd.ParameterValues = b;
            pd.ParameterIsFixed = true;
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
            m = raylstat(this.B);
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            [~,v] = raylstat(this.B);
        end
       function ci = paramci(this,varargin)
           [varargin{:}] = convertStringsToChars(varargin{:});
           
           requireScalar(this)
           if isscalar(varargin) && isnumeric(varargin{1})
               % Support syntax of older ProbDistUnivParam/paramci method
               varargin = {'alpha' varargin{1}};
           end
           okargs =   {'alpha' 'type' 'parameter' 'logflag'};
           defaults = {0.05, 'exact', 1:this.NumParameters false};
           [alpha,citype] = internal.stats.parseArgs(okargs,defaults,varargin{:});
           if ~(isscalar(alpha) && isnumeric(alpha) && alpha>0 && alpha<1)
               error(message('stats:ecdf:BadAlpha'))
           end
           if isequal(citype,'exact')  && ~all(this.ParameterIsFixed)
               freq = this.InputData.freq;
               % Number of observations
               if isempty(freq) || isequal(freq,1)
                   N = length(this.InputData.data);
               else
                   N = sum(freq);
               end
               % The exact confidence interval is based on chi-square
               p_int = [1-alpha/2; alpha/2];
               ci = sqrt(2 * N * (this.B)^2 ./ chi2inv(p_int, 2*N));
           else
               ci = paramci@prob.ToolboxFittableParametricDistribution(this,varargin{:});
           end
       end
    end
    methods
        function this = set.B(this,b)
            checkargs(b);
            this.ParameterValues(1) = b;
            this = invalidateFit(this);
        end
        function b = get.B(this)
            b = this.ParameterValues(1);
        end
        function this = set.b(this,b)
            checkargs(b);
            this.ParameterValues(1) = b;
            this = invalidateFit(this);
        end
        function b = get.b(this)
            b = this.ParameterValues(1);
        end
    end
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the rayleigh distribution to data.
%    Fitting requires the Statistics and Machine Learning Toolbox. You should call the FITDIST
%    function instead of calling this method directly.
%
%    See also FITDIST.

            [xOriginal,cens,freq] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            x = prob.ToolboxFittableParametricDistribution.removeCensoring(xOriginal,cens,freq,'rayleigh');
            params = raylfit(x,0.05);
            [nll,cov] = prob.RayleighDistribution.likefunc(params,xOriginal,cens,freq);
            pd = prob.RayleighDistribution.makeFitted(params,nll,cov,xOriginal,cens,freq);
        end
        function varargout = likefunc(varargin)
           % Supports [NLOGL,AVAR] = likefunc(PARAMS,DATA,CENSORING,FREQUENCY)
           [varargout{1:nargout}] = rayllike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = raylcdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = raylpdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = raylinv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = raylrnd(varargin{:});
        end
        function pd = makeFitted(p,nll,cov,x,cens,freq)
            pd = prob.RayleighDistribution(p);
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = false;
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.RayleighDistribution');
            info.name = prob.RayleighDistribution.DistributionName;
            info.code = 'rayleigh';
            info.islocscale = true;
            info.logci = false;
            info.support = [0 Inf];
        end
    end
end % classdef

function checkargs(b)
if ~(isscalar(b) && isnumeric(b) && isreal(b) && b>0 && isfinite(b))
    error(message('stats:probdists:PositiveParameter','B'))
end
end

function [nlogL,acov] = rayllike(b,data,cens,freq)
% Analytical negative log-likelihood and covariance for the Rayleigh
% Distribution.  ACOV = (2nd derivative(log-likelihood)^-1.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end
nlogL = rayl_nloglf(b,data,cens,freq);
if nargout>1
    z2 = (data/b).^2;
    % log''(f) if uncensored
    d2 = (2 - 3 * z2) / b^2;
    censored_obs = (cens==1);
    % log''(S) if censored, S is the survival function
    d2(censored_obs) = - 3 * z2(censored_obs) / b^2;
    acov = - 1 / dot(freq,d2);
end
end

function [nll,ngrad] = rayl_nloglf(b, x, cens, freq)
% Analytical negative log-likelihood and its derivative for the Rayleigh
% Distribution.
z2 = (x/b).^2;
% log(f) and log'(f) if uncensored
logfz = -2*log(b) - z2/2 + log(x);
dlogfz = (z2 - 2) / b;
% log(S) and log'(S) if censored, S is survival function.
logS = - z2/2;
dlogS = - z2 * 3 / b^2;

L = logfz;
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end
censored_obs = (cens==1);
L(censored_obs) = logS(censored_obs);
nll = - dot(freq,L);

if nargout>1
    G = dlogfz;
    G(censored_obs) = dlogS(censored_obs);
    ngrad = - dot(freq,G);
end
end
