classdef GammaDistribution < prob.ToolboxFittableParametricDistribution
%GammaDistribution Gamma probability distribution.
%    An object of the GammaDistribution class represents an gamma
%    probability distribution with specific values of the parameters A and B.
%    This distribution object can be created directly 
%    using the MAKEDIST function or fit to data using the FITDIST function.
%
%    GammaDistribution methods:
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
%    GammaDistribution properties:    
%       DistributionName      - Name of the distribution
%       a                     - Value of the a parameter
%       b                     - Value of the b parameter
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

%    Copyright 2012 The MathWorks, Inc.

    properties(Dependent=true)
%A  The A property represents the shape parameter A of the gamma distribution.
%
%    See also ParameterValues.
        a
%B  The B property represents the scale parameter B of the gamma distribution.
%
%    See also ParameterValues.
        b
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameGamma'));

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
        ParameterNames = {'a' 'b'};

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
%    See also MU.
        ParameterValues
    end
    
    methods(Hidden)
        function pd = GammaDistribution(a,b)
            if nargin==0
                a = 1;
                b = 1;
            end
            checkargs(a,b)
          
            pd.ParameterValues = [a b];
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
            m = gamstat(this.a,this.b);
        end

        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            [~,v] = gamstat(this.a,this.b);
        end
        function this = set.a(this,a)
            checkargs(a,this.b)
            this.ParameterValues(1) = a;
            this = invalidateFit(this);
        end
        function a = get.a(this)
            a = this.ParameterValues(1);
        end
        function this = set.b(this,b)
            checkargs(this.a,b)
            this.ParameterValues(2) = b;
            this = invalidateFit(this);
        end
        function b = get.b(this)
            b = this.ParameterValues(2);
        end
    end
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the beta distribution to data. You
%    should call the FITDIST function instead of calling this method
%    directly.
%
%    See also FITDIST.
            [x,cens,freq,opt] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            p = gamfit(x,0.05,cens,freq,opt);
            [nll,cov] = gamlike(p,x,cens,freq);
            pd = prob.GammaDistribution.makeFitted(p,nll,cov,x,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = gamlike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = gamcdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = gampdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = gaminv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = gamrnd(varargin{:});
        end
        function pd = makeFitted(p,nll,cov,x,cens,freq)
            pd = prob.GammaDistribution(p(1),p(2));
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.GammaDistribution');
            info.name = prob.GammaDistribution.DistributionName;
            info.code = 'gamma';
            info.hasconfbounds = true;
            info.censoring = true;
            info.support = [0 Inf];
            info.optimopts = true;
            info.logci = [true true];
        end
    end
end % classdef

function checkargs(a,b)
if ~(isscalar(a) && isnumeric(a) && isreal(a) && ~isnan(a) && a>0)
    error(message('stats:probdists:PositiveParameter','A'))
end
if ~(isscalar(b) && isnumeric(b) && isreal(b) && isfinite(b) && b>=0)
    error(message('stats:probdists:PositiveParameter','B'))
end
end
