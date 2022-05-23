classdef ExponentialDistribution < prob.ToolboxFittableParametricDistribution
%ExponentialDistribution Exponential probability distribution.
%    An object of the ExponentialDistribution class represents an exponential
%    probability distribution with specific values of the parameter MU.
%    This distribution object can be created directly 
%    using the MAKEDIST function or fit to data using the FITDIST function.
%
%    ExponentialDistribution methods:
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
%    ExponentialDistribution properties:    
%       DistributionName      - Name of the distribution
%       mu                    - Value of the mu parameter (mean)
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
%MU Defining parameter of exponential distribution
%    The MU property represents the defining parameter, and the mean, 
%    of the exponential distribution
%
%    See also ParameterValues.
        mu
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameExponential'));

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
        ParameterNames = {'mu'};

%ParameterDescription Parameter description.
%    ParameterNames is a cell array of strings containing short
%    descriptions of the parameters of the probability distribution.
%
%    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionMean'))};
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
        function pd = ExponentialDistribution(mu)
            if nargin==0
                mu = 1;
            end
            checkargs(mu)
          
            pd.ParameterValues = mu;
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
            m = this.mu;
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            v = this.mu^2;
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
                ci = dfswitchyard('statexpci', this.mu, this.ParameterCovariance, alpha, ...
                    this.InputData.data, this.InputData.cens, this.InputData.freq);
            else
                ci = paramci@prob.ToolboxFittableParametricDistribution(this,varargin{:});
            end
        end
        function this = set.mu(this,mu)
            checkargs(mu)
            this.ParameterValues(1) = mu;
            this = invalidateFit(this);
        end
        function mu = get.mu(this)
            mu = this.ParameterValues(1);
        end
    end
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the extreme value distribution to
%    data. You should call the FITDIST function instead of calling this
%    method directly.
%
%    See also FITDIST.
            [x,cens,freq,~] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            p = expfit(x,0.05,cens,freq);
            [nll,cov] = explike(p,x,cens,freq);
            pd = prob.ExponentialDistribution.makeFitted(p,nll,cov,x,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = explike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = expcdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = exppdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = expinv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = exprnd(varargin{:});
        end
        function pd = makeFitted(mu,nll,cov,x,cens,freq)
            pd = prob.ExponentialDistribution(mu);
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = false;
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.ExponentialDistribution');
            info.name = prob.ExponentialDistribution.DistributionName;
            info.code = 'exponential';
            info.hasconfbounds = true;
            info.censoring = true;
            info.support = [0 Inf];
            info.closedbound = [true false];
            info.islocscale = true;
            info.optimopts = false;
            info.logci = true;
        end
    
        function name = matlabCodegenRedirect(~) % redirect to the codegen class
             name = 'prob.coder.ExponentialDistribution';
         end    
    
    end
end % classdef

function checkargs(mu)
if ~(isscalar(mu) && isnumeric(mu) && isreal(mu) && isfinite(mu) && mu>0)
    error(message('stats:probdists:NonnegativeParameter','MU'))
end
end
