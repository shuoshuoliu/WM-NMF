classdef PoissonDistribution < prob.ToolboxFittableParametricDistribution
%PoissonDistribution Poisson probability distribution.
%    An object of the PoissonDistribution class represents a Poisson
%    probability distribution with a specific mean LAMBDA. This distribution
%    object can be created directly using the MAKEDIST function or fit to
%    data using the FITDIST function.
%
%    PoissonDistribution methods:
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
%    PoissonDistribution properties:    
%       DistributionName      - Name of the distribution
%       lambda                - Value of the lambda parameter (mean)
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
%LAMBDA Mean of Poisson distribution
%    The LAMBDA property represents the parameter that is the mean of the
%    Poisson distribution.
%
%    See also ParameterNames.
        lambda;
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NamePoisson'));

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
        ParameterNames = {'lambda'};

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
%    See also LAMBDA.
        ParameterValues
    end
    
    methods(Hidden)
        function pd = PoissonDistribution(lambda)
            if nargin==0
                lambda = 1;
            end
            checkargs(lambda)
            
            pd.ParameterValues = lambda;
            pd.ParameterIsFixed = true;
            pd.ParameterCovariance = zeros(pd.NumParameters);
        end
    end
    methods
        function m = mean(this)
            requireScalar(this)
            m = this.lambda;
            if this.IsTruncated
                m = truncatedMeanDiscrete(this,m);
            end
        end
        function v = var(this)
            requireScalar(this)
            v = this.lambda;
            if this.IsTruncated
                v = truncatedVarDiscrete(this,v,v);
            end
        end
        function ci = paramci(this,varargin)
            [varargin{:}] = convertStringsToChars(varargin{:});
            
            requireScalar(this)
            if isscalar(varargin) && isnumeric(varargin{1}) && isnumeric(varargin{1})
                % Support syntax of older ProbDistUnivParam/paramci method
                varargin = {'alpha' varargin{1}};
            end
            okargs =   {'alpha' 'parameter' 'type' 'logflag'};
            defaults = {0.05, 1:this.NumParameters 'exact' []};
            [alpha,pnums,citype] = internal.stats.parseArgs(okargs,defaults,varargin{:});
            if ~(isscalar(alpha) && isnumeric(alpha) && alpha>0 && alpha<1)
                error(message('stats:ecdf:BadAlpha'))
            end
            if ~isempty(citype)
                citype = internal.stats.getParamVal(citype,{'wald' 'exact' 'lr'},'''type''');
            end
            if isequal(citype,'exact')
                if isempty(this.InputData.freq)
                    m = numel(this.InputData.data);
                else
                    m = sum(this.InputData.freq);
                end
                ci = dfswitchyard('statpoisci',m,this.lambda,alpha);
                ci = ci(:,pnums);
            else
                ci = paramci@prob.ToolboxFittableParametricDistribution(this,varargin{:});
            end
        end
    end
    methods
        function this = set.lambda(this,lambda)
            checkargs(lambda)
            this.ParameterValues(1) = lambda;
            this = invalidateFit(this);
        end
        function lambda = get.lambda(this)
            lambda = this.ParameterValues(1);
        end
    end
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the Poisson distribution to data. You
%    should call the FITDIST function instead of calling this method directly.
%
%    See also FITDIST.
            [x,cens,freq] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            if any(cens)
                error(message('stats:ProbDistUnivParam:fit:CensoringNotAllowed', 'poisson'));
            elseif any(x<0)
                error(message('stats:poissfit:InvalidX'))
            end
            if isempty(freq)
                freq = ones(size(x));
            end
            n = sum(freq);
            lambda = sum(freq.*x)/n;
            nll = - sum(freq.*log(poisspdf(x,lambda)));
            cov = lambda/n;
            pd = prob.PoissonDistribution.makeFitted(lambda,nll,cov,x,cens,freq);
        end
        function [nll,acov] = likefunc(lambda,data)
            n = length(data);
            nll = -sum(log(poisspdf(data,lambda)));
            acov = lambda/n;
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = poisscdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = poisspdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = poissinv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = poissrnd(varargin{:});
        end
        function pd = makeFitted(lam,nll,cov,x,cens,freq)
            pd = prob.PoissonDistribution(lam);
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = false;
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.PoissonDistribution');
            info.name = getString(message('stats:dfittool:NamePoisson'));
            info.code = 'poisson';
            info.support = [0 Inf];
            info.closedbound = [true false];
            info.iscontinuous = false;
            info.plim = [0; Inf];
        end
    end
end % classdef

function checkargs(lam)
if ~(isscalar(lam) && isnumeric(lam) && isreal(lam) && lam>=0 && isfinite(lam))
    error(message('stats:probdists:NonnegativeParameter','LAMBDA'))
end
end
