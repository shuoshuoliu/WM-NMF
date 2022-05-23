classdef BinomialDistribution < prob.ToolboxFittableParametricDistribution
%BinomialDistribution Binomial probability distribution.
%    An object of the BinomialDistribution class represents a binomial
%    probability distribution with a specific number of trials N and
%    success probability P. This distribution object can be created directly
%    using the MAKEDIST function or fit to data using the FITDIST function.
%
%    BinomialDistribution methods:
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
%    BinomialDistribution properties:    
%       DistributionName      - Name of the distribution
%       N                     - Value of the N parameter (number of trials)
%       p                     - Value of the p parameter (success probability)
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
%N Number of trials
%    The N property represents the parameter that is the number of trials
%    for a binomial distribution.
%
%    See also P.
        N

%P Success probability
%    The P property represents the parameter that is the success probability
%    for a binomial distribution.
%
%    See also N.
        p
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameBinomial'));

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
        ParameterNames = {'N' 'p'};

%ParameterDescription Parameter description.
%    ParameterNames is a cell array of strings containing short
%    descriptions of the parameters of the probability distribution.
%
%    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionNumTrials')) ...
                                getString(message('stats:probdists:ParameterDescriptionSuccessProbability'))};
    end
    
    properties(GetAccess='public',SetAccess='protected')
%ParameterValues Parameter values.
%    ParameterVales is a vector containing the values of the parameters of
%    the probability distribution.
%
%    See also N,P.
        ParameterValues
    end
    
    methods(Hidden)
        function pd = BinomialDistribution(N,p)
            if nargin==0
                N = 1;
                p = 0.5;
            end
            checkargs(N,p,0)
            
            pd.ParameterValues = [N p];
            pd.ParameterIsFixed = [true true];
            pd.ParameterCovariance = zeros(pd.NumParameters);
        end
    end
    methods
        function m = mean(this)
            requireScalar(this)
            m = this.N*this.p;
            if this.IsTruncated
                m = truncatedMeanDiscrete(this,m);
            end
        end
        function v = var(this)
            requireScalar(this)
            v = this.N * this.p * (1-this.p);
            if this.IsTruncated
               v = truncatedVarDiscrete(this, this.N*this.p, v);
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
            if isequal(citype,'exact')  && ~all(this.ParameterIsFixed)
                pnums = pNamesToNums(this,pnums);
                id = this.InputData;
                freq = id.freq;
                if isempty(freq)
                    freq = ones(size(id.data));
                end
                ci = dfswitchyard('statbinoci',sum(freq.*id.data),this.N*sum(freq),alpha);
                ci = [[this.N;this.N],ci(:)];
                ci = ci(:,pnums);
            else
                ci = paramci@prob.ToolboxFittableParametricDistribution(this,varargin{:});
            end
        end
        function this = set.N(this,N)
            checkargs(N,this.p,0);
            this.ParameterValues(1) = N;
            this = invalidateFit(this);
        end
        function this = set.p(this,p)
            checkargs(this.N,p,0);
            this.ParameterValues(2) = p;
            this = invalidateFit(this);
        end
        function N = get.N(this)
            N = this.ParameterValues(1);
        end
        function p = get.p(this)
            p = this.ParameterValues(2);
        end
    end
    methods(Static,Hidden)
        function pd = fit(x,varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the binomial distribution to data.
%    You should call the FITDIST function instead of calling this method
%    directly.
%
%    See also FITDIST.
            [N,~,varargin] = internal.stats.parseArgs({'ntrials'},{1},varargin{:});
            [x,cens,freq] = prob.ToolboxFittableParametricDistribution.processFitArgs(x,varargin{:});
            if any(cens)
                 error(message('stats:ProbDistUnivParam:fit:CensoringNotAllowed', 'binomial'));
            end
            checkargs(N,.5,x);
            if isempty(freq)
                freq = ones(size(x));
            end
            n = N*sum(freq);
            p = sum(freq.*x)/n;
            nll = -sum(freq.*log(binopdf(x,N,p)));
            cov = [0 0;0 p*(1-p)/n];
            pd = prob.BinomialDistribution.makeFitted([N,p],nll,cov,x,cens,freq);
        end
        function [nll,acov] = likefunc(params,data)
            N = params(1);
            p = params(2);
            nx = length(data);
            nll = -sum(log(binopdf(data,N,p)));
            acov = [0 0;0 p*(1-p)/(N*nx)];
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = binocdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = binopdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = binoinv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = binornd(varargin{:});
        end
    end
    methods(Hidden,Access=protected)
        function s = fixSupport(~,s)
            s.closedbound = [true true]; % x may be at either extreme
        end
    end
    methods(Static,Hidden)
        function pd = makeFitted(Np,nll,cov,x,cens,freq)
            pd = prob.BinomialDistribution(Np(1),Np(2));
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [true false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.BinomialDistribution');
            info.name = prob.BinomialDistribution.DistributionName;
            info.code = 'binomial';
            info.prequired = [true false];
            info.support = [0 Inf];
            info.closedbound = [true false];
            info.iscontinuous = false;
            info.logci = [false false];
            info.plim = [1 0;Inf 1];
        end
    end
end % classdef

function checkargs(N,p,x)
if ~internal.stats.isScalarInt(N,0)
    error(message('stats:probdists:NonnegativeIntegerParameter','N'))
elseif ~(isscalar(p) && isnumeric(p) && isreal(p) && p>=0 && p<=1)
    error(message('stats:probdists:Parameter0To1','P'))
elseif any(x<0)
    error(message('stats:binofit:InvalidX'))
elseif any(x>N)
    error(message('stats:binofit:InvalidN'))
end

end
