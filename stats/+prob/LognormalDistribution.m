classdef LognormalDistribution < prob.ToolboxFittableParametricDistribution
%LognormalDistribution Lognormal probability distribution.
%    An object of the LognormalDistribution class represents a lognormal
%    probability distribution with specific parameters MU and SIGMA.
%    This distribution object can be created directly using the MAKEDIST
%    function or fit to data using the FITDIST function.
%
%    LognormalDistribution methods:
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
%    LognormalDistribution properties:    
%       DistributionName      - Name of the distribution
%       mu                    - Value of the mu parameter (log mean)
%       sigma                 - Value of the sigma parameter (log standard deviation)
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

    properties(GetAccess='public',Constant=true)
    end
    properties(Dependent=true)
%MU Value of MU parameter
%    The MU property represents the parameter that is the mean of the
%    log of X, when X has a lognormal distribution.
%
%    See also SIGMA.
        mu
        
%SIGMA Value of SIGMA parameter
%    The SIGMA property represents the parameter that is the standard
%    deviation of the log of X, when X has a lognormal distribution.
%
%    See also MU.
        sigma
    end
    properties(GetAccess='public',Constant=true)
%NumParameter Number of parameters.
%    NumParameters is the number of parameters in the distribution.
%
%    See also ParameterValues.
        NumParameters = 2;

%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameLognormal'));

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
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionLogLocation')) ...
                                getString(message('stats:probdists:ParameterDescriptionLogScale'))};
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
        function pd = LognormalDistribution(mu,sigma)
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
            m = lognstat(this.mu, this.sigma);
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            [~,v]= lognstat(this.mu, this.sigma);
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
            [alpha,pnums,citype] = internal.stats.parseArgs(okargs,defaults,varargin{:});
            if ~(isscalar(alpha) && isnumeric(alpha) && alpha>0 && alpha<1)
                error(message('stats:ecdf:BadAlpha'))
            end
            if isequal(citype,'exact')  && ~all(this.ParameterIsFixed)
                pnums = pNamesToNums(this,pnums);
                ci = dfswitchyard('statnormci',[this.mu this.sigma],this.ParameterCovariance,alpha,...
                            this.InputData.data,this.InputData.cens,this.InputData.freq);
                ci = ci(:,pnums);
            else
                ci = paramci@prob.ToolboxFittableParametricDistribution(this,varargin{:});
            end
        end
    end
    methods
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
%    FIT is a static method that fits the lognormal distribution to data.
%    Fitting requires the Statistics and Machine Learning Toolbox. You should call the FITDIST
%    function instead of calling this method directly.
%
%    See also FITDIST.
            
            [x,cens,freq,opt] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            p = lognfit(x,0.05,cens,freq,opt);
            [nll,cov] = lognlike(p,x,cens,freq);
            pd = prob.LognormalDistribution.makeFitted(p,nll,cov,x,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = lognlike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = logncdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = lognpdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = logninv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = lognrnd(varargin{:});
        end
        function pd = makeFitted(p,nll,cov,x,cens,freq)
            pd = prob.LognormalDistribution(p(1),p(2));
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.LognormalDistribution');
            info.name = prob.LognormalDistribution.DistributionName;
            info.code = 'lognormal';
            info.hasconfbounds = true;
            info.censoring = true;
            info.islocscale = true;
            info.uselogpp = true;
            info.logci = [false true];
        end
        function name = matlabCodegenRedirect(~) % redirect to the codegen class
            name = 'prob.coder.LognormalDistribution';
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
