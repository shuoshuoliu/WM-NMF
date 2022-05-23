classdef ExtremeValueDistribution < prob.ToolboxFittableParametricDistribution
%ExtremeValueDistribution Extreme Value probability distribution.
%    An object of the ExtremeValueDistribution class represents an extreme
%    value probability distribution with specific values of the MU and SIGMA
%    parameters. This distribution object can be created directly 
%    using the MAKEDIST function or fit to data using the FITDIST function.
%
%    ExtremeValueDistribution methods:
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
%    ExtremeValueDistribution properties:    
%       DistributionName      - Name of the distribution
%       mu                    - Value of the mu parameter (location)
%       sigma                 - Value of the sigma parameter (scale)
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

%   Copyright 2012-2019 The MathWorks, Inc.

    properties(Dependent=true)
%MU Location parameter of the extreme value distribution
%    The MU property represents the location parameter of the extreme value
%    distribution
%
%    See also ParameterValues.
        mu
%SIGMA Scale parameter of the extreme value distribution
%    The SIGMA property represents the scale parameter of the extreme value
%    distribution
%
%    See also ParameterValues.
        sigma
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameExtremeValue'));

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
        function pd = ExtremeValueDistribution(mu,sigma)
            if nargin==0
                mu = 0;
                sigma = 1;
            end
            checkargs(mu, sigma)
          
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
            m = evstat(this.mu,this.sigma);
        end

        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            [~,v] = evstat(this.mu,this.sigma);
        end
        function this = set.mu(this,mu)
            checkargs(mu,this.sigma)
            this.ParameterValues(1) = mu;
            this = invalidateFit(this);
        end
        function mu = get.mu(this)
            mu = this.ParameterValues(1);
        end
        function this = set.sigma(this,sigma)
            checkargs(this.mu,sigma)
            this.ParameterValues(2) = sigma;
            this = invalidateFit(this);
        end
        function mu = get.sigma(this)
            mu = this.ParameterValues(2);
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
            [x,cens,freq,options] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            p = evfit(x,0.05,cens,freq,options);
            [nll,cov] = evlike(p,x,cens,freq);
            pd = prob.ExtremeValueDistribution.makeFitted(p,nll,cov,x,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = evlike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = evcdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = evpdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = evinv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = evrnd(varargin{:});
        end
        function pd = makeFitted(p,nll,cov,x,cens,freq)
            pd = prob.ExtremeValueDistribution(p(1),p(2));
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.ExtremeValueDistribution');
            info.name = prob.ExtremeValueDistribution.DistributionName;
            info.code = 'extreme value';
            info.hasconfbounds = true;
            info.censoring = true;
            info.islocscale = true;
            info.optimopts = true;
            info.logci = [false true];
        end
        function name = matlabCodegenRedirect(~) % redirect to the codegen class
            name = 'prob.coder.ExtremeValueDistribution';
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
