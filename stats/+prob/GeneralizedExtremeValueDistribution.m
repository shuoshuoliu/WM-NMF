classdef GeneralizedExtremeValueDistribution < prob.ToolboxFittableParametricDistribution 
%GeneralizedExtremeValueDistribution Generalized extreme value probability distribution.
%    An object of the GeneralizedExtremeValueDistribution class represents
%    generalized extreme value probability distribution with specific values
%    of the K, SIGMA, and MU parameters. This distribution object can be created
%    directly using the MAKEDIST function or fit to data using the FITDIST function.
%
%    GeneralizedExtremeValueDistribution methods:
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
%    GeneralizedExtremeValueDistribution properties:    
%       DistributionName      - Name of the distribution
%       k                     - Value of the k parameter (shape)
%       sigma                 - Value of the sigma parameter (scale)
%       mu                    - Value of the mu parameter (location)
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

%    Copyright 2012-2014 The MathWorks, Inc.

    properties(Dependent=true)
%K Shape parameter
%    The K property represents the shape parameter of the generalized
%    extreme value distribution.
%
%    See also SIGMA, MU.
        k

%SIGMA Scale parameter
%    The SCALE property represents the scale parameter of the generalized
%    extreme value distribution.
%
%    See also K, MU.
        sigma

%MU Location parameter
%    The MU property represents the location parameter of the generalized
%    extreme value distribution.
%
%    See also K, SIGMA.
        mu
    
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameGeneralizedExtremeValue'));

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
        ParameterNames = {'k' 'sigma' 'mu'};

%ParameterDescription Parameter description.
%    ParameterNames is a cell array of strings containing short
%    descriptions of the parameters of the probability distribution.
%
%    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionShape')) ...
                                getString(message('stats:probdists:ParameterDescriptionScale')) ...
                                getString(message('stats:probdists:ParameterDescriptionLocation'))};
    end
    properties(GetAccess='public',SetAccess='protected')
%ParameterValues Parameter values.
%    ParameterVales is a vector containing the values of the parameters of
%    the probability distribution.
%
%    See also K, SIGMA, MU.
        ParameterValues
    end
    methods(Hidden)
        function pd = GeneralizedExtremeValueDistribution(k,sigma,mu)
            if nargin==0
                k = 0;
                sigma = 1;
                mu = 0;
            end
            checkargs(k,sigma,mu)

            pd.ParameterValues = [k sigma mu];
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
            m = gevstat(this.k, this.sigma, this.mu);
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            [~,v] = gevstat(this.k, this.sigma, this.mu);
        end
        function this = set.k(this,k)
            checkargs(k,this.sigma,this.mu);
            this.ParameterValues(1) = k;
            this = invalidateFit(this);
        end
        function this = set.sigma(this,sigma)
            checkargs(this.k,sigma,this.mu);
            this.ParameterValues(2) = sigma;
            this = invalidateFit(this);
        end
        function this = set.mu(this,mu)
            checkargs(this.k,this.sigma,mu);
            this.ParameterValues(3) = mu;
            this = invalidateFit(this);
        end
        function mu = get.k(this)
            mu = this.ParameterValues(1);
        end
        function sigma = get.sigma(this)
            sigma = this.ParameterValues(2);
        end
        function sigma = get.mu(this)
            sigma = this.ParameterValues(3);
        end

    end
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the generalized extreme value
%    distribution to data. You should call the FITDIST function
%    instead of calling this method directly.
%
%    See also FITDIST.
            [xOriginal,cens,freq,opt] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            x = prob.ToolboxFittableParametricDistribution.removeCensoring(xOriginal,cens,freq,'generalized extreme value');
            p = gevfit(x,0.05,opt);
            [nll,cov] = gevlike(p,x);
            pd = prob.GeneralizedExtremeValueDistribution.makeFitted(p,nll,cov,xOriginal,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = gevlike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = gevcdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = gevpdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = gevinv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = gevrnd(varargin{:});
        end
        function pd = makeFitted(p,nll,cov,x,cens,freq)
            pd = prob.GeneralizedExtremeValueDistribution(p(1),p(2),p(3));
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.GeneralizedExtremeValueDistribution');
            info.name = getString(message('stats:dfittool:NameGeneralizedExtremeValue'));
            info.code = 'generalized extreme value';
            info.logci = [false true false];
        end
    end
end % classdef

function checkargs(k,sigma,mu)
if ~(isscalar(k) && isnumeric(k) && isreal(k) && isfinite(k))
    error(message('stats:probdists:ScalarParameter','K'))
end
if ~(isscalar(sigma) && isnumeric(sigma) && isreal(sigma) && sigma>=0 && isfinite(sigma))
    error(message('stats:probdists:NonnegativeParameter','SIGMA'))
end
if ~(isscalar(mu) && isnumeric(mu) && isreal(mu) && isfinite(mu))
    error(message('stats:probdists:ScalarParameter','MU'))
end
end

% function [range,closed] = localgevsupport(params)
% k = params(1);
% sigma = params(2);
% mu = params(3);
% 
% if k==0
%     range = [-Inf Inf];
%     closed = [false false];
% elseif k>0
%     range = [mu-sigma/k, Inf];
%     closed = [true false];
% else
%     range = [-Inf, mu-sigma/k];
%     closed = [false true];
% end
% end
