classdef GeneralizedParetoDistribution < prob.ToolboxFittableParametricDistribution 
%GeneralizedParetoDistribution Generalized pareto probability distribution.
%    An object of the GeneralizedParetoDistribution class represents the
%    generalized extreme value probability distribution with specific values
%    of the K, SIGMA, and THETA parameters. This distribution object can be created
%    directly using the MAKEDIST function or fit to data using the FITDIST function.
%
%    GeneralizedParetoDistribution methods:
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
%    GeneralizedParetoDistribution properties:    
%       DistributionName      - Name of the distribution
%       k                     - Value of the k parameter (shape)
%       sigma                 - Value of the sigma parameter (scale)
%       theta                 - Value of the theta parameter (location)
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

%    Copyright 2012-2016 The MathWorks, Inc.

    properties(Dependent=true)
%K Shape parameter
%    The K property represents the shape parameter of the generalized
%    pareto distribution.
%
%    See also SIGMA, THETA.
        k

%SIGMA Scale parameter
%    The SCALE property represents the scale parameter of the generalized
%    pareto distribution.
%
%    See also K, THETA.
        sigma

%MU Location parameter
%    The THETA property represents the threshold parameter of the generalized
%    pareto distribution.
%
%    See also K, SIGMA.
        theta
    
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameGeneralizedPareto'));

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
        ParameterNames = {'k' 'sigma' 'theta'};

%ParameterDescription Parameter description.
%    ParameterNames is a cell array of strings containing short
%    descriptions of the parameters of the probability distribution.
%
%    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionShape')) ...
                                getString(message('stats:probdists:ParameterDescriptionScale')) ...
                                getString(message('stats:probdists:ParameterDescriptionThreshold'))};
    end
    properties(GetAccess='public',SetAccess='protected')
%ParameterValues Parameter values.
%    ParameterVales is a vector containing the values of the parameters of
%    the probability distribution.
%
%    See also K, SIGMA, THETA.
        ParameterValues
    end
    methods(Hidden)
        function pd = GeneralizedParetoDistribution(k,sigma,theta)
            if nargin==0
                k = 1;
                sigma = 1;
                theta = 1;
            end
            checkargs(k,sigma,theta)

            pd.ParameterValues = [k sigma theta];
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
            m = gpstat(this.k, this.sigma, this.theta);
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            [~,v] = gpstat(this.k, this.sigma, this.theta);
        end
        function this = set.k(this,k)
            checkargs(k,this.sigma,this.theta);
            this.ParameterValues(1) = k;
            this = invalidateFit(this);
        end
        function this = set.sigma(this,sigma)
            checkargs(this.k,sigma,this.theta);
            this.ParameterValues(2) = sigma;
            this = invalidateFit(this);
        end
        function this = set.theta(this,theta)
            checkargs(this.k,this.sigma,theta);
            this.ParameterValues(3) = theta;
            this = invalidateFit(this);
        end
        function theta = get.k(this)
            theta = this.ParameterValues(1);
        end
        function sigma = get.sigma(this)
            sigma = this.ParameterValues(2);
        end
        function sigma = get.theta(this)
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

%    PD = PROB.GENERALIZEDPARETODISTRIBUTION.FIT(X) creates a 
%    GeneralizedParetoDistribution object with THETA=0 and with parameters
%    K and SIGMA fit by maximum likelihood to the data in X.
%
%    PD = PROB.GENERALIZEDPARETODISTRIBUTION.FIT(X,NAME,VALUE) specifies
%    one or more of the following optional parameters:
%
%       'Frequency'    A vector of the same size as X, containing
%                      non-negative integer frequencies for the
%                      corresponding elements in X.  Default is one
%                      observation per element of X.
%       'Theta'        A scalar numeric value giving the threshold
%                      parameter of the generalized pareto distribution.
%                      Theta is not estimated from the data.

            [theta,~,fitargs] = internal.stats.parseArgs({'theta'},{0},varargin{2:end});
            [xOriginal,cens,freq,opt] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{1}, fitargs{:});
            x = prob.ToolboxFittableParametricDistribution.removeCensoring(xOriginal,cens,freq,'generalized pareto');
            params = localgpfit(x,theta,0.05,opt);
            [nll,cov] = localgplike(params,x);
            pd = prob.GeneralizedParetoDistribution.makeFitted(params,nll,cov,xOriginal,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = localgplike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = gpcdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = gppdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = gpinv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = gprnd(varargin{:});
        end
        function pd = makeFitted(params,nll,cov,x,cens,freq)
            c = num2cell(params);
            pd = prob.GeneralizedParetoDistribution(c{:});
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false true];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.GeneralizedParetoDistribution');
            info.name = prob.GeneralizedParetoDistribution.DistributionName;
            info.code = 'generalized pareto';
            info.prequired = [false false true];
            info.optimopts = true;
            info.logci = [false true false];
        end
    end
end % classdef

function checkargs(k,sigma,theta)
if ~(isscalar(k) && isnumeric(k) && isreal(k) && isfinite(k))
    error(message('stats:probdists:ScalarParameter','K'))
end
if ~(isscalar(sigma) && isnumeric(sigma) && isreal(sigma) && sigma>=0 && isfinite(sigma))
    error(message('stats:probdists:NonnegativeParameter','SIGMA'))
end
if ~(isscalar(theta) && isnumeric(theta) && isreal(theta) && isfinite(theta))
    error(message('stats:probdists:ScalarParameter','THETA'))
end
end

% ------------ generalized Pareto functions are a special case
function [phat,pci] = localgpfit(x,theta,alpha,varargin)
%LOCALGPFIT Version of gpfit that handles a fixed threshold param

if any(x<=theta)
    error(message('stats:gpfit:BadDataThreshold'));
end

if nargout < 2
    phat = [gpfit(x-theta,alpha,varargin{:}) theta];
else
    [phat,pci] = gpfit(x-theta,alpha);
    phat = [phat theta];
    pci = [pci [theta; theta]];
end
end

function [nlogL,acov] = localgplike(params,data)
%LOCALGPLIKE Version of gplike that handles a fixed threshold param

theta = params(3);
params = params(1:2);
if nargout < 2
    nlogL = gplike(params,data-theta);
else
    [nlogL,acov] = gplike(params,data-theta);
    acov = [acov [0; 0]; [0 0 0]];
end
end

% function [range,closed] = localgpsupport(params)
% k = params(1);
% theta = params(3);
% if k<0
%     sigma = params(2);
%     range = sort([theta, theta-sigma/k]);
%     closed = [true true];
% else
%     range = [theta Inf];
%     closed = [false false];
% end
% end
