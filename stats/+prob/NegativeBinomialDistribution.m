classdef NegativeBinomialDistribution < prob.ToolboxFittableParametricDistribution
%NegativeBinomialDistribution Negative binomial probability distribution.
%    An object of the NegativeBinomialDistribution class represents a Negative
%    Binomial probability distribution with specific values of the parameters
%    R and P. This distribution object can be created directly using the 
%    MAKEDIST function or fit to data using the FITDIST function.
%
%    NegativeBinomialDistribution methods:
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
%    NegativeBinomialDistribution properties:    
%       DistributionName      - Name of the distribution
%       R                     - Value of the R parameter (#successes)
%       P                     - Value of the P parameter (Prob(success))
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
%R Number of Successes parameter of Negative Binomial distribution
%    The R property represents the parameter of the Negative Binomial distribution
%    that corresponds to the number of successes.
%
%    See also ParameterNames.
        R;
%P Probability of Success parameter of Negative Binomial distribution
%    The P property represents the parameter of the Negative Binomial distribution
%    that corresponds to the probability of success of an individual trial.
%
%    See also ParameterNames.
        P;
    end
    properties(Hidden,Dependent=true)
        r
        p
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameNegativeBinomial'));

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
        ParameterNames = {'R', 'P'};

%ParameterDescription Parameter description.
%    ParameterNames is a cell array of strings containing short
%    descriptions of the parameters of the probability distribution.
%
%    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionSuccesses')) ...
                                getString(message('stats:probdists:ParameterDescriptionProbability'))};
    end
    
    properties(GetAccess='public',SetAccess='protected')
%ParameterValues Parameter values.
%    ParameterVales is a vector containing the values of the parameters of
%    the probability distribution.
%
%    See also R, P.
        ParameterValues
    end
    
    methods(Hidden)
        function pd = NegativeBinomialDistribution(r,p)
            if nargin==0
                r = 1;
                p = 0.5;
            end
            checkargs(r,p)
            
            pd.ParameterValues = [r, p];
            pd.ParameterIsFixed = [true, true];
            pd.ParameterCovariance = zeros(pd.NumParameters);
        end
    end
    methods
        function m = mean(this)
            requireScalar(this)
            m = nbinstat(this.R,this.P);
            mu = m;
            if this.IsTruncated
                % Determine the range of integers that the truncation
                % points contain.
                [plower, pupper, lower, upper] = tailprobs(this);
                lower = ceil(lower);
                upper = floor(upper);
                if isequal(upper,Inf)
                    % Can't sum infinite series.  Instead, subtract off
                    % the left tail from the adjusted untruncated mean.
                    w = 1 / (pupper - plower);
                    x = 1:(lower-1);
                    m = w * (mu - sum(x .* pdffun(this,x)));
                else
                    x = lower:upper;
                    p = pdf(this,x);
                    m = sum(p .* x);
                end
                return
            end
        end
        function v = var(this)
            requireScalar(this)
            [mu,v] = nbinstat(this.R,this.P);
            if this.IsTruncated
                m = mean(this);
                % Determine the range of integers that the truncation
                % points contain.
                [plower, pupper, lower, upper] = tailprobs(this);
                lower = ceil(lower);
                upper = floor(upper);
                if isequal(upper,Inf)
                    % Can't sum infinite series.  Instead, subtract off
                    % the left tail from the adjusted untruncated variance.
                    w = 1 / (pupper - plower);
                    x = 0:(lower-1);
                    v = w * (v + (mu-m)^2 - sum( (x-m).^2 .* pdffun(this,x) ));
                else
                    x = lower:upper;
                    p = pdf(this,x);
                    v = sum(p .* (x-m).^2);
                end
                return
            end
        end
    end
    methods
        function this = set.R(this,R)
            checkargs(R,this.P)
            this.ParameterValues(1) = R;
            this = invalidateFit(this);
        end
        function R = get.R(this)
            R = this.ParameterValues(1);
        end
        function this = set.P(this,P)
            checkargs(this.R,P)
            this.ParameterValues(2) = P;
            this = invalidateFit(this);
        end
        function P = get.P(this)
            P = this.ParameterValues(2);
        end
        function this = set.r(this,R)
            checkargs(R,this.P)
            this.ParameterValues(1) = R;
            this = invalidateFit(this);
        end
        function R = get.r(this)
            R = this.ParameterValues(1);
        end
        function this = set.p(this,P)
            checkargs(this.R,P)
            this.ParameterValues(2) = P;
            this = invalidateFit(this);
        end
        function P = get.p(this)
            P = this.ParameterValues(2);
        end
    end
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the Negative Binomial distribution to data. You
%    should call the FITDIST function instead of calling this method directly.
%
%    See also FITDIST.
            [xOriginal,cens,freq,opt] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            x = prob.ToolboxFittableParametricDistribution.removeCensoring(xOriginal,cens,freq,'nbinon');
            params = nbinfit(x,0.05,opt);
            [nll,cov] = nbinlike(params,x);
            pd = prob.NegativeBinomialDistribution.makeFitted(params,nll,cov,xOriginal,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = nbinlike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = nbincdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = nbinpdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = nbininv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = nbinrnd(varargin{:});
        end
        function pd = makeFitted(params,nll,cov,x,cens,freq)
            c = num2cell(params);
            pd = prob.NegativeBinomialDistribution(c{:});
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.NegativeBinomialDistribution');
            info.name = getString(message('stats:dfittool:NameNegativeBinomial'));
            info.code = 'negative binomial';
            info.support = [0 Inf];
            info.closedbound = [true false];
            info.iscontinuous = false;
            info.optimopts = true;
            info.plim = [0 0;Inf 1];
        end
    end
end % classdef

function checkargs(r,p)
if ~(isscalar(r) && isnumeric(r) && isreal(r) && r>0 && isfinite(r))
    error(message('stats:probdists:PositiveParameter','R'))
end
if ~(isscalar(p) && isnumeric(p) && isreal(p) && p>0 && p<=1)
    error(message('stats:probdists:PositiveParameterLE1','P','P'))
end
end
