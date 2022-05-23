classdef WeibullDistribution < prob.ToolboxFittableParametricDistribution
%WeibullDistribution Weibull probability distribution.
%    An object of the WeibullDistribution class represents a Weibull
%    probability distribution with specific scale parameter A and shape
%    parameter B.  This distribution object can be created directly using 
%    the MAKEDIST function or fit to data using the FITDIST function.
%
%    WeibullDistribution methods:
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
%    WeibullDistribution properties:    
%       DistributionName      - Name of the distribution
%       A                     - Value of the A parameter (scale)
%       B                     - Value of the B parameter (shape)
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
%A Value of A parameter
%    The A property represents the scale parameter of the Weibull distribution.
%
        A 
%B Value of B parameter
%    The B property represents the shape parameter of the Weibull distribution.
%
        B 
    end
    properties(Hidden,Dependent=true)
        a
        b
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameWeibull'));

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
        ParameterNames = {'A' 'B'};

%ParameterDescription Parameter description.
%    ParameterNames is a cell array of strings containing short
%    descriptions of the parameters of the probability distribution.
%
%    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionScale')) ...
                                getString(message('stats:probdists:ParameterDescriptionShape'))};;
    end
    properties(GetAccess='public',SetAccess='protected')
%ParameterValues Parameter values.
%    ParameterVales is a vector containing the values of the parameters of
%    the probability distribution.
%
%    See also A,B.
        ParameterValues
    end
    methods(Hidden)
        function pd = WeibullDistribution(a,b)
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
            m = wblstat(this.A, this.B);
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            [~,v] = wblstat(this.A, this.B);
        end
        function ci = paramci(this,varargin)
            [varargin{:}] = convertStringsToChars(varargin{:});
            
            requireScalar(this)
            if isscalar(varargin) && isnumeric(varargin{1})
                % Support syntax of older ProbDistUnivParam/paramci method
                varargin = {'alpha' varargin{1}};
            end
            logci = this.getInfo.logci;
            okargs =   {'alpha' 'parameter'          'type' 'logflag'};
            defaults = {0.05,   1:this.NumParameters 'wald' logci};
            [alpha,pnums,citype] = internal.stats.parseArgs(okargs,defaults,varargin{:});
            
            if all(this.ParameterIsFixed)
                pnums = pNamesToNums(this,pnums);
                ci = repmat(this.ParameterValues(pnums),2,1);
            elseif isequal(citype,'wald') && all(logci==this.getInfo.logci)
                % Convert parameter values to the ev values used in their
                % computation
                %
                parmhat = [log(this.ParameterValues(1)), 1/(this.ParameterValues(2))];
                
                % Calculate ci values on the ev scale
                %
                probs = [alpha/2; 1-alpha/2];
                [~, acov] = evlike(parmhat, log(this.InputData.data), ...
                    this.InputData.cens, this.InputData.freq);
                % The first parameter is already on a log scale.  Also
                % work on a log scale for the second parameter.
                transfhat = [parmhat(1) log(parmhat(2))];
                se = sqrt(diag(acov))';
                se(2) = se(2) ./ parmhat(2); % se(log(sigmahat)) = se(sigmahat) / sigmahat
                ciEV = norminv([probs, probs], [transfhat; transfhat], [se; se]);
                % Reverse the log transform the second parameter.
                ciEV(:,2) = exp(ciEV(:,2));
                
                % Convert back to weibull ci values
                %
                ci = [exp(ciEV(:,1)), 1./ciEV([2 1],2)];

                % Return confidence intervals for the parameters for which they were
                % requested
                pnums = pNamesToNums(this,pnums);
                ci = ci(:,pnums);
            else
                ci = paramci@prob.ToolboxFittableParametricDistribution(this,varargin{:});
            end
        end
    end
    methods
        function this = set.A(this,a)
            checkargs(a,this.B);
            this.ParameterValues(1) = a;
            this = invalidateFit(this);
        end
        function a = get.A(this)
            a = this.ParameterValues(1);
        end
        function this = set.B(this,b)
            checkargs(this.A,b);
            this.ParameterValues(2) = b;
            this = invalidateFit(this);
        end
        function b = get.B(this)
            b = this.ParameterValues(2);
        end
        function this = set.a(this,a)
            checkargs(a,this.B);
            this.ParameterValues(1) = a;
            this = invalidateFit(this);
        end
        function a = get.a(this)
            a = this.ParameterValues(1);
        end
        function this = set.b(this,b)
            checkargs(this.A,b);
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
%    FIT is a static method that fits the Weibull distribution to data.
%    Fitting requires the Statistics and Machine Learning Toolbox. You should call the FITDIST
%    function instead of calling this method directly.
%
%    See also FITDIST.

            [xOriginal,cens,freq,opts] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            params = wblfit(xOriginal,0.05,cens,freq,opts);
            [nll,cov] = prob.WeibullDistribution.likefunc(params,xOriginal,cens,freq);
            pd = prob.WeibullDistribution.makeFitted(params,nll,cov,xOriginal,cens,freq);
        end
        function varargout = likefunc(varargin)
           [varargout{1:nargout}] = wbllike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = wblcdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = wblpdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = wblinv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = wblrnd(varargin{:});
        end
        function pd = makeFitted(params,nll,cov,x,cens,freq)
            pd = prob.WeibullDistribution(params(1),params(2));
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.WeibullDistribution');
            info.name = getString(message('stats:dfittool:NameWeibull'));
            info.code = 'weibull';
            info.hasconfbounds = true;
            info.censoring = true;
            info.support = [0 Inf];
            info.islocscale = true;
            info.uselogpp = true;
            info.optimopts = true;
            info.logci = [true true];
        end
        function name = matlabCodegenRedirect(~) % redirect to the codegen class
                         name = 'prob.coder.WeibullDistribution';
        end
    end
end % classdef

function checkargs(a,b)
if ~(isscalar(a) && isnumeric(a) && isreal(a) && a>0 && isfinite(a))
    error(message('stats:probdists:PositiveParameter','A'))
end
if ~(isscalar(b) && isnumeric(b) && isreal(b) && b>0 && ~isnan(b))
    error(message('stats:probdists:PositiveParameter','B'))
end
end

