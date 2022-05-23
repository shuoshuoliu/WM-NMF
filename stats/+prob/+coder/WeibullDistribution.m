classdef WeibullDistribution <  prob.coder.FittableParametricDistribution %#codegen
    %WeibullDistribution
    %    An object of the WeibullDistribution class represents a Weibull
    %    probability distribution with specific scale parameter A and shape
    %    parameter B.
    %    This distribution object can be fit to data using the FITDIST function.
    
    %   Copyright 2019 The MathWorks, Inc.
    
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
    
    properties(GetAccess='public',Constant=true)
        
        %NumParameter Number of parameters.
        %    NumParameters is the number of parameters in the distribution.
        %
        %    See also ParameterValues.
        NumParameters = 2;
        
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
            if nargin==0 % for debugging purpose only
                a = 1;
                b = 1;
            end
            pd.ParameterValues = [a b];
            pd.Truncation = zeros(1,2,'like',a); % We are initializing the Truncation property to maintain a type and size throughout.
            % All zeros is not an issue because this is a read only
            % property for the user of probability distribution objects.
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
        function a = get.A(this)
            a = this.ParameterValues(1);
        end
        
        function b = get.B(this)
            b = this.ParameterValues(2);
        end
    end
    
    methods(Static,Hidden)
        function pd = fit(x,cens,freq,varargin)
            %FIT Fit distribution to data.
            %    FIT is a static method that fits the Weibull distribution to data.
            %    Fitting requires the Statistics and Machine Learning Toolbox. You should call the FITDIST
            %    function instead of calling this method directly.
            %
            %    See also FITDIST.
            
            
            if ~isempty(varargin)
                opt = varargin{1}; % add more parameters here in the future
            else
                opt = statset('wblfit');
            end
            
            params = wblfit(x,0.05,cens,freq,opt);
            
            pd = prob.coder.WeibullDistribution(params(1),params(2));
            
            if coder.target('MEX')
                pd.InputData = struct('Data',x,'Censoring',cens,'Frequency',freq);
            end
        end
        
        function coderobj = matlabCodegenToRedirected(mlobj)  % runs in MATLAB
            if isempty(mlobj.InputData) % if the object was not returned from fitdist
                error(message('stats:probdists:CodegenArgsNotFromFitdist'));
            end
            coderobj = prob.coder.WeibullDistribution(mlobj.A, mlobj.B) ;
            
            % There are two scenarios of reconstruction of the object is needed
            % for MEX targets
            % (1) output from fitdist : this returns a distribution object in
            % MATLAB which contains a nonempty InputData
            % (2) output from truncate : this returns a truncated distribution
            % object in MATLAB, which contains empty InputData
            % We are handling both the cases in the
            % matlabCodegenFromRedirected, the InputData in the coderobj has
            % nonempty fields in case (1) and empty fields in case (2).
        end
        
        function mlobj = matlabCodegenFromRedirected(coderobj)  % runs in MATLAB
            % if the target is MATLAB return the entire object
            if ~coderobj.IsTruncated % an untruncated distribution means it is the output of fitdist, not truncate()
                x = coderobj.InputData.Data;
                cens = coderobj.InputData.Censoring;
                freq = coderobj.InputData.Frequency;
                
                [nll, cov] = wbllike(coderobj.ParameterValues,x,cens,freq);
                mlobj = prob.WeibullDistribution.makeFitted(coderobj.ParameterValues ,nll,cov,...
                    x,cens,freq);
            else
                mlobj = prob.WeibullDistribution(coderobj.ParameterValues(1), ...
                    coderobj.ParameterValues(2));
            end
            if coderobj.IsTruncated
                mlobj = truncate(mlobj,coderobj.Truncation(1),coderobj.Truncation(2));
            end
        end
        function y = cdffunc(varargin)
            y = wblcdf(varargin{:});
            %             [varargout{1:nargout}] = wblcdf(varargin{:});
        end
        function y = pdffunc(varargin)
            y =  wblpdf(varargin{:});
            %             [varargout{1:nargout}] = wblpdf(varargin{:});
        end
        function y = invfunc(varargin)
            y = wblinv(varargin{:});
            %             [varargout{1:nargout}] = wblinv(varargin{:});
        end
        
    end % static, hidden methods
end % classdef


