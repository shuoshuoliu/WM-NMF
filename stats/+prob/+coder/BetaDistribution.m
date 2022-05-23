classdef BetaDistribution < prob.coder.FittableParametricDistribution %#codegen %#internal
    %
    
    %   Copyright 2019 The MathWorks, Inc.
    
    % Code generation version of the BetaDistribution class.
    %    An object of the BetaDistribution class represents a beta
    %    probability distribution with specific values of the A and P
    %    parameters.
    %    This distribution object can be fit to data using the FITDIST function.
    
    %
    %    BetaDistribution methods:
    %       cdf                   - Cumulative distribution function
    %       icdf                  - Inverse cumulative distribution function
    %       iqr                   - Interquartile range
    %       mean                  - Mean
    %       median                - Median
    %       pdf                   - Probability density function
    %       random                - Random number generation
    %       std                   - Standard deviation
    %       truncate              - Truncation distribution to an interval
    %       var                   - Variance
    %
    %    BetaDistribution properties:
    %       a                     - Value of the a parameter
    %       b                     - Value of the b parameter
    %       NumParameters         - Number of parameters
    %       ParameterValues       - Vector of values of parameters
    %       Truncation            - Two-element vector indicating truncation limits
    %       IsTruncated           - Boolean flag indicating if distribution is truncated
    %       InputData             - Structure containing data used to fit the distribution
    %
    %    See also fitdist.
    
    
    properties(Dependent=true)
        %A First parameter of beta distribution
        %    The A property represents the first of the two shape parameters of the
        %    beta distribution.
        %
        %    See also B.
        a
        
        %B Second parameter of beta distribution
        %    The B property represents the second of the two shape parameters of
        %    the beta distribution.
        %
        %    See also A.
        b
    end
    properties(GetAccess='public',Constant=true)
        %DistributionName Distribution name.
        %    The DistributionName property indicates the name of the probability
        %    distribution.
        %
        %    See also ParameterValues.
        
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
        %    See also A, B.
        ParameterValues
    end
    
    methods(Hidden)
        function pd = BetaDistribution(a,b)
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
            A = this.a;
            B = this.b;
            m = A / (A+B);
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            A = this.a;
            B = this.b;
            v = A*B / ((A+B)^2 * (A+B+1));
        end
        function a = get.a(this)
            a = this.ParameterValues(1);
        end
        function b = get.b(this)
            b = this.ParameterValues(2);
        end
    end
    methods(Static,Hidden)
        function pd = fit(x,cens,freq,~)
            %FIT Fit distribution to data.
            %    FIT is a static method that fits the beta distribution to data. You
            %    should call the FITDIST function instead of calling this method
            %    directly.
            %
            %    See also FITDIST.
            coder.internal.errorIf(any(cens), 'stats:ProbDistUnivParam:fit:CensoringNotAllowed', 'beta');
            if ~isempty(freq)
                coder.internal.warning('stats:probdists:FrequencyInputIgnored','beta');
            end
            
            p = betafit(x,0.05);
            pd = prob.coder.BetaDistribution(p(1),p(2));
            if coder.target('MEX') % for MEX targets, this data will be needed in reconstruction
                pd.InputData = struct('Data',x,'Censoring',cens,'Frequency',freq);
            end
        end
        
        function coderobj = matlabCodegenToRedirected(mlobj) % runs in MATLAB
            if isempty(mlobj.InputData) % if the object was not returned from fitdist
                error(message('stats:probdists:CodegenArgsNotFromFitdist'));
            end
            coderobj = prob.coder.BetaDistribution(mlobj.a, mlobj.b) ;
            
            % There are two scenarios where reconstruction of the object is needed
            % for MEX targets
            % (1) output from fitdist : this returns a distribution object in
            % MATLAB which contains a nonempty InputData
            % (2) output from truncate : this returns a truncated distribution
            % object in MATLAB, which contains empty InputData
            % We are handling both the cases in the
            % matlabCodegenFromRedirected, the InputData in the coderobj has
            % nonempty fields in case (1) and empty fields in case (2).
        end
        
        function mlobj = matlabCodegenFromRedirected(coderobj) % runs in MATLAB
            if ~coderobj.IsTruncated % an untruncated distribution means it is the output of fitdist, not truncate()
                x = coderobj.InputData.Data;
                cens = coderobj.InputData.Censoring;
                freq = coderobj.InputData.Frequency;
                
                [nll, cov] = betalike(coderobj.ParameterValues,x);
                mlobj = prob.BetaDistribution.makeFitted(coderobj.ParameterValues ,nll,cov,...
                    x,cens,freq);
            else
                mlobj = prob.BetaDistribution(coderobj.ParameterValues(1), ...
                    coderobj.ParameterValues(2));
            end
            
            if coderobj.IsTruncated
                mlobj = truncate(mlobj,coderobj.Truncation(1),coderobj.Truncation(2));
            end
        end
        
        function y = cdffunc(varargin)
            y = betacdf(varargin{:});
            %                         [varargout{1:nargout}] = betacdf(varargin{:});
            
        end
        function y = pdffunc(varargin)
            y = betapdf(varargin{:});
            %             [varargout{1:nargout}] = betapdf(varargin{:});
        end
        function y = invfunc(varargin)
            y = betainv(varargin{:});
            %             [varargout{1:nargout}] = betainv(varargin{:});
        end
    end % static, hidden methods
end % classdef
