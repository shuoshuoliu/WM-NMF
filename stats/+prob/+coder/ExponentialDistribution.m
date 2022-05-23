classdef ExponentialDistribution < prob.coder.FittableParametricDistribution %#codegen %#internal
    %
    
    %   Copyright 2019 The MathWorks, Inc.
    
    % Code generation version of the ExponentialDistribution class.
    %    An object of the ExponentialDistribution class represents an exponential
    %    probability distribution with specific values of the parameter MU.
    %    This distribution object can be fit to data using the FITDIST function.
    
    %
    %    ExponentialDistribution methods:
    %       cdf                   - Cumulative distribution function
    %       icdf                  - Inverse cumulative distribution function
    %       iqr                   - Interquartile range
    %       mean                  - Mean
    %       median                - Median
    %       pdf                   - Probability density function
    %       std                   - Standard deviation
    %       truncate              - Truncation distribution to an interval
    %       var                   - Variance
    %
    %    ExponentialDistribution properties:
    %       DistributionName      - Name of the distribution
    %       mu                    - Value of the mu parameter (mean)
    %       NumParameters         - Number of parameters
    %       ParameterValues       - Vector of values of parameters
    %       Truncation            - Two-element vector indicating truncation limits
    %       IsTruncated           - Boolean flag indicating if distribution is truncated
    %       InputData             - Structure containing data used to fit the distribution
    %
    %    See also fitdist.
    
    properties (Dependent=true)
        %MU Defining parameter of exponential distribution
        %    The MU property represents the defining parameter, and the mean,
        %    of the exponential distribution
        %
        %    See also ParameterValues.
        mu
    end
    properties (Access='public', Constant = true)
        NumParameters = 1;
    end
    properties(GetAccess='public',SetAccess='protected')
        %ParameterValues Parameter values.
        %    ParameterVales is a vector containing the values of the parameters of
        %    the probability distribution.
        ParameterValues
    end
    
    methods (Hidden)
        function pd = ExponentialDistribution(mu)
            if nargin==0 % for debugging purpose, only
                mu = 1;
            end
            pd.ParameterValues = mu;
            pd.Truncation = zeros(1,2,'like',mu); % We are initializing the Truncation property to maintain a type and size throughout.
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
            m = this.mu;
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            v = this.mu^2;
        end
        
        function mu = get.mu(this)
            mu = this.ParameterValues(1);
        end
    end
    
    methods (Static, Hidden)
        
        function pd = fit(x,cens,freq,~)
            %FIT Fit distribution to data.
            %    FIT is a static method that fits the exponential distribution to data. You
            %    should call the FITDIST function instead of calling this method
            %    directly.
            %
            %    See also FITDIST.
            mu = expfit(x,0.05,cens,freq);
            pd = prob.coder.ExponentialDistribution(mu);
            if coder.target('MEX') % for MEX targets, this data will be needed in reconstruction
                pd.InputData = struct('Data',x,'Censoring',cens,'Frequency',freq);
            end
        end
        
        function coderobj = matlabCodegenToRedirected(mlobj) % runs in MATLAB
            if isempty(mlobj.InputData) % if the object was not returned from fitdist
                error(message('stats:probdists:CodegenArgsNotFromFitdist'));
            end
            coderobj = prob.coder.ExponentialDistribution(mlobj.mu) ;
            
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
            if ~coderobj.IsTruncated  % an untruncated distribution means it is the output of fitdist, not truncate()
                x = coderobj.InputData.Data;
                cens = coderobj.InputData.Censoring;
                freq = coderobj.InputData.Frequency;
                [nll, cov] = explike(coderobj.ParameterValues,x,cens,freq);
                mlobj = prob.ExponentialDistribution.makeFitted(coderobj.ParameterValues ,nll,cov,...
                    x,cens,freq);
            else
                mlobj = prob.ExponentialDistribution(coderobj.ParameterValues(1));
            end
            if coderobj.IsTruncated
                mlobj = truncate(mlobj,coderobj.Truncation(1),coderobj.Truncation(2));
            end
        end
        
        function y = invfunc(varargin)
            y = expinv(varargin{:});
            %             [varargout{1:nargout}] = expinv(varargin{:});
        end
        function y = cdffunc(varargin)
            y = expcdf(varargin{:});
            %             [varargout{1:nargout}] = expcdf(varargin{:});
        end
        function y = pdffunc(varargin)
            y = exppdf(varargin{:});
            %             [varargout{1:nargout}] = exppdf(varargin{:});
        end
        
    end % static, hidden methods
end  % classdef
