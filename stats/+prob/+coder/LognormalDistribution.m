classdef LognormalDistribution <  prob.coder.FittableParametricDistribution %#codegen
    %
    
    %   Copyright 2019 The MathWorks, Inc.
    
    % Code generation version of the LognormalDistribution class.
    %    An object of the LognormalDistribution class represents a lognormal
    %    probability distribution with specific parameters MU and SIGMA.
    %    This distribution object can be fit to data using the FITDIST function.
    %
    %    LognormalDistribution methods:
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
    %    LognormalDistribution properties:
    %       mu                    - Value of the mu parameter (log mean)
    %       sigma                 - Value of the sigma parameter (log standard deviation)
    %       NumParameters         - Number of parameters
    %       ParameterValues       - Vector of values of parameters
    %       Truncation            - Two-element vector indicating truncation limits
    %       IsTruncated           - Boolean flag indicating if distribution is truncated
    %       InputData             - Structure containing data used to fit the distribution
    %
    %    See also fitdist.
    
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
            if nargin==0 % for debugging purpose only
                mu = 0;
                sigma = 1;
            end
            pd.ParameterValues = [mu sigma];
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
    end
    methods
        
        function mu = get.mu(this)
            mu = this.ParameterValues(1);
        end
        function sigma = get.sigma(this)
            sigma = this.ParameterValues(2);
        end
        
    end
    methods(Static,Hidden)
        function pd = fit(x,cens,freq,varargin)
            %FIT Fit distribution to data.
            %    FIT is a static method that fits the lognormal distribution to data.
            %    Fitting requires the Statistics and Machine Learning Toolbox. You should call the FITDIST
            %    function instead of calling this method directly.
            %
            %    See also FITDIST.
            if ~isempty(varargin)
                opt = varargin{1};
            else
                opt = statset('lognfit');
            end
            p = lognfit(x,0.05,cens,freq,opt);
            pd = prob.coder.LognormalDistribution(p(1),p(2));
            
            if coder.target('MEX') %  for MEX targets, this data will be needed in reconstruction
                pd.InputData = struct('Data',x,'Censoring',cens,'Frequency',freq);
            end
        end
        function coderobj = matlabCodegenToRedirected(mlobj)
            if isempty(mlobj.InputData) % if the object was not returned from fitdist
                error(message('stats:probdists:CodegenArgsNotFromFitdist'));
            end
            coderobj = prob.coder.LognormalDistribution(mlobj.mu, mlobj.sigma) ;
            
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
        
        function mlobj = matlabCodegenFromRedirected(coderobj)
            if ~coderobj.IsTruncated % an untruncated distribution means it is the output of fittin
                x = coderobj.InputData.Data;
                cens = coderobj.InputData.Censoring;
                freq = coderobj.InputData.Frequency;
                
                [nll, cov] = lognlike(coderobj.ParameterValues,x,cens,freq);
                mlobj = prob.LognormalDistribution.makeFitted(coderobj.ParameterValues ,nll,cov,...
                    x,cens,freq);
            else
                mlobj = prob.LognormalDistribution(coderobj.ParameterValues(1),...
                    coderobj.ParameterValues(2));
            end
            if coderobj.IsTruncated
                mlobj = truncate(mlobj,coderobj.Truncation(1),coderobj.Truncation(2));
                
            end
        end
        
        function y = cdffunc(varargin)
            y = logncdf(varargin{:});
            %             [varargout{1:nargout}] = logncdf(varargin{:});
        end
        function y = pdffunc(varargin)
            y = lognpdf(varargin{:});
            %             [varargout{1:nargout}] = lognpdf(varargin{:});
        end
        function y = invfunc(varargin)
            y =  logninv(varargin{:});
            %             [varargout{1:nargout}] = logninv(varargin{:});
        end
        
    end % static, hidden methods
end % classdef

