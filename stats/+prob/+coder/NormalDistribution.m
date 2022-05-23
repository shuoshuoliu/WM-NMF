classdef NormalDistribution < prob.coder.FittableParametricDistribution %#codegen %#internal
%

%   Copyright 2019 The MathWorks, Inc.
    
    %  % Code generation version of the NormalDistribution object
    %    An object of the NormalDistribution class represents a normal
    %    probability distribution with a specific mean MU and standard
    %    deviation SIGMA.
    %    This distribution object can be fit to data using the FITDIST function.
    
    %
    %    NormalDistribution methods:
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
    %    NormalDistribution properties:
    %       mu                    - Value of the mu parameter (mean)
    %       sigma                 - Value of the sigma parameter (standard deviation)
    %       NumParameters         - Number of parameters
    %       ParameterValues       - Vector of values of parameters
    %       Truncation            - Two-element vector indicating truncation limits
    %       IsTruncated           - Boolean flag indicating if distribution is truncated
    %       InputData             - Structure containing data used to fit the distribution
    %
    %    See also fitdist.
        
    properties(Dependent=true)
        %MU Mean of normal distribution
        %    The MU property represents the parameter that is the mean of the
        %    normal distribution.
        %
        %    See also SIGMA.
        mu
        
        %SIGMA Standard deviation of normal distribution
        %    The SIGMA property represents the parameter that is the standard
        %    deviation of the normal distribution.
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
        function pd = NormalDistribution(mu,sigma)
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
            m = this.mu;
        end

        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            v = this.sigma^2;
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
            %    FIT is a static method that fits the normal distribution to data.
            %    You should call the FITDIST function instead of calling this method
            %    directly.
            %
            %    See also FITDIST.
            if ~isempty(varargin)
                opt = varargin{1};
            else
                opt = statset('normfit');
            end
            [m,s] = normfit(x,0.05,cens,freq,opt);
            p = [m s];
            pd = prob.coder.NormalDistribution(p(1),p(2));
            if coder.target('MEX')
                pd.InputData = struct('Data',x,'Censoring',cens,...
                    'Frequency',freq);
            end
        end
        function coderobj = matlabCodegenToRedirected(mlobj)  % runs in MATLAB
            if isempty(mlobj.InputData) % if the object was not returned from fitdist
                error(message('stats:probdists:CodegenArgsNotFromFitdist'));
            end
            coderobj = prob.coder.NormalDistribution(mlobj.mu, mlobj.sigma) ;
            
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
            if ~coderobj.IsTruncated % an untruncated distribution means it is the output of fitdist, not truncate()
                x = coderobj.InputData.Data;
                cens = coderobj.InputData.Censoring;
                freq = coderobj.InputData.Frequency;
                [nll, cov] = normlike(coderobj.ParameterValues,x,cens,freq);
                mlobj = prob.NormalDistribution.makeFitted(coderobj.ParameterValues ,nll,cov,...
                    x,cens,freq);
            else
                mlobj = prob.NormalDistribution(coderobj.ParameterValues(1),...
                    coderobj.ParameterValues(2));
            end
            if coderobj.IsTruncated
                mlobj = truncate(mlobj,coderobj.Truncation(1),coderobj.Truncation(2));
            end
        end
        
        function p = cdffunc(x,varargin)
            if nargin>1 && strcmpi(varargin{end},'upper')
                varargin(end) = [];
            else
                coder.internal.errorIf(nargin>1 && ischar(varargin{end})&&...
                    ~strcmpi(varargin{end},'upper'),...
                    'stats:cdf:UpperTailProblem');
                
                uflag=false;
            end
            p = localnormcdf(uflag,x,varargin{:});
            
            function p = localnormcdf(uflag,x,mu,sigma)
                coder.inline('always'); % inline this internal function
                z = (x-mu) ./ sigma;
                if uflag
                    z = -z;
                end
                
                p = coder.internal.nan(size(z),class(z));
                
                % Set edge case sigma=0
                if sigma ==0
                    if uflag
                        p(x<mu) = 1;
                        p(x>=mu) = 0;
                    else
                        p( x<mu) = 0;
                        p( x>=mu) = 1;
                    end
                else
                    % Normal cases
                    
                    % Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
                    % to produce accurate near-zero results for large negative x.
                    p  = 0.5 * erfc(-z ./ sqrt(2));
                end
                
            end
        end
        function y = pdffunc(x,mu,sigma)
            y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
        end
        function x = invfunc(p,mu,sigma,varargin)
            %             if nargin>3
            %                 x = norminv(p,mu,sigma,varargin{:});
            %                 return
            %             end
            x0 = -sqrt(2).*erfcinv(2*p);
            x = mu + sigma.*x0;
        end
    end % static, hidden methods
end % classdef
