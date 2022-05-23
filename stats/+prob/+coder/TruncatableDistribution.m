classdef TruncatableDistribution %#codegen
%

%   Copyright 2019 The MathWorks, Inc.
    
    % Code generation version of the TruncatableDistribution
    % Interface for truncatable distributions.
    % In most cases of extending support to other distributions for fitdist,
    % it should not require any modification.
    
    % This classdef is significantly different from its MATLAB version.
    % This class borrows its properties and methods from the following
    % classes under the "prob" package:
    % 1. ProbabilityDistribution
    % 2. UnivariateDistribution
    % 3. TruncatableDistribution
    
    % 1. prob.ProbabilityDistribution methods

    methods
        function s = std(this)
            %STD Standard deviation of the distribution.
            %   S = STD(PD) returns the standard deviation S of the probability
            %   distribution PD.
            %
            %   See also VAR, MEAN.
            s = sqrt(var(this));
        end
    end
    
    % 2. prob.UnivariateDistribution methods
    methods
        function y = median(this)
            %MEDIAN Median.
            %    M = MEDIAN(P) computes the median M for the probability distribution
            %    P. M is the value such that half of the probabilty is below M and half
            %    is above.
            %
            %    See also IQR, ICDF.
            y = icdf(this,.5);
        end
        function y = iqr(this)
            %IQR Interquartile range.
            %    R = IQR(P) computes the interquartile range R for the probability
            %    distribution P. R is the difference between the 75th and 25th
            %    percentage points of P.
            %
            %    See also MEDIAN, ICDF.
            y = icdf(this,.75) - icdf(this,.25);
        end
    end
    
    % 3. prob.TruncatableDistribution properties and methods
    
    properties(GetAccess='public',SetAccess='protected')
        %Truncation - Two-element vector defining the truncation interval.
        %    Truncation is a two-element vector of the form [LOWER,UPPER]
        %    indicating that the distribution is truncated to the interval with
        %    lower limit LOWER and upper limit UPPER. For an untruncated
        %    distribution, the vector is [-Inf,Inf] theoretically (for
        %    codegen, it is initialized with all zeros for convenience). This information
        %    is used for computing the IsTruncated property.
        %  The TRUNCATE method sets the Truncation property.
        %
        %
        %    See also truncate, IsTruncated.
        
        Truncation
    end
    properties(Dependent=true,GetAccess='public',SetAccess='protected')
        %IsTruncated - Boolean indicating if the distribution is truncated.
        %    IsTruncated is true if the distribution is truncated to a subset of
        %    the real line, or false if it is not.
        %
        %    See also truncate, Truncation.
        IsTruncated
    end
    
    
    
    methods(Access=protected,Abstract = true)
        % The TruncatableDistribution class implements the methods below
        % and provides support for truncation. These implementations use a
        % template pattern that expects a derived class to implement a hook
        % or callback to perform the basic (untruncated) calculation:
        %    Method        Hook
        %    pdf           pdffun
        %    cdf           cdffun
        %    icdf          icdffun
        %    random        randomfun
        y = icdffun(this,x)
        y = cdffun(this,x)
        y = pdffun(this,x)
    end
    
    methods
        
        
        function y = pdf(this,xin)
            %PDF Probability density function.
            %    Y = PDF(P,X) computes the probability density function of the
            %    probability distribution P at the values in X, and returns the result
            %    in the array Y.
            %
            %    See also FITDIST.
            narginchk(2,2); % add a nargin check so that the error message matches matlab
            checkData(xin,'x','pdf');
            requireScalar(this);
            % the type of the distribution parameters determines the output of the final type
            if isa(xin,'double') && isa(this.ParameterValues,'single')
                x = cast(xin,'single'); % cast this to a single value if the fit was performed on singles
            else
                x = xin; % if both singles: output is single, in all other cases, output will be double
            end
            y = zeros(size(x),'like',x);
            
            if ~this.IsTruncated
                y = pdffun(this,x);
                for idx = 1:numel(x)
                    if isnan(x(idx)) % if x was a compile-time constant, we could eliminate this check from the code for non nan inputs
                        y(idx) = coder.internal.nan('like',x);
                    end
                end
            else
                [plower,pupper,Lower,Upper] = tailprobs(this);
                
                for idx = 1:numel(x)
                    if isnan(x(idx)) % if x was a compile-time constant, we could eliminate this check from the code for non nan inputs
                        y(idx) = coder.internal.nan('like',x);
                        continue;
                    end
                    
                    if (x(idx)>=Lower && x(idx) <= Upper)
                        y(idx) = pdffun(this,x(idx)) / (pupper-plower);
                    end
                end
            end
        end
        
        function y = cdf(this,xin,varargin)
            %CDF Cumulative distribution function.
            %    Y = CDF(P,X) computes the cumulative distribution function of the
            %    probability distribution P at the values in X, and returns the result
            %    in the array Y.
            %
            %    Y = CDF(P,X,'upper') computes the upper tail probability of the
            %    probability distribution P at the values in X, and returns the result
            %    in the array Y.
            %
            %    See also FITDIST.
            narginchk(2,3); % add a nargin check so that the error message matches matlab
            checkData(xin,'x','cdf');
            requireScalar(this);
            
            % the type of the distribution parameters determines the output of the final type
            if isa(xin,'double') && isa(this.ParameterValues,'single')
                x = cast(xin,'single'); % cast this to a single value if the fit was performed on singles
            else
                x = xin; % if both singles: output is single, in all other cases, output will be double
            end
            
            %             wasnan = isnan(x);
            y = coder.nullcopy(zeros(size(x),'like',x));
            
            
            if ~this.IsTruncated
                y = cdffun(this,x,varargin{:});
            else
                [plower,pupper,Lower,Upper] = tailprobs(this);
                
                % move this check out of the for loop if possible
                if ~isempty(varargin)
                    if strcmpi(varargin{end},'upper')
                        uflag = true;
                    else
                        uflag = false;
                    end
                    coder.internal.errorIf(~strcmpi(varargin{end},'upper'),...
                        'stats:cdf:UpperTailProblem');
                    
                else
                    uflag = false;
                end
                
                for idx = 1:numel(x)
                    if isnan(x(idx)) % if x was a compile-time constant, we could eliminate this check from the code for non nan inputs
                        y(idx) = coder.internal.nan;
                        continue;
                    end
                    x(idx) = min(x(idx),Upper);
                    t = (x(idx)>=Lower & x(idx)<=Upper);
                    t1 = (x(idx)<Lower);
                    
                    
                    
                    if uflag
                        if (t||t1)
                            y(idx) = (pupper-cdffun(this,x(idx))) / (pupper-plower);
                            y(idx) = 1;
                        end
                    else
                        y(idx) = (cdffun(this,x(idx))-plower) / (pupper-plower);
                    end
                end
                
            end
            
            
            
            
        end
        
        function y = icdf(this,xin)
            %ICDF Inverse cumulative distribution function.
            %    Y = ICDF(P,PROB) computes the inverse cumulative distribution function
            %    of the probability distribution P at the values in PROB, and returns
            %    the result in the array Y.
            %
            %
            %    See also FITDIST.
            narginchk(2,2);
            checkData(xin,'p','icdf');
            
            % the type of the distribution parameters determines the output of the final type
            if isa(xin,'double') && isa(this.ParameterValues,'single')
                pin = cast(xin,'single'); % cast this to a single value if the fit was performed on singles
            else
                pin = xin; % if both singles: output is single, in all other cases, output will be double
            end
            
            requireScalar(this);
            
            if this.IsTruncated
                % Make p relative to the truncation limit cdf values
                [plower,pupper,Lower,Upper] = tailprobs(this);
                p = plower + pin*(pupper-plower);
            else
                Lower = cast(0,'like',this.ParameterValues);
                Upper = cast(0, 'like',this.ParameterValues);
                p = pin;
            end
            
            y = icdffun(this,p);
            
            % All distributions from here are assumed continuous
            if this.IsTruncated
                SL = max(Lower,icdffun(this,0));
                SU = min(Upper,icdffun(this,1));
                for idx=1:numel(y)
                    if any(y(idx)<SL)
                        y(idx) = SL;
                    end
                    if any(y(idx)>SU)
                        y(idx) = SU;
                    end
                end
            end
            
            
            for idx=1:numel(pin)
                if pin(idx)<0 || pin(idx)>1
                    y(idx) = coder.internal.nan('like',pin);
                end
            end
            
            
        end
        
        
        function td = truncate(this,lower,upper)
            %TRUNCATE Truncate probability distribution to an interval.
            %    T = TRUNCATE(P,LOWER,UPPER) takes a probability distribution P and
            %    returns another probability distribution T that represents P
            %    truncated to the interval with lower limit LOWER and upper limit
            %    UPPER. The pdf of T is zero outside the interval. Inside the
            %    interval it is equal to the pdf of P, but divided by the probability
            %    assigned to that interval by P.
            %
            %    See also Truncation, IsTruncated.
            requireScalar(this);
            checkTruncationArgs(this,lower,upper);
            
            td = this; % is it being copied?
            td.Truncation = cast([lower upper], 'like', this.ParameterValues);
        end
        
        function tf = get.IsTruncated(this)
            tr = this.Truncation;
            tf = ~(nnz(tr)==0 || (numel(tr)==2 && tr(1)==-coder.internal.inf && tr(2)==coder.internal.inf));
        end
        
    end
    
    
    
    
    
    
    methods (Hidden)
        % Hidden but public so interface classes can use this utility
        function requireScalar(this)
            % essentially a validation function for the first input to methods, same error as MATLAB
            coder.internal.errorIf(~isscalar(this), 'stats:probdists:RequiresScalar');
            % do further type validation if indeed scalar, only independent, non-constant properties need to be
            % validated. The idea behind it is to restrict options available to adventurous customers, through coder.typeof, in some ways but not all.
            
            validateattributes(this.ParameterValues,{'double', 'single'},{'size', [1 this.NumParameters]},mfilename,'ParameterValues');
            validateattributes(this.Truncation,{'double', 'single'},{'size', [1 2]},mfilename,'Truncation');
            validateattributes(this.InputData, {'struct'},{'scalar'} ,mfilename, 'InputData');
        end
        
    end
    
    methods(Access=protected,Hidden=true)
        
        function TM = truncatedMoment(this,degree)
            %TRUNCATEDMOMENT Calculates the truncated central moment of degree DEGREE.
            %
            % This function should be called for continuous functions distributions only.
            %
            % TM           Mean or variance of the truncated distributiion, depending
            %              on DEGREE.
            %
            % degree       The degree of the moment. Must be 1 (mean) or 2 (variance).
            %
            % See also quadgk.
            
            
            % Create waypoints for the numerical integration.
            % Use equi-spaced quantiles.
            waypoints = icdffun(this,(0.01:0.01:0.99)');
            varargin = {'Waypoints',waypoints};
            
            % The following quantities are:
            % TL     Lower truncation limit
            % TU     Upper truncation limit
            % SL     Lower extreme of support, untruncated distribution
            % SU     Upper extreme of support, untruncated distribution
            % L      Endpoint for semi-infinite integrals
            % U      Endpoint for semi-infinite integrals
            
            [~,~,TL,TU] = tailprobs(this);
            SL = icdffun(this,0);
            SU = icdffun(this,1);
            L = max(TL,SL);
            U = min(TU,SU);
            
            if isequal(degree,1)
                f = @(x) x .* pdf(this,x);
            else
                theMean = truncatedMoment(this,1);
                f = @(x) (x-theMean).^2 .* pdf(this,x);
            end
            
            if isequal(SU,coder.internal.inf)
                % Infinite right support (untruncated). Upper truncation
                % point is unbounded, and (L,U) may be very large,though
                % finite.  The integral() function may fail on such
                % intervals, so give it a semi-infinite integral.
                [TM] = doIntegral(f,L,coder.internal.inf,varargin{:});
            elseif isequal(SL,-coder.internal.inf)
                % Likewise for infinite left tail.
                [TM] = doIntegral(f,-coder.internal.inf,U,varargin{:});
            else
                % Untruncated support is finite
                [TM] = doIntegral(f,L,U,varargin{:});
            end
        end
        
        function checkTruncationArgs(this,lower,upper)
            coder.internal.errorIf(~(isscalar(lower) && isnumeric(lower) && isreal(lower) && ~isnan(lower)), ...
                'stats:probdists:BadTruncationParameter','LOWER');
            coder.internal.errorIf(~(isscalar(upper) && isnumeric(upper) && isreal(upper) && ~isnan(upper)),...
                'stats:probdists:BadTruncationParameter','UPPER');
            coder.internal.errorIf(~(lower<upper),...
                'stats:probdists:LowerLTUpper','LOWER','UPPER');
            
            lo = protectiveLowBound(lower);
            plower = cdffun(this,lo);
            pupper = cdffun(this,upper);
            coder.internal.errorIf(isequal(plower,pupper), ...
                'stats:probdists:ZeroMassTruncation');
        end
        
        
        function [plower,pupper,Lower,Upper] = tailprobs(this)
            % Get tail probabilities and truncation bounds. Take pains to
            % compute the cdf just below the lower truncation limit, in
            % case this distribution gives positive probability exactly at
            % that point.
            if ~this.IsTruncated
                Lower = -coder.internal.inf('like', this.ParameterValues);
                Upper = coder.internal.inf('like', this.ParameterValues);
                lo = -coder.internal.inf('like', this.ParameterValues);
            else
                Lower = this.Truncation(1);
                Upper = this.Truncation(2);
                lo = protectiveLowBound(Lower);
            end
            plower = cdffun(this,lo);
            pupper = cdffun(this,Upper);
        end
    end
    
    
end % classdef

function checkData(x, methodinput, methodname) % ensure the inputs are either scalars or array of scalars
coder.internal.errorIf(~(isnumeric(x) && isreal(x)),... % && ismatrix(x)
    'stats:probdists:BadInputsToMethods',methodinput, methodname);
end

function lo = protectiveLowBound(low)
if isfinite(low)
    lo = low - eps(low);
    if lo+eps(lo) < low
        lo = lo + eps(lo);
    end
else
    lo = low;
end
end


function [q] = doIntegral(f,L,U,varargin)
q = quadgk(f,L,U,varargin{:}); % integral is not supported for codegen, use quadgk instead
end



