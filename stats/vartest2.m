function [h,p,ci,stats] = vartest2(x,y,varargin)
%VARTEST2 Two-sample F test for equal variances.
%   H = VARTEST2(X,Y) performs an F test of the hypothesis that two
%   independent samples, in the vectors X and Y, come from normal
%   distributions with the same variance, against the alternative that
%   they come from normal distributions with different variances.
%   The result is H=0 if the null hypothesis ("variances are equal")
%   cannot be rejected at the 5% significance level, or H=1 if the null
%   hypothesis can be rejected at the 5% level.  X and Y can have
%   different lengths.
%
%   X and Y can also be matrices or N-D arrays.  For matrices, VARTEST2
%   performs separate tests along each column, and returns a vector of
%   results.  X and Y must have the same number of columns.  For N-D
%   arrays, VARTEST2 works along the first non-singleton dimension.  X and Y
%   must have the same size along all the remaining dimensions.
%
%   VARTEST2 treats NaNs as missing values, and ignores them.
%
%   [H,P] = VARTEST2(...) returns the p-value, i.e., the probability of
%   observing the given result, or one more extreme, by chance if the null
%   hypothesis is true.  Small values of P cast doubt on the validity of
%   the null hypothesis.
%
%   [H,P,CI] = VARTEST2(...) returns a 100*(1-ALPHA)% confidence interval for
%   the true ratio var(X)/var(Y).
%
%   [H,P,CI,STATS] = VARTEST2(...) returns a structure with the following
%   fields:
%      'fstat'  -- the value of the test statistic
%      'df1'    -- the numerator degrees of freedom of the test
%      'df2'    -- the denominator degrees of freedom of the test
%
%  [...] = VARTEST2(X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies one or
%  more of the following name/value pairs:
%
%       Parameter       Value
%       'alpha'         A value ALPHA between 0 and 1 specifying the
%                       significance level as (100*ALPHA)%. Default is
%                       0.05 for 5% significance.
%       'dim'           Dimension DIM to work along. For example, specifying
%                       'dim' as 1 tests the column variances. Default is the
%                       first non-singleton dimension.
%       'tail'          A string specifying the alternative hypothesis:
%           'both'  -- "variances are not equal" (two-tailed test)
%           'right' -- "variance of X is greater than variance of Y"
%                      (right-tailed test)
%           'left'  -- "variance of X is less than variance of Y" (left-
%                      tailed test)
%   [...] = VARTEST2(X,Y,ALPHA,TAIL,DIM) works along dimension DIM
%   of X and Y.
%
%   Example:  Is the variance significantly different for two model years,
%             and what is a confidence interval for the ratio of these
%             variances?
%      load carsmall
%      [h,p,ci] = vartest2(MPG(Model_Year==82),MPG(Model_Year==76))
%
%   See also ANSARIBRADLEY, VARTEST, VARTESTN, TTEST2.

%   Copyright 2005-2015 The MathWorks, Inc.


if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

narginchk(2,Inf);

% Process remaining arguments
alpha = 0.05;
tail = 0;    % code for two-sided;
dim = '';

if nargin>=3
    if isnumeric(varargin{1})
        % Old syntax
        %    VARTEST2(X,Y,ALPHA,TAIL,DIM)
        alpha = varargin{1};
        if nargin>=4
            tail = varargin{2};
              if nargin>=5
                  dim = varargin{3};
              end
        end
        
    elseif nargin==3
            error(message('stats:vartest2:BadAlpha'));
            
    else
        % Calling sequence with named arguments
        okargs =   {'alpha' 'tail' 'dim'};
        defaults = {0.05    'both'  ''};
        [alpha, tail, dim] = ...
                         internal.stats.parseArgs(okargs,defaults,varargin{:});
    end
end

if isempty(alpha)
    alpha = 0.05;
elseif ~isscalar(alpha) || alpha <= 0 || alpha >= 1
    error(message('stats:vartest2:BadAlpha'));
end

if isempty(tail)
    tail = 0;
elseif isnumeric(tail) && isscalar(tail) && ismember(tail,[-1 0 1])
    % OK, grandfathered
else
    [~,tail] = internal.stats.getParamVal(tail,{'left','both','right'},'''tail''');
    tail = tail - 2;
end

if isempty(dim)
    % Figure out which dimension to work along by looking at x.
    % y will have to be compatible.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
    
    % If we haven't been given an explicit dimension, and we have two
    % vectors, then make y the same orientation as x.
    if isvector(x) && isvector(y)
        if (size(x,1)==1)
            y = y(:)'; 
        else
            y = y(:);
        end
    end
end
if ~isscalar(dim) || ~ismember(dim,1:ndims(x))
    error(message('stats:vartest2:BadDim', ndims( x )));
end        

% Make sure all of x's and y's non-working dimensions are identical.
sizex = size(x); sizex(dim) = 1;
sizey = size(y); sizey(dim) = 1;
if ~isequal(sizex,sizey)
    error(message('stats:vartest2:InputSizeMismatch'));
end

% Compute statistics for each sample
[df1,varx] = getstats(x,dim);
[df2,vary] = getstats(y,dim);

% Compute F statistic
F = NaN(size(varx),'like',internal.stats.dominantType(varx,vary));
t1 = (vary>0);
F(t1) = varx(t1) ./ vary(t1);
t2 = (varx>0) & ~t1;
F(t2) = Inf;

% Compute the correct p-value for the test, and confidence intervals
% if requested.
if tail == 0 % two-tailed test
    p = 2*min(fcdf(F,df1,df2),fpval(F,df1,df2));
    if nargout > 2
        % Avoid precision loss from subtracting alpha from one
        ci = cat(dim, F.*finv(alpha/2,df2,df1), ... % == F./finv(1-alpha/2,df1,df2)
                      F./finv(alpha/2,df1,df2));
    end
elseif tail == 1 % right one-tailed test
    p = fpval(F,df1,df2);
    if nargout > 2
        ci = cat(dim, F.*finv(alpha,df2,df1), ... % == F./finv(1-alpha,df1,df2)
                      Inf(size(F)));
    end
elseif tail == -1 % left one-tailed test
    p = fcdf(F,df1,df2);
    if nargout > 2
        ci = cat(dim, zeros(size(F)), ...
                      F./finv(alpha,df1,df2));
    end
end

% Determine if the actual significance exceeds the desired significance
h = cast(p <= alpha, 'like', p);
h(isnan(p)) = NaN; % p==NaN => neither <= alpha nor > alpha

if nargout > 3
    stats = struct('fstat', F, 'df1', cast(df1,'like',F), ...
                               'df2', cast(df2,'like',F));
    if isscalar(df1) && ~isscalar(F)
        stats.df1 = repmat(stats.df1,size(F));
    end
    if isscalar(df2) && ~isscalar(F)
        stats.df2 = repmat(stats.df2,size(F));
    end
end

% -----------------------------------
function [df,varx] = getstats(x,dim)
%GETSTATS Compute statistics for one sample
  
% Get sample sizes and df
xnans = isnan(x);
nx = sum(~xnans,dim);
df = max(nx-1, 0);

% Get means
x(xnans) = 0;
xmean = sum(x,dim) ./ max(1,nx);

% Get variances
if isscalar(xmean)
   xcntr = x - xmean;
else
   rep = ones(1,ndims(x));
   rep(dim) = size(x,dim);
   xcntr = x - repmat(xmean,rep);
end
xcntr(xnans) = 0;
varx = sum(abs(xcntr).^2,dim);
t = (df>0);
varx(t) = varx(t) ./ df(t);
varx(~t) = NaN;

% Make df a scalar if possible, better for later finv call
if numel(df)>1 && all(df(:)==df(1))
   df = df(1);
end
