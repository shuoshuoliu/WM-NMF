function [h,p,ci,stats] = vartest(x,v,varargin)
%VARTEST  One-sample test of variance.
%  H = VARTEST(X,V) performs a chi-square test of the hypothesis that the
%  data in the vector X come from a normal distribution with variance V,
%  against the alternative that X comes from a normal distribution with a
%  different variance.  The result is H=0 if the null hypothesis ("variance
%  is V") cannot be rejected at the 5% significance level, or H=1 if the
%  null hypothesis can be rejected at the 5% level.
%
%  X may also be a matrix or an N-D array.  For matrices, VARTEST performs
%  separate tests along each column of X, and returns a vector of
%  results.  For N-D arrays, VARTEST works along the first non-singleton
%  dimension of X.  V must be a scalar.
%
%  VARTEST treats NaNs as missing values, and ignores them.
%
%  [H,P] = VARTEST(...) returns the p-value, i.e., the probability of
%  observing the given result, or one more extreme, by chance if the null
%  hypothesis is true.  Small values of P cast doubt on the validity of
%  the null hypothesis.
%
%  [H,P,CI] = VARTEST(...) returns a 100*(1-ALPHA)% confidence interval for
%  the true variance.
%
%  [H,P,CI,STATS] = VARTEST(...) returns a structure with the following
%  fields:
%     'chisqstat' -- the value of the test statistic
%     'df'        -- the degrees of freedom of the test
%
%  [...] = VARTEST(X,V,'PARAM1',val1,'PARAM2',val2,...) specifies one or
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
%          'both'  -- "variance is not V" (two-tailed test, default)
%          'right' -- "variance is greater than V" (right-tailed test)
%          'left'  -- "variance is less than V" (left-tailed test)
%
%  Example:  Is the standard deviation significantly different from 7?
%      load carsmall
%      [h,p,ci] = vartest(MPG, 7^2)
%
%  See also TTEST, ZTEST, VARTEST2.

%   Copyright 2005-2015 The MathWorks, Inc.

if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

narginchk(2,Inf);

if ~isscalar(v) || ~isnumeric(v) || ~isreal(v) || v<0
    error(message('stats:vartest:BadVar'));
end

% Process remaining arguments
alpha = 0.05;
tail = 0;    % code for two-sided;
dim = '';

if nargin>=3
    if isnumeric(varargin{1})
        % Old syntax
        %    VARTEST(X,V,ALPHA,TAIL,DIM)
        alpha = varargin{1};
        if nargin>=4
            tail = varargin{2};
              if nargin>=5
                  dim = varargin{3};
              end
        end
    elseif nargin==3
            error(message('stats:vartest:BadAlpha'))
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
elseif ~isscalar(alpha) || ~isnumeric(alpha)|| alpha <= 0 || alpha >= 1
    error(message('stats:vartest:BadAlpha'));
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
    % Figure out which dimension mean will work along
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end
if ~isscalar(dim) || ~ismember(dim,1:ndims(x))
    error(message('stats:vartest:BadDim', ndims( x )));
end        

nans = isnan(x);
dims = ndims(x);
x(nans) = 0;
if any(nans(:))
    samplesize = sum(~nans,dim);
else
    samplesize = size(x,dim); % a scalar, => a scalar call to tinv
end

df = max(samplesize - 1,0);
xmean = sum(x,dim) ./ max(1,samplesize);
if isscalar(xmean)
   xcntr = x - xmean;
else
   rep = ones(1,dims);
   rep(dim) = size(x,dim);
   xcntr = x - repmat(xmean,rep);
end
xcntr(nans) = 0;
sumsq = sum(abs(xcntr).^2,dim);
if v>0
   chisqstat = sumsq ./ v;
else
   chisqstat = Inf( size(sumsq), 'like', sumsq );
   chisqstat(sumsq==0) = NaN;
end

% Compute the correct p-value for the test, and confidence intervals
% if requested.
if tail == 0 % two-tailed test
    p = chi2cdf(chisqstat, df);
    p = 2*min(p, 1-p);
    if nargout > 2
        ci = cat(dim, sumsq ./ chi2inv(1 - alpha/2, df), ...
                      sumsq ./ chi2inv(alpha/2, df));
    end
elseif tail == 1 % right one-tailed test
    p = chi2pval(chisqstat, df);
    if nargout > 2
        ci = cat(dim, sumsq./chi2inv(1 - alpha, df), Inf(size(p)));
    end
elseif tail == -1 % left one-tailed test
    p = chi2cdf(chisqstat, df);
    if nargout > 2
        ci = cat(dim, zeros(size(p)), sumsq./chi2inv(alpha, df));
    end
end
  
% Determine if the actual significance exceeds the desired significance
h = cast(p <= alpha, 'like', p);
h(isnan(p)) = NaN; % p==NaN => neither <= alpha nor > alpha

if nargout >= 4
    stats = struct('chisqstat', chisqstat, 'df', cast(df,'like',chisqstat));
    if isscalar(df) && ~isscalar(chisqstat)
        stats.df = repmat(stats.df,size(chisqstat));
    end
end

