function [h,p,ci,zval] = ztest(x,m,sigma,varargin)
%ZTEST  One-sample Z-test.
%   H = ZTEST(X,M,SIGMA) performs a Z-test of the hypothesis that the data
%   in the vector X come from a distribution with mean M, and returns the
%   result of the test in H.  H=0 indicates that the null hypothesis
%   ("mean is M") cannot be rejected at the 5% significance level.  H=1
%   indicates that the null hypothesis can be rejected at the 5% level.  The
%   data are assumed to come from a normal distribution with standard
%   deviation SIGMA.
%
%   X may also be a matrix or an N-D array.  For matrices, ZTEST performs
%   separate Z-tests along each column of X, and returns a vector of
%   results.  For N-D arrays, ZTEST works along the first non-singleton
%   dimension of X.  M and SIGMA must be scalars.
%
%   ZTEST treats NaNs as missing values, and ignores them.
%
%   [H,P] = ZTEST(...) returns the p-value, i.e., the probability of
%   observing the given result, or one more extreme, by chance if the null
%   hypothesis is true.  Small values of P cast doubt on the validity of
%   the null hypothesis.
%
%   [H,P,CI] = ZTEST(...) returns a 100*(1-ALPHA)% confidence interval for
%   the true mean.
%
%   [H,P,CI,ZVAL] = ZTEST(...) returns the value of the test statistic.
%
%   [...] = ZTEST(X,M,SIGMA,'PARAM1',val1,'PARAM2',val2,...) specifies one
%   or more of the following name/value pairs:
%
%       Parameter       Value
%       'alpha'         A value ALPHA between 0 and 1 specifying the
%                       significance level as (100*ALPHA)%. Default is
%                       0.05 for 5% significance.
%       'dim'           Dimension DIM to work along. For example, specifying
%                       'dim' as 1 tests the column means. Default is the
%                       first non-singleton dimension.
%       'tail'          A string specifying the alternative hypothesis:
%           'both'  -- "mean is not M" (two-tailed test)
%           'right' -- "mean is greater than M" (right-tailed test)
%           'left'  -- "mean is less than M" (left-tailed test)
%
%   See also TTEST, SIGNTEST, SIGNRANK, VARTEST.

%   References:
%      [1] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, page 206.

%   Copyright 1993-2017 The MathWorks, Inc.

if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

narginchk(3,Inf);
if ~isscalar(m)
    error(message('stats:ztest:NonScalarM'));
elseif ~isscalar(sigma) || (sigma < 0)
    error(message('stats:ztest:NonScalarSigma'));
end

% Process remaining arguments
alpha = 0.05;
tail = 0;    % code for two-sided;
dim = '';

if nargin>=4
    if isnumeric(varargin{1})
        % Old syntax
        %   ZTEST(X,M,SIGMA,ALPHA,TAIL,DIM)
        alpha = varargin{1};
        if nargin>=5
            tail = varargin{2};
              if nargin>=6
                  dim = varargin{3};
              end
        end
        
    elseif nargin==4
            error(message('stats:ztest:BadAlpha'));
 
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
    error(message('stats:ztest:BadAlpha'));
end

if isempty(tail)
    tail = 0;
elseif isnumeric(tail) && isscalar(tail) && ismember(tail,[-1 0 1])
    % OK, grandfathered
else
    [~,tail] = internal.stats.getParamVal(tail,{'left','both','right'},'TAIL'); tail = tail - 2;
end

if isempty(dim)
    % Figure out which dimension mean will work along
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

nans = isnan(x);
if any(nans(:))
    samplesize = sum(~nans,dim);
else
    samplesize = size(x,dim); % make this a scalar if possible
end
xmean = nanmean(x,dim);
ser = sigma ./ sqrt(samplesize);
zval = (xmean - m) ./ ser;

% Compute the correct p-value for the test, and confidence intervals
% if requested.
if tail == 0 % two-tailed test
    p = 2 * normcdf(-abs(zval),0,1);
    if nargout > 2
        crit = norminv(1 - alpha/2, 0, 1) .* ser;
        ci = cat(dim, xmean-crit, xmean+crit);
    end
elseif tail == 1 % right one-tailed test
    p = normcdf(-zval,0,1);
    if nargout > 2
        crit = norminv(1 - alpha, 0, 1) .* ser;
        ci = cat(dim, xmean-crit, Inf(size(p)));
    end
elseif tail == -1 % left one-tailed test
    p = normcdf(zval,0,1);
    if nargout > 2
        crit = norminv(1 - alpha, 0, 1) .* ser;
        ci = cat(dim, -Inf(size(p)), xmean+crit);
    end
end

% Determine if the actual significance exceeds the desired significance
h = cast(p <= alpha, 'like', p);
h(isnan(p)) = NaN; % p==NaN => neither <= alpha nor > alpha
