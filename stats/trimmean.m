function m = trimmean(x,percent,flag,dim)
%TRIMMEAN Trimmed mean.
%   M = TRIMMEAN(X,PERCENT) calculates the trimmed mean of the values in X.
%   For a vector input, M is the mean of X, excluding the highest and
%   lowest K data values, where K=N*(PERCENT/100)/2 and where N is the
%   number of values in X.  For a matrix input, M is a row vector
%   containing the trimmed mean of each column of X.  For N-D arrays,
%   TRIMMEAN operates along the first non-singleton dimension.  PERCENT is
%   a scalar between 0 and 100.
%
%   M = TRIMMEAN(X,PERCENT,FLAG) controls how to trim when K is not an
%   integer.  FLAG can be chosen from the following:
%
%      'round'    Round K to the nearest integer (round to a smaller
%                 integer if K is a half integer).  This is the default.
%      'floor'    Round K down to the next smaller integer.
%      'weight'   If K=I+F where I is the integer part and F is the
%                 fraction, compute a weighted mean with weight (1-F) for
%                 the (I+1)th and (N-I)th values, and full weight for the
%                 values between them.
%
%   M = TRIMMEAN(X,PERCENT,FLAG,'all') is the trimmed mean of all the
%   elements of X.
%
%   M = TRIMMEAN(X,PERCENT,FLAG,DIM) takes the trimmed mean along dimension
%   DIM of X.
%
%   M = TRIMMEAN(X,PERCENT,FLAG,VECDIM) finds the trimmed mean of the 
%   elements of X based on the dimensions specified in the vector VECDIM.
%
%   The trimmed mean is a robust estimate of the sample location.
%
%   TRIMMEAN treats NaNs as missing values, and removes them.
%
%   See also MEAN, NANMEAN, IQR.

%   References:
%     [1] Wilcox, Rand R. "Introduction to Robust Estimation and
%         Hypothesis Testing." New York: Academic Press. 2005.
%     [2] Stigler, Stephen M. "Do robust estimators work with real data?"
%         Annals of Statistics, Vol. 5, No. 6, 1977, pp. 1055-1098.    

%   Copyright 1993-2018 The MathWorks, Inc.


if nargin > 2
    flag = convertStringsToChars(flag);
end

if nargin < 2
    error(message('stats:trimmean:TooFewInputs'));
elseif ~isscalar(percent) || percent >= 100 || percent < 0
    error(message('stats:trimmean:InvalidPercent'));
end

% Be flexible about syntax, so dim/flag may be in reverse order
if nargin<3
    flag = 'round';
    dim = [];
elseif nargin<4
    if isnumeric(flag) || strncmpi(flag,{'all'},max(strlength(flag), 1))
        dim = flag;
        flag = 'round';
    else
        dim = [];
    end
end
if isnumeric(flag)
    temp = dim;
    dim = flag;
    flag = temp;
end

dimSet = true;
if isempty(dim)
    % The output size for [] is a special case, handle it here.
    if isequal(x,[]), m = NaN('like',x); return; end

    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
    dimSet = false;
end

% Checking validity of dim
if dimSet
    if isnumeric(dim)
        if ~isreal(dim) || any(floor(dim) ~= ceil(dim)) || any(dim < 1) || any(~isfinite(dim))
            error(message('MATLAB:getdimarg:invalidDim'));
        end
        if ~isscalar(dim) && ~all(diff(sort(dim)))
            error(message('MATLAB:getdimarg:vecDimsMustBeUniquePositiveIntegers'));
        end
    elseif ((ischar(dim) && isrow(dim)) || ...
     (isstring(dim) && isscalar(dim) && (strlength(dim) > 0))) && ...
     strncmpi(dim,'all',max(strlength(dim), 1))
            x = x(:);
            dim = 1;
    else
        error(message('MATLAB:getdimarg:invalidDim'));
    end
end

if ischar(flag) && (isempty(flag) || isequal(lower(flag),'round'))
    F = @roundtrim;
elseif ischar(flag) && (isempty(flag) || isequal(lower(flag),'weighted'))
    F = @wtdtrim;
elseif ischar(flag) && isequal(lower(flag),'floor')
    F = @unwtdtrim;
else
    error(message('stats:trimmean:BadFlag'))
end

% Keep track of columns that were all missing data, or length zero.
allmissing = all(isnan(x),dim);

% Permute dimensions so we are working along columns

xdims = ndims(x);
if numel(dim)>1 || dim>1
    perm = [dim setdiff(1:max(xdims,max(dim)),dim)];       % New perm for 
    x = permute(x, perm);                                  % handling dim>xdims
    dimArgGiven = true;
else
    perm = [];
    dimArgGiven = false;
end
sz = size(x);
% Sort each column, get desired output size before inverse permutation
if dimArgGiven
    work_dim = 1:numel(dim);
else
    work_dim = dim;
end

if ~isempty(x)
    nrows = prod(sz(work_dim));             % Converting N-Dimensional Array
    ncols = numel(x) ./ nrows;              % into 2-D and then so that 
    x = reshape(x, nrows, ncols);           % operand dimensions come in 1st dim
end

sizey = size(x);
sizey(work_dim) = 1;   % Handling of multiple zero dimension array

x = sort(x,1);
if ~any(isnan(x(:)))
    % No missing data, operate on all columns at once
    if isempty(x)            % When one of the dimension of x is 0
        n = 0;
    else
        n = size(x,1);
    end
    [m,alltrimmed] = F(x,n,percent,sizey);
else
    % Need to loop over columns
    m = NaN(sizey,'like',x);
    alltrimmed = false(sizey);
    for j = 1:prod(sizey(2:end))
        n = find(~isnan(x(:,j)),1,'last');
        [m(j),alltrimmed(j)] = F(x(:,j),n,percent,[1 1]);
    end
end

% Permute back
szout = sz;
szout(work_dim) = 1;        %Using new dimension i.e work_dim

m = reshape(m,szout);
if length(alltrimmed) > 1 && all(szout>0)
    alltrimmed = reshape(alltrimmed,szout);
end
if ~isempty(perm)
    m = ipermute(m,perm);
    alltrimmed = ipermute(alltrimmed,perm);
end

% Warn if everything was trimmed, but not if all missing to begin with.
alltrimmed = (alltrimmed & ~allmissing);
if any(alltrimmed(:))
    if all(alltrimmed(:))
        warning(message('stats:trimmean:NoDataRemaining'));
    else
        warning(message('stats:trimmean:NoDataRemainingSomeColumns'));
    end
end

% --- Trim complete observations only, no weighting, rounding
function [m,alltrimmed] = roundtrim(x,n,percent,sizey)
k = n*percent/200;
k0 = round(k - eps(k));
if ~isempty(n) && n>0 && k0<n/2
    m = mean(x((k0+1):(n-k0),:),1);
    alltrimmed = false;
else
    m = NaN(sizey,'like',x);
    alltrimmed = true;
end

% --- Trim complete observations only, no weighting
function [m,alltrimmed] = unwtdtrim(x,n,percent,sizey)
k0 = floor(n*percent/200);
if ~isempty(n) && n>0 && k0<n/2
    m = mean(x((k0+1):(n-k0),:),1);
    alltrimmed = false;
else
    m = NaN(sizey,'like',x);
    alltrimmed = true;
end

% --- Weight observations to achieve desired percentage trimming
function [m,alltrimmed] = wtdtrim(x,n,percent,sizey)
k = n*percent/200;     % desired k to trim
k0 = floor(k);         % integer version
f = 1+k0-k;            % fraction to use for weighting
if ~isempty(n) && n>0 && (k0<n/2 || f>0)
    m = (sum(x((k0+2):(n-k0-1),:),1) + f*x(k0+1,:) + f*x(n-k0,:)) ...
        / (max(0,n-2*k0-2) + 2*f);
    alltrimmed = false;
else
    m = NaN(sizey,'like',x);
    alltrimmed = true;
end
