function D = heteropdist2(X,Y,dist,varargin)
% HETEROPDIST2 Pairwise distance between two sets of heterogeneous data.
%
%   INPUTS AND OUTPUTS
%   D = HETEROPDIST2(X,Y,DIST) returns a matrix D containing the distances
%   specified by DIST, between each pair of observations in the MX-by-N
%   data matrix X and MY-by-N data matrix Y. Rows of X and Y correspond to
%   observations, and columns correspond to variables. D is an MX-by-MY
%   matrix, with the (I,J) entry equal to distance between observation I in
%   X and observation J in Y.
%
%   For DIST, the following choices exist-
%   +------------+------------------------------+
%   | DIST       | DESCRIPTION                  |
%   +------------+------------------------------+
%   | 'ofd'      | Occurrence Frequency Distance|
%   +------------+------------------------------+
%   | 'goodall3' | Modified Goodall distance    |
%   +------------+------------------------------+
%
%   Categorical predictors are expected to be encoded as integers for this
%   function. This function uses accumarray for speed. To work well with
%   accumarray, it is recommended that categorical predictors are encoded
%   as consecutive integers using ismember() just like classreg. For
%   example, if there are only two categories- 1 and 100, it is advisable
%   to encode them as 1 and 2.
%   The optional argument at the fourth place is used to specify
%   the categorical predictor indices.  The indices themselves can be
%   specified either as integers OR logicals.
%   If not specified, all predictors are treated as categorical.
%
%   ALGORITHM DETAILS
%   All supported distances are data-driven metrics and require frequency
%   computation for categorical predictors. For continuous predictors, an
%   ECDF is calculated. As a convention, the first input is used as the
%   “reference” for frequency and ECDF computation. It might so happen that
%   the categorical values in Y do not appear in X for the same feature. To
%   avoid such zero probabilities, an “add-one” Laplace correction is added
%   for frequency calculation. For continuous features in Y, ECDF value is
%   taken as the average of values in X. For example, if X has feature A
%   with values 1 and 3, and Y has feature A with value 2. ECDF(Y(A) == 2)
%   = (ECDF(X(A)==1) + ECDF(X(A)==3))/2. Values in Y outside the range of
%   values in X are assigned a 0 or 1 value depending on where they lie.
%
%   Per feature distance contribution is calculated in the following way
%   for categorical predictors (formulas are listed in the reference).
%   1. OfD (Occurrence Frequency Distance) For a match, occurrence
%   frequency assigns zero distance. For mismatches on less frequent
%   values, 'ofd' assigns a higher distance and mismatches on more frequent
%   values are assigned a lower distance.
%
%   2. Goodall3 (Modified Goodall distance)
%   Defined in the reference below, the Goodall3 measure is a variant of
%   the "Goodall" distance which "assigns a high similarity if the matching
%   values are infrequent regardless of the frequencies of the other
%   values." A high similarity implies a small distance and reverse is also
%   true. For mismatches, the unweighted distance contribution of the
%   predictor is 1.
%
%   The contribution of each predictor to the total distance is weighted
%   uniformly by (1/TotalNumberOfPredictors).
%
%   EXAMPLE USAGE:
%   X = randi(10,10,3);
%   Y = randi(10,10,3);
%   Dof = internal.stats.heteropdist2(X, Y, 'ofd', [1 0 1]);
%   Dgood = internal.stats.heteropdist2(X, Y, 'goodall3', [1 0 1]);
%
%   Reference: Boriah, Shyam, Varun Chandola, and Vipin Kumar. "Similarity
%   measures for categorical data: A comparative evaluation." Proceedings
%   of the 2008 SIAM international conference on data mining. Society for
%   Industrial and Applied Mathematics, 2008.

%   Copyright 2019-2020 The MathWorks, Inc.

if nargin < 2
    error(message('stats:pdist2:TooFewInputs'));
end

if ~ismatrix(X) || ~ismatrix(Y)
    error(message('stats:pdist2:UnsupportedND'));
end

[nx,p] = size(X);
[ny,py] = size(Y);

% integer/logical/char/anything data is converted to float.

try
    outClass = superiorfloat(X,Y);
catch
    if isfloat(X)
        outClass = class(X);
    elseif isfloat(Y)
        outClass = class(Y);
    else
        outClass = 'double';
    end
end

if py ~= p
    error(message('stats:pdist2:SizeMismatch'));
end

% Degenerate case, just return an empty of the proper size.
if (nx == 0) || (ny == 0)
    D = zeros(nx,ny,outClass); % X and Y were single/double, or cast to double
    return;
end

methods = {'OfD';'Goodall3'}; % OF metric elongated to OfD because the C++ code expects 3 letter metrics
% Lin is not included here because a weighting scheme for continuous
% distances that can be combined with Lin remains to be found

i = find(strncmpi(dist,methods,length(dist)));
if length(i) > 1
    error(message('stats:pdist2:AmbiguousDistance', dist));
elseif isempty(i)
    error(message('stats:pdist2:UnrecognizedDistance', dist));
else
    dist = lower(methods{i}(1:3));
end

% process the optional argument of categorical predictor indices
if ~isempty(varargin)
    categCols = varargin{1};
    if islogical(categCols)
        categCols = find(categCols);
    end
else
    categCols = 1:p;
end
assert(numel(categCols)>0); % guard against empty input from lime()

if ~strcmp(class(X),outClass) || ~strcmp(class(Y),outClass)
    warning(message('stats:pdist2:DataConversion', outClass));
end

% add the same casting as pdist2
X = cast(X,outClass);
Y = cast(Y,outClass);

if  ~isreal(X) || ~isreal(Y)
    error(message('stats:pdist2:ComplexData'));
end

freqAndECDFX = zeros(size(X),'like',X); % getting the types right
freqAndECDFY = zeros(size(Y), 'like', Y); % getting the types right

for ind = 1:size(X,2) % for every column (in X and Y)
    if any(ind==categCols(:)) % if it is a categorical column, calculate frequency
        % categorical columns are encoded as integers so
        % accumarray can provide the counts
        uniqueX = unique(X(:,ind)); % unique categories in X
        F = accumarray(X(:,ind),1); % frequency for integers in X, essentially a map of values to their counts (starting with 1)
        freqAndECDFX(:,ind) = F(X(:,ind));
        isAMatchOfYinX = ismember(Y(:,ind),uniqueX);
        freqAndECDFY(isAMatchOfYinX,ind) = F(Y(isAMatchOfYinX,ind));
        freqAndECDFX(:,ind) =  freqAndECDFX(:,ind) + 1; % laplace correction
        freqAndECDFY(:, ind) = freqAndECDFY(:, ind) + 1; % laplace correction
    else
        [~,~,idxx] = unique(X(:,ind)); % indices used to map X values to their ECDF value
        [F,xsorted] = ecdf((X(:,ind))); % xsorted is a vector of sorted, unique values in x (the first value is repeated)
        freqAndECDFX(:,ind) = F(idxx(1:size(X,1))+1); % an offset of 1 because the first two x values are repeated in ecdf(x)
        
        % for values in Y we use the ECDF for X
        % if a value in Y lies between two values in X, its ECDF value is
        % the average of the two,
        % if a values in Y lies outside of the range of values in X,
        y = Y(:,ind);
        minX = min(X(:,ind));
        maxX = max(X(:,ind));
        yValuesLessThanMin = (y<minX);
        yValuesGreaterThanMax = (y>maxX);
        yValuesInXRange = ~yValuesLessThanMin & ~yValuesGreaterThanMax;
        freqAndECDFY(yValuesGreaterThanMax, ind) = ones(sum(yValuesGreaterThanMax),1); % set this ECDF to 1
        
        valuesToBeFilled = sum(yValuesInXRange);
        low_indices = zeros(valuesToBeFilled,1);
        high_indices = zeros(valuesToBeFilled, 1);
        y_fill = y(yValuesInXRange);
        for k = 1:valuesToBeFilled % use binary search to find the location of the point in Y in log(N) time
            [low_indices(k), high_indices(k)] = binarySearchIdx(1, numel(xsorted(2:end)), xsorted(2:end), y_fill(k)); % ECDF offset of 1, use x from the second index
        end
        freqAndECDFY(yValuesInXRange,ind) = (F(low_indices+1) + F(high_indices+1))/2;  % take the average of the ECDF (it is offset by 1)
    end
end
% Call the builtin file to compute distances for the build-in distance measures
% on non-sparse real float (double or single) data.
if ~(issorted(categCols))
    categCols = sort(categCols);
end

D = internal.stats.pdist2HeterogeneousMEX(X',Y',dist,freqAndECDFX', freqAndECDFY',uint64(categCols));  % do not forget the transpose and the casting


function [lowidx,highidx] = binarySearchIdx(lowidx, highidx, xsorted, key)
% This function is used for continuous features to determine where a value
% in Y might lie in X. Calculates the interval in logN time

% x is a sorted vector, key is a scalar whose position needs to be
% determined
% lowidx and highidx are indices
if (highidx-lowidx)==1 || (highidx-lowidx)==0 % add the edge case when there is only one element in the sorted array
    return;
else
    medianidx = floor((lowidx+highidx)/2);
    if key <= xsorted(medianidx)
        [lowidx, highidx] = binarySearchIdx(lowidx, medianidx,xsorted,key);
    else
        [lowidx, highidx] = binarySearchIdx(medianidx, highidx, xsorted, key);
    end
end
