function [idx,dist]=knnsearch(obj,Y,varargin)
%KNNSEARCH Find K nearest neighbors using a KDTreeSearcher object.
%   IDX = KNNSEARCH(NS,Y) finds the nearest neighbor (closest point) in
%   NS.X for each point in Y. Rows of Y correspond to observations and
%   columns correspond to variables. Y must have the same number of columns
%   as NS.X. IDX is a column vector with NY rows, where NY is the number
%   of rows in Y. Each row in IDX contains the index of the observation in
%   NS.X that has the smallest distance to the corresponding observation
%   in Y.
%
%   [IDX, D] = KNNSEARCH(NS,Y) returns a column vector D containing the
%   distances between each observation in Y and its corresponding closest
%   observation in NS.X. That is, D(I) is the distance between
%   NS.X(IDX(I),:) and Y(I,:).
%
%   [IDX, D] = KNNSEARCH(NS,Y,'NAME1',VALUE1,...,'NAMEN',VALUEN)
%   specifies optional argument name/value pairs:
%
%     Name          Value
%     'K'           A positive integer, K, specifying the number of nearest
%                   neighbors in NS.X to find for each point in Y. Default
%                   is 1. IDX and D are NY-by-K matrices. D sorts the
%                   distances in each row in ascending order. Each row in
%                   IDX contains the indices of K closest neighbors in X
%                   corresponding to the K smallest distances in D.
%
%    'IncludeTies'  A logical value indicating whether KNNSEARCH will
%                   include all the neighbors whose distance values are
%                   equal to the Kth smallest distance. Default is false.
%                   If the value is true, KNNSEARCH includes all these
%                   neighbors. In this case, IDX and D are NY-by-1 cell
%                   arrays. Each row in IDX and D contains a vector with at
%                   least K numeric numbers. D sorts the distances in each
%                   vector in ascending order. Each row in IDX contains the
%                   indices of the closest neighbors corresponding to these
%                   smallest distances in D.
%
%    'Distance'     A string specifying the distance metric. The value can
%                   be one of the following:
%                     'euclidean'   - Euclidean distance
%                     'cityblock'   - City Block distance
%                     'chebychev'   - Chebychev distance (maximum
%                                     coordinate difference)
%                     'minkowski'   - Minkowski distance
%                     Default is NS.Distance
%
%    'P'            A positive scalar indicating the exponent of Minkowski
%                   distance. This argument is only valid when KNNSEARCH
%                   uses the 'minkowski' distance metric. Default is
%                   NS.DistParameter if NS.Distance is 'minkowski', or 2
%                   otherwise.
%
%    'SortIndices'  A flag to indicate if output distances and the
%                   corresponding indices should be sorted in the order of 
%                   distances ranging from the smallest to the largest 
%                   distance.  Default is true.
%
%   Example:
%      % Create a KDTreeSearcher object for data X with the 'euclidean'
%      % distance:
%      X = randn(100,5);
%      ns = createns(X,'nsmethod','kdtree');
%
%      % Find 5 nearest neighbors in X and the corresponding distance
%      % values for each point in Y:
%      Y = randn(25, 5);
%      [idx, dist] = knnsearch(ns,Y,'k',5);
%
%  See also KNNSEARCH, KDTreeSearcher, CREATENS, RANGESEARCH, PDIST2.

%   Copyright 2009-2018 The MathWorks, Inc.



[nX,nDims] = size(obj.X);
[nY,nDims2]= size(Y);

if nDims2 ~= nDims
    error(message('stats:KDTreeSearcher:knnsearch:SizeMisMatch', nDims));
end

[varargin{:}] = convertStringsToChars(varargin{:});
pnames = { 'k'  'distance', 'p','includeties'  'sortindices'};
dflts =  { 1    []          [] , false          true};
[numNN,distMetric, minExp,includeTies,doSort] = ...
     internal.stats.parseArgs(pnames, dflts, varargin{:});

validateattributes(doSort,{'logical','numeric'},{'scalar'},'','sortindices');
doSort = logical(doSort);
 
if ~isscalar(numNN) || ~isnumeric(numNN) ||  ...
        numNN <1 || numNN~=round(numNN)
    error(message('stats:KDTreeSearcher:knnsearch:BadK'));
end

if ~islogical(includeTies)
    error(message('stats:KDTreeSearcher:knnsearch:InvalidIncludeTies'));
end

if ~isempty(distMetric)
    methods = {'euclidean'; 'cityblock'; 'chebychev'; 'minkowski'};
    distMetric = internal.stats.getParamVal(distMetric,methods,'Distance');
else
    distMetric = obj.Distance; % use the default distance saved in obj.dist
end

if strncmp(distMetric,'min',3) % 'minkowski'
    if  ~isempty(minExp)
        if ~(isscalar(minExp) && isnumeric(minExp) && minExp > 0 )
            error(message('stats:KDTreeSearcher:knnsearch:BadMinExp'))
        end
    elseif strncmp(obj.Distance,'minkowski',3) && ...
            ~isempty(obj.DistParameter)
        minExp = obj.DistParameter;
    else
        minExp = 2;
    end
else% 'euclidean', 'cityblock' or 'chebychev' distance
    if ~isempty(minExp)
        error(message('stats:KDTreeSearcher:knnsearch:InvalidMinExp'));
    end
    switch distMetric(1:3)
        case 'euc'
            minExp = 2;
        case 'cit'
            minExp = 1;
        case 'che'
            minExp = inf;
    end
end


% Integer/logical/char/anything data will be converted to float. Complex
% floating point data can't be handled.

try
    outClass = superiorfloat(obj.X,Y);
catch
    outClass = class(obj.X);
    warning(message('stats:KDTreeSearcher:knnsearch:DataConversion', class( Y ), outClass));
end

Y = cast(Y,outClass);
if ~strcmp(outClass, class(obj.X))
    %only happens when X is double and Y is single
    X2 = cast(obj.X, outClass)';
else
    X2 = obj.X';
end

if ~isreal(Y)
    error(message('stats:KDTreeSearcher:knnsearch:ComplexData'));
end

if issparse(Y)
    warning(message('stats:KDTreeSearcher:knnsearch:DataConversion', 'sparse', 'full'));
    Y = full(Y);
end

numNN = min(numNN,nX);
% Degenerate case, just return an empty matrix or cell of the proper size.
if (nY == 0 || numNN ==0)
    if ~includeTies
        idx = zeros(nY, numNN);
        dist = zeros(nY, numNN, outClass);
    else
        idx =repmat({zeros(1,0)},nY,1);
        dist=repmat({zeros(1,0,outClass)},nY,1);
    end
    return;
end

numNN2 = min(numNN, nX);
wasNaNY= any(isnan(Y),2);

if numNN2 > 0
    if nargout < 2
        idx = internal.stats.KDTreeSearcher.knnsearchmex(X2, Y', numNN2, minExp, obj.cutDim, obj.cutVal, ...
            obj.lowerBounds', obj.upperBounds',obj.leftChild, obj.rightChild, ...
            obj.leafNode,obj.idx,obj.nx_nonan,wasNaNY, includeTies,[],doSort);
    else
        [idx, dist]= internal.stats.KDTreeSearcher.knnsearchmex(X2, Y',numNN2,minExp, obj.cutDim, obj.cutVal, ...
            obj.lowerBounds', obj.upperBounds',obj.leftChild, obj.rightChild, ...
            obj.leafNode, obj.idx,obj.nx_nonan,wasNaNY, includeTies,[],doSort);
      
    end
    if  ~includeTies
        idx =idx';
        if nargout > 1
            dist = dist';
        end
    end
else %numNN2 ==0
      if  ~includeTies
            idx = zeros(nY,0,outClass);
            if nargout > 1
            dist = zeros(nY,0,outClass);
            end
      else
          idx=cell(nY,1);
          if nargout > 1
              dist = cell(nY,1);
          end
      end
 end

if includeTies && any(wasNaNY)
    %Process the point in Y having NaN values
    NaNYIdx = find(wasNaNY);
    lenNaNY = numel(NaNYIdx);
    idxTemp = repmat((1:numNN), lenNaNY,1);
    idx(NaNYIdx) = mat2cell(idxTemp,ones(lenNaNY,1),numNN);
    if nargout > 1
       distTemp = nan(lenNaNY,numNN, outClass);
       dist(NaNYIdx) = mat2cell(distTemp,ones(lenNaNY,1),numNN);
    end
end
end


