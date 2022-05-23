function [idx,dist]=rangesearch(obj,Y,radius,varargin)
%RANGESEARCH Radius search using a KDTreeSearcher object.
%   IDX = RANGESEARCH(NS,Y,RADIUS) finds all the points in NS.X that are
%   within distance RADIUS for points in Y. Rows of Y correspond to
%   observations and columns correspond to variables. Y must have the same
%   number of columns as NS.X. RADIUS is a numeric non-negative number
%   specifying the radius threshold. IDX is NY-by-1 cell array, where NY is
%   the number of rows in Y. IDX{I} contains the indices of points in NS.X
%   whose distance to Y(I,:) are not greater than RADIUS, and these indices 
%   are sorted in the ascending order of the corresponding distance values.
%
%   [IDX, D] = RANGESEARCH(NS,Y,RADIUS) returns a NY-by-1 cell array D. D{I}
%   contains the distance values between Y(I,:) and the corresponding
%   points returned in IDX{I}.
%
%   [IDX, D] = RANGESEARCH(NS,Y,RADIUS,'NAME1',VALUE1,...,'NAMEN',VALUEN)
%   specifies optional argument name/value pairs:
%
%     Name          Value
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
%                   distance. Default is true.
%
%   Example:
%      % Create a KDTreeSearcher object for data X using the default
%      % 'euclidean' distance:
%      X = randn(100,5);
%      ns = createns(X,'nsmethod','kdtree');
%      % Find the points in X whose distance are not greater than 1.5 to
%      % the points in Y
%      Y = randn(10, 5);
%      [idx, dist] = rangesearch(ns,Y,1.5)
%
%  See also  RANGESEARCH, KNNSEARCH, KDTreeSearcher, CREATENS,
%  ExhaustiveSearcher.

%   Copyright 2009-2018 The MathWorks, Inc.



if  ~(isscalar(radius) && isnumeric(radius) && radius >= 0 )
    error(message('stats:KDTreeSearcher:rangesearch:BadRadius'));
end

[nX,nDims] = size(obj.X);
[nY,nDims2]= size(Y);

if nDims2 ~= nDims
    error(message('stats:KDTreeSearcher:rangesearch:SizeMisMatch', nDims));
end

[varargin{:}] = convertStringsToChars(varargin{:});
% to do, check whether radius >0 and is numeric
pnames = { 'distance', 'p' , 'sortindices'};
dflts =  {   []         []  true};
[distMetric, minExp, doSort] = ...
     internal.stats.parseArgs(pnames, dflts, varargin{:});

validateattributes(doSort,{'logical','numeric'},{'scalar'},'','sortindices');
doSort = logical(doSort);

if ~isempty(distMetric)
    methods = {'euclidean'; 'cityblock'; 'chebychev'; 'minkowski'};
    distMetric = internal.stats.getParamVal(distMetric,methods,'Distance');
else
    distMetric = obj.Distance; % use the default distance saved in obj.dist
end

if strncmp(distMetric,'min',3) % 'minkowski'
    if  ~isempty(minExp)
        if ~(isscalar(minExp) && isnumeric(minExp) && minExp > 0 )
            error(message('stats:KDTreeSearcher:rangesearch:BadMinExp'))
        end
    elseif strncmp(obj.Distance,'minkowski',3) && ...
            ~isempty(obj.DistParameter)
        minExp = obj.DistParameter;
    else
        minExp = 2;
    end
else% 'euclidean', 'cityblock' or 'chebychev' distance
    if ~isempty(minExp)
        error(message('stats:KDTreeSearcher:rangesearch:InvalidMinExp'));
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
    warning(message('stats:KDTreeSearcher:rangesearch:DataConversion', class( Y ), outClass));
end

Y = cast(Y,outClass);
if ~strcmp(outClass, class(obj.X))
    %only happens when X is double and Y is single
    X2 = cast(obj.X, outClass)';
else
    X2 = obj.X';
end

if ~isreal(Y)
    error(message('stats:KDTreeSearcher:rangesearch:ComplexData'));
end

if issparse(Y)
    warning(message('stats:KDTreeSearcher:rangesearch:DataConversion', 'sparse', 'full'));
    Y = full(Y);
end
% Degenerate case, just cell of the proper size.
if(nX == 0|| nY== 0)
    idx =repmat({zeros(1,0,outClass)},nY,1);
    if nargout > 1
        dist=repmat({zeros(1,0,outClass)},nY,1);
    end
    return;
end

wasNaNY= any(isnan(Y),2);

if nargout < 2
    idx= internal.stats.KDTreeSearcher.knnsearchmex(X2, Y', [], minExp, obj.cutDim, obj.cutVal, ...
        obj.lowerBounds', obj.upperBounds',obj.leftChild, obj.rightChild, ...
        obj.leafNode,obj.idx,obj.nx_nonan,wasNaNY, [],radius,doSort);
else
    [idx, dist]= internal.stats.KDTreeSearcher.knnsearchmex(X2, Y',[],minExp, obj.cutDim, obj.cutVal, ...
        obj.lowerBounds', obj.upperBounds',obj.leftChild, obj.rightChild, ...
        obj.leafNode, obj.idx,obj.nx_nonan, wasNaNY, [],radius,doSort);
    
end
