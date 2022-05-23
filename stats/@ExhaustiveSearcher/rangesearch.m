function [idx,dist]=rangesearch(obj,Y,radius,varargin)
%RANGESEARCH Radius search using an ExhaustiveSearcher object.
%   IDX = RANGESEARCH(NS,Y,RADIUS) finds all the points in NS.X that are
%   within distance RADIUS for points in Y. Rows of Y correspond to
%   observations and columns correspond to variables. Y must have the same
%   number of columns as NS.X. RADIUS is a numeric non-negative number
%   specifying the radius threshold. IDX is NY-by-1 cell array, where NY is
%   the number of rows in Y. IDX{I} contains the indices of points in NS.X
%   whose distance to Y(I,:) are not greater than RADIUS. These indices 
%   are sorted in the ascending order of the corresponding distance values.
%
%   [IDX, D] = RANGESEARCH(NS,Y,RADIUS) returns a NY-by-1 cell array D.
%   D{I} contains the distance values between Y(I,:) and the corresponding
%   points returned in IDX{I}.
%
%   [IDX, D] = RANGESEARCH(NS,Y,RADIUS,'NAME1',VALUE1,...,'NAMEN',VALUEN)
%   specifies optional argument name/value pairs.
%
%     Name          Value
% 
%     'Distance'    A string or a function handle specifying the distance
%                   metric. The value can be one of the following (default
%                   is NS.Distance):
%                     'euclidean'   - Euclidean distance.
%                     'seuclidean'  - Standardized Euclidean distance. Each
%                                     coordinate difference between X and a
%                                     query point is scaled by dividing by
%                                     a scale value S. The default value of
%                                     S is NS.DistParameter if NS.Distance
%                                     is 'seuclidean', otherwise the
%                                     default is the standard deviation
%                                     computed from X, S=NANSTD(X). To
%                                     specify another value for S, use the
%                                     'Scale' argument.
%                     'cityblock'   - City Block distance.
%                     'chebychev'   - Chebychev distance (maximum
%                                     coordinate difference).
%                     'minkowski'   - Minkowski distance. Default is
%                                     NS.DistParameter if NS.Distance is
%                                     'minkowski', or 2 otherwise. To specify
%                                     a different exponent, use the 'P'
%                                     argument.
%                     'mahalanobis' - Mahalanobis distance, computed using
%                                     a positive definite covariance matrix
%                                     C. Default is NS.DistParameter if
%                                     NS.Distance is 'mahalanobis', or
%                                     NANCOV(X) otherwise. To specify
%                                     another value for C, use the 'Cov'
%                                     argument.
%                     'cosine'      - One minus the cosine of the included
%                                     angle between observations (treated
%                                     as vectors).
%                     'correlation' - One minus the sample linear
%                                     correlation between observations
%                                     (treated as sequences of values).
%                     'spearman'    - One minus the sample Spearman's rank
%                                     correlation between observations
%                                     (treated as sequences of values).
%                     'hamming'     - Hamming distance, percentage of
%                                     coordinates that differ.
%                     'jaccard'     - One minus the Jaccard coefficient,
%                                     the percentage of nonzero coordinates
%                                     that differ.
%                     function      - A distance function specified using @
%                                     (for example @DISTFUN). A distance
%                                     function must be of the form:
%
%                                     function D2 = DISTFUN(ZI, ZJ),
%
%                                     taking as arguments a 1-by-N vector
%                                     ZI containing a single row of X or Y,
%                                     and an M2-by-N matrix ZJ containing
%                                     multiple rows of X or Y, and
%                                     returning an M2-by-1 vector of
%                                     distances D2, whose Jth element is
%                                     the distance between the observations
%                                     ZI and ZJ(J,:).
%
%     'P'           A positive scalar P indicating the exponent for
%                   Minkowski distance. This argument is only valid when
%                   RANGESEARCH uses the 'minkowski' distance metric. Default
%                   is NS.DistParameter if NS.Distance is 'minkowski', or 2
%                   otherwise.
%
%     'Cov'         A positive definite matrix indicating the covariance
%                   matrix when computing the Mahalanobis distance. This
%                   argument is only valid when RANGESEARCH uses the
%                   'mahalanobis' distance metric. Default is
%                   NS.DistParameter if NS.Distance is 'mahalanobis', or
%                   NANCOV(X) otherwise.
%
%     'Scale'       A vector S containing non-negative values, with length
%                   equal to the number of columns in X. Each coordinate
%                   difference between X and a query point is scaled by the
%                   corresponding element of S when computing the
%                   standardized Euclidean distance. This argument is only
%                   valid when 'Distance' is 'seuclidean'. Default is
%                   NS.DistParameter if NS.Distance is 'seuclidean', or
%                   NANSTD(X) otherwise.
%
%    'SortIndices'  A flag to indicate if output distances and the
%                   corresponding indices should be sorted in the order of 
%                   distances ranging from the smallest to the largest 
%                   distance.  Default is true.
% 
%   Example:
%      % Create an ExhaustiveSearcher object for data X using the cosine 
%      % distance. The kd-tree method does not support radius search for
%      % the cosine distance, therefore CREATENS creates an 
%      % ExhaustiveSearcher object.
%      X = randn(100,5);
%      Y = randn(25,5);
%      ns = createns(X, 'distance', 'cosine');
%
%      % Find the points in X whose cosine distance are not greater than
%      % 0.3 to the points in Y
%      [idx, dist] = rangesearch(ns,Y,0.3);
%
%      % Choose a different distance metric from the one used when creating
%      % this ExhaustiveSearcher object
%      [idx2, dist2] = rangesearch(ns,Y,0.3, 'distance','correlation');
%
%   See also RANGESEARCH, ExhaustiveSearcher, CREATENS, KNNSEARCH
%   KDTreeSearcher, PDIST2.

%   Copyright 2011-2018 The MathWorks, Inc.


% if nargin < 3
%     error(message('stats:ExhaustiveSearcher:rangesearch:TooFewInputs'));
% end

[~,nDims] = size(obj.X);
[~,nDims2]= size(Y);

if nDims2 ~= nDims
    error(message('stats:ExhaustiveSearcher:rangesearch:SizeMisMatch', nDims));
end

[varargin{:}] = convertStringsToChars(varargin{:});
pnames = {  'distance'  'p', 'cov', 'scale', 'sortindices'};
dflts =  {    []           []   []     []    true};
[distMetric, minExp, mahaCov, seucInvWgt, doSort] = ...
    internal.stats.parseArgs(pnames, dflts, varargin{:});

validateattributes(doSort,{'logical','numeric'},{'scalar'},'','sortindices');
doSort = logical(doSort);

if (~isnumeric(radius) || ~isscalar(radius) || radius < 0)
    error(message('stats:ExhaustiveSearcher:rangesearch:BadRadius'));
end
 
if isempty(distMetric)
      distMetric = obj.Distance;
elseif ischar(distMetric)
     methods = {'euclidean'; 'seuclidean'; 'cityblock'; 'chebychev'; ...
        'mahalanobis'; 'minkowski'; 'cosine'; 'correlation'; ...
        'spearman'; 'hamming'; 'jaccard'};
    i =  find(strncmpi(distMetric, methods, length(distMetric)));
    if length(i) > 1
        error(message('stats:ExhaustiveSearcher:rangesearch:BadDistance', distMetric));
    elseif isempty(i)
        error(message('stats:ExhaustiveSearcher:rangesearch:UnrecognizedDistance', distMetric));
    else
        distMetric = methods{i};
    end

end

arg ={};
if ischar(distMetric)
 
    checkExtraArg(distMetric,minExp,seucInvWgt,mahaCov);
    
    switch distMetric(1:3)
        case 'min'
            if  isempty(minExp) && ~isempty(obj.DistParameter) && ....
                    strncmp(obj.Distance,'minkowski',3)
                minExp = obj.DistParameter;
            end
            arg = {minExp};

        case 'seu'
            if isempty(seucInvWgt) && ~isempty(obj.DistParameter) &&...
                    strncmp(obj.Distance,'seu',3)
                seucInvWgt = obj.DistParameter;
            end
            arg = {seucInvWgt};
            
        case 'mah'
            if isempty(mahaCov) && ~isempty(obj.DistParameter) &&...
                    strncmp(obj.Distance,'mah',3)
                mahaCov = obj.DistParameter;
            end
            arg = {mahaCov};
    end
    
else %no built-in distance
    checkExtraArg('userDist',minExp,seucInvWgt,mahaCov);
end

if nargout < 2
    [~,idx]= pdist2(obj.X,Y, distMetric, arg{:},'radius',radius,'sortindices',doSort);
else
    [dist,idx] = pdist2(obj.X,Y, distMetric, arg{:}, 'radius',radius,'sortindices',doSort);
    dist = dist';
end
idx = idx';

end %method RANGESEARCH

%Give an error if an extra input is provided
function checkExtraArg(distMetric,minExp,seucInvWgt,mahaCov)
if ~isempty(minExp) && ~strncmp(distMetric,'min',3)
    error(message('stats:ExhaustiveSearcher:rangesearch:InvalidMinExp'));
end
if ~isempty(seucInvWgt) && ~strncmp(distMetric,'seu',3)
    error(message('stats:ExhaustiveSearcher:rangesearch:InvalidSeucInvWgt'));
end
if ~isempty(mahaCov) && ~strncmp(distMetric,'mah',3)
    error(message('stats:ExhaustiveSearcher:rangesearch:InvalidMahaCov'));
end

end


