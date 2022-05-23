function [idx, dist] = rangesearch(X,Y,radius,varargin)
%RANGESEARCH Radius search.
%   IDX = RANGESEARCH(X,Y,RADIUS) finds all the points in X that are
%   within distance RADIUS for points in Y. Rows of X and Y correspond to
%   observations, and columns correspond to variables. Y must have the same
%   number of columns as X. RADIUS is a numeric non-negative number
%   specifying the radius threshold. IDX is NY-by-1 cell array, where NY is
%   the number of rows in Y. IDX{I} contains the indices of points in X
%   whose distance to Y(I,:) are not greater than RADIUS, and these indices 
%   are sorted in the ascending order of the corresponding distance values.
%
%   [IDX, D] = RANGESEARCH(X,Y,RADIUS) returns a NY-by-1 cell array D. D{I}
%   contains the distance values between Y(I,:) and the corresponding
%   points returned in IDX{I}.
%
%   [IDX, D]= RANGESEARCH(X,Y,'NAME1',VALUE1,...,'NAMEN',VALUEN) specifies
%   optional argument name/value pairs:
%
%     Name          Value
%
%     'NSMethod'    Nearest neighbors search method. Value is either:
%                   'kdtree'    - Creates and uses a kd-tree to find
%                                 nearest neighbors. 'kdtree' is only valid
%                                 when the distance metric is one of the
%                                 following metrics:
%                                   - 'euclidean'
%                                   - 'cityblock'
%                                   - 'minkowski'
%                                   - 'chebychev'
%                   'exhaustive' - Uses the exhaustive search algorithm.
%                                  The distance values from all the points
%                                  in X to each point in Y are computed to
%                                  find nearest neighbors.
%                   Default is 'kdtree' when the number of columns of X is
%                   not greater than 10, X is not sparse, and the distance
%                   metric is one of the above 4 metrics; otherwise,
%                   default is 'exhaustive'.
%
%    'Distance'     A string or a function handle specifying the distance
%                   metric. The value can be one of the following:
%                   'euclidean'   - Euclidean distance (default).
%                   'seuclidean'  - Standardized Euclidean distance. Each
%                                   coordinate difference between X and a
%                                   query point is scaled by dividing by a
%                                   scale value S. The default value of S
%                                   is the standard deviation computed from
%                                   X, S=NANSTD(X). To specify another
%                                   value for S, use the 'Scale' argument.
%                   'cityblock'   - City Block distance.
%                   'chebychev'   - Chebychev distance (maximum coordinate
%                                   difference).
%                   'minkowski'   - Minkowski distance. The default
%                                   exponent is 2. To specify a different
%                                   exponent, use the 'P' argument.
%                   'mahalanobis' - Mahalanobis distance, computed using a
%                                   positive definite covariance matrix C.
%                                   The default value of C is the sample
%                                   covariance matrix of X, as computed by
%                                   NANCOV(X). To specify another value for
%                                   C, use the 'Cov' argument.
%                   'cosine'      - One minus the cosine of the included
%                                   angle between observations (treated as
%                                   vectors).
%                   'correlation' - One minus the sample linear
%                                   correlation between observations
%                                   (treated as sequences of values).
%                   'spearman'    - One minus the sample Spearman's rank
%                                   correlation between observations
%                                  (treated as sequences of values).
%                   'hamming'     - Hamming distance, percentage of
%                                   coordinates that differ.
%                   'jaccard'     - One minus the Jaccard coefficient, the
%                                   percentage of nonzero coordinates that
%                                   differ.
%                   function      - A distance function specified using @
%                                   (for example @DISTFUN). A distance
%                                   function must be of the form
%  
%                                   function D2 = DISTFUN(ZI, ZJ),
%  
%                                   taking as arguments a 1-by-N vector ZI
%                                   containing a single row of X or Y, an
%                                   M2-by-N matrix ZJ containing multiple
%                                   rows of X or Y, and returning an
%                                   M2-by-1 vector of distances D2, whose
%                                   Jth element is the distance between the
%                                   observations ZI and ZJ(J,:).
%
%    'P'            A positive scalar indicating the exponent of Minkowski
%                   distance. This argument is only valid when 'Distance'
%                   is 'minkowski'. Default is 2.
%  
%    'Cov'          A positive definite matrix indicating the covariance
%                   matrix when computing the Mahalanobis distance. This
%                   argument is only valid when 'Distance' is
%                  'mahalanobis'. Default is NANCOV(X).
%  
%    'Scale'        A vector S containing non-negative values, with length
%                   equal to the number of columns in X. Each coordinate
%                   difference between X and a query point is scaled by the
%                   corresponding element of S. This argument is only valid
%                   when 'Distance' is 'seuclidean'. Default is NANSTD(X).
%  
%    'BucketSize'   The maximum number of data points in the leaf node of the
%                   kd-tree (default is 50). This argument is only
%                   meaningful when kd-tree is used for finding nearest
%                   neighbors.
%
%    'SortIndices'  A flag to indicate if output distances and the
%                   corresponding indices should be sorted in the order of 
%                   distances ranging from the smallest to the largest 
%                   distance.  Default is true.
% 
%   Example:
%      % Find the points in X whose distance are not greater than 1.5 to
%      % the points in Y, using the default distance metric 'euclidean'.
%      X = randn(100,5);
%      Y = randn(10, 5);
%      [idx, dist] = rangesearch(X,Y,1.5);
%
%   See also CREATENS, ExhaustiveSearcher, KDTreeSearcher, KNNSEARCH,
%   PDIST2.

%   Copyright 2011-2018 The MathWorks, Inc.



if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

pnames = { 'nsmethod' 'bucketsize' 'sortindices'};
dflts =  { []         []        true};
[nsmethod, bSize, doSort, ~,args] =...
     internal.stats.parseArgs(pnames, dflts, varargin{:});

O=createns(X,args{:},'nsmethod', nsmethod,'bucketSize',bSize);
if nargout < 2
    idx = rangesearch(O,Y,radius,'sortindices',doSort);
else
    [idx, dist] = rangesearch(O,Y,radius,'sortindices',doSort);
end
end
