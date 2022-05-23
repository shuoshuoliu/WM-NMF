classdef ExhaustiveSearcher < stats.coder.searcher.NeighborSearcher
    %#codegen
    %ExhaustiveSearcher Codegen compatible Neighbor search object using exhaustive search.
    %   An ExhaustiveSearcher object performs KNN (K-nearest-neighbor) search
    %   or radius search using exhaustive search.
    %
    %   Rows of X correspond to observations and columns correspond to
    %   variables. When one of the above methods creates an ExhaustiveSearcher
    %   object, it saves X. The KNNSEARCH method and RANGESEARCH method compute
    %   the distance values from all the points in X to the query points to
    %   find requested points.
    %
    %   ExhaustiveSearcher properties:
    %       X               - Data used to create the object.
    %       Distance        - The distance metric.
    %       DistParameter   - The additional parameter for the distance metric.
    %
    %   ExhaustiveSearcher methods:
    %       ExhaustiveSearcher - Construct an ExhaustiveSearcher object
    %       knnsearch          - Find nearest neighbors for query points.
    %       rangesearch        - Find points within a given radius for query
    %                            points.

    %   Copyright 2017 The MathWorks, Inc.
    properties
        X
        Distance
        DistParameter
    end
    
    properties(GetAccess=protected, SetAccess=private)
        %checkNegativeDistance = false;
    end
    
    methods(Access = protected)
        function obj = ExhaustiveSearcher(X,varargin)
            %ExhaustiveSearcher Construct an ExhaustiveSearcher object.
            %   NS = ExhaustiveSearcher(X,'NAME1',VALUE1,...,'NAMEN',VALUEN) creates
            %   an ExhaustiveSearcher object, specifying the following optional
            %   argument name/value pairs:
            %
            %     Name         Value
            %     'Distance'   A string or a function handle specifying the default
            %                  distance metric when you call the KNNSEARCH method.
            %                  The value can be one of the following:
            %                  'euclidean'   - Euclidean distance (default).
            %                  'seuclidean'  - Standardized Euclidean distance. Each
            %                                  coordinate difference between X and a
            %                                  query point is scaled by dividing by a
            %                                  scale value S. The default value
            %                                  of S is the standard deviation computed
            %                                  from X, S=NANSTD(X). To specify another
            %                                  value for S, use the 'Scale' argument.
            %                 'cityblock'   -  City Block distance.
            %                 'chebychev'   -  Chebychev distance (maximum coordinate
            %                                  difference).
            %                 'minkowski'   -  Minkowski distance. The default exponent
            %                                  is 2. To specify a different exponent,
            %                                  use the 'P' argument.
            %                 'mahalanobis' -  Mahalanobis distance, computed using a
            %                                  positive definite covariance matrix C.
            %                                  The default value of C is the sample
            %                                  covariance matrix, as computed by
            %                                  NANCOV(X). To specify another value for
            %                                  C, use the 'COV' argument.
            %                 'cosine'      -  One minus the cosine of the included
            %                                  angle between observations (treated as
            %                                  vectors).
            %                 'correlation' -  One minus the sample linear
            %                                  correlation between observations
            %                                  (treated as sequences of values).
            %                 'spearman'    -  One minus the sample Spearman's rank
            %                                  correlation between observations
            %                                  (treated as sequences of values).
            %                 'hamming'     -  Hamming distance, percentage of
            %                                  coordinates that differ.
            %                 'jaccard'     -  One minus the Jaccard coefficient, the
            %                                  percentage of nonzero coordinates that
            %                                  differ.
            %                  distance     -  A distance function specified using @
            %                                  (for example @DISTFUN). A distance
            %                                  function must be of the form:
            %
            %                                  function D2 = DISTFUN(ZI, ZJ),
            %
            %                                  taking as arguments a 1-by-N vector ZI
            %                                  containing a single row from X or from
            %                                  the query points Y, and an M2-by-N
            %                                  matrix ZJ containing multiple rows of X
            %                                  or Y, and returning an M2-by-1 vector of
            %                                  distances D2, whose Jth element is the
            %                                  distance between the observations ZI and
            %                                  ZJ(J,:).
            %
            %    'P'          A positive scalar indicating the exponent for Minkowski
            %                 distance. This argument is only valid when 'Distance' is
            %                 'minkowski'. Default is 2.
            %
            %    'Cov'        A positive definite matrix indicating the covariance
            %                 matrix when computing the Mahalanobis distance. This
            %                 argument is only valid when 'Distance' is 'mahalanobis'.
            %                 Default is NANCOV(X).
            %
            %    'Scale'      A vector S containing non-negative values, with length
            %                 equal to the number of columns in X. Each coordinate
            %                 difference between X and a query point is scaled by the
            %                 corresponding element of S when computing the
            %                 standardized Euclidean distance. This argument is only
            %                 valid when 'Distance' is 'seuclidean'. Default is
            %                 NANSTD(X).
            %
            narginchk(5,5);
            nx = size(X,2);
            
            % Parse PV pairs and validate
            [~,~,~,distance,distanceparam] = stats.coder.distutils.validateNSParams('XSize',nx,varargin{:});
            coder.internal.errorIf(~coder.internal.isConst(ismatrix(X)) || ~ismatrix(X),'stats:pdist2:UnsupportedND');

            coder.internal.errorIf( ~isreal(X),'stats:pdist2:ComplexData');
            classreg.learning.coderutils.checkSupportedNumeric('X',X,false,false,false);            
            obj.Distance = distance;
            obj.X = X;
            obj.DistParameter = distanceparam;
            
        end %  ExhaustiveSearcher constructor
    end
    methods
        function [Idx,D] = knnsearch(obj,Y,varargin)
            %KNNSEARCH Find K nearest neighbors using an ExhaustiveSearcher object.
            %   IDX = KNNSEARCH(NS,Y) finds the nearest neighbor (closest point) in
            %   X=NS.X for each point in Y. Rows of Y correspond to observations and
            %   columns correspond to variables. Y must have the same number of columns
            %   as X. IDX is a column vector with NY rows, where NY is the number of
            %   rows in Y. Each row in IDX contains the index of the observation in X
            %   that has the minimum distance to the corresponding row in Y. The
            %   KNNSEARCH method computes the distance values from all the points in X
            %   to the query points to find nearest neighbors.
            %
            %   [IDX, D] = KNNSEARCH(NS,Y) returns a column vector D containing
            %   the distance between each row of Y and its closest point in X.
            %   That is, D(I) is the distance between X(IDX(I),:) and Y(I,:).
            %
            %   [IDX, D]= KNNSEARCH(NS,Y,'NAME1',VALUE1,...,'NAMEN',VALUEN) specifies
            %   optional argument name/value pairs.
            %
            %     Name          Value
            %     'K'           A positive integer, K, specifying the number of nearest
            %                   neighbors in X for each point in Y. Default is 1. IDX
            %                   and D are NY-by-K matrices. D sorts the distances in
            %                   each row in ascending order. Each row in IDX contains
            %                   the indices of K closest neighbors in X corresponding
            %                   to the K smallest distances in D.
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
            %                   smallest distances in D. It must be a
            %                   compile-time constant for code generation.
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
            %                     Custom distance functions are not
            %                     supported.
            %
            %     'P'           A positive scalar P indicating the exponent for
            %                   Minkowski distance. This argument is only valid when
            %                   KNNSEARCH uses the 'minkowski' distance metric. Default
            %                   is NS.DistParameter if NS.Distance is 'minkowski', or 2
            %                   otherwise.
            %
            %     'Cov'         A positive definite matrix indicating the covariance
            %                   matrix when computing the Mahalanobis distance. This
            %                   argument is only valid when KNNSEARCH uses the
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
            narginchk(2,14);
            coder.internal.errorIf(~coder.internal.isConst(obj),'stats:coder:pdist2:ExpectedConstant','ExhaustiveSearcher Object');
            nx = size(obj.X,2);
            coder.internal.prefer_const(obj);
            [includeTies,K,Distance,P,Cov,Scale] =  stats.coder.distutils.extractNSParams(varargin{:});
            [includeTies,K,~,Distance,AdditionalArg] = stats.coder.distutils.validateNSParams('XSize',nx,...
                                                       'includeTies',includeTies,'Distance',Distance,...
                                                       'K',K,'P',P,'Cov',Cov,'Scale',Scale,...
                                                       'DistExtra',obj.Distance,'DistParamExtra',obj.DistParameter);
            
            coder.internal.prefer_const(Distance,AdditionalArg,K,includeTies);
            [D,Idx] = stats.coder.distutils.pdist2ties(obj.X,Y,Distance,AdditionalArg,K,includeTies);
            
        end
        
        function [Idx,D] = rangesearch(obj,Y,r,varargin)
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
            %                     Custom distance functions are not
            %                     supported.
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
            narginchk(3,11);
            coder.internal.errorIf(~coder.internal.isConst(obj),'stats:coder:pdist2:ExpectedConstant','ExhaustiveSearcher Object');
            nx = size(obj.X,2);
            coder.internal.prefer_const(obj);
            [~,~,Distance,P,Cov,Scale] = stats.coder.distutils.extractNSParams(varargin{:});
            
            [~,~,~,Distance,AdditionalArg] = stats.coder.distutils.validateNSParams('XSize',nx,...
                                         'Distance',Distance,'P',P,'Cov',Cov,'Scale',Scale,...
                                         'DistExtra',obj.Distance,'DistParamExtra',obj.DistParameter);
            
            if ~isempty(AdditionalArg)
                [D, Idx] = pdist2(obj.X,Y,Distance,AdditionalArg,'Radius',r);
            else
                [D, Idx] = pdist2(obj.X,Y,Distance,'Radius',r);
            end
        end
        
    end %for methods
    methods(Static, Hidden)
        
        function outObj = matlabCodegenToRedirected(inObj)
            % static method to return equivalent the target
            % MCOS instance class i.e. this class object,
            % for the given source MCOS instance.
            outObj = stats.coder.searcher.ExhaustiveSearcher(inObj.X,'DistExtra',inObj.Distance,'DistParamExtra',inObj.DistParameter);
        end

        function obj = fromStruct(str)
            % static method to return a
            % stats.coder.searcher.ExhaustiveSearcher object from a given
            % struct
            coder.inline('always');
            coder.internal.prefer_const(str);
            obj =  stats.coder.searcher.ExhaustiveSearcher(str.X,'DistExtra',str.Distance,'DistParamExtra',str.DistParameterValue);
        end
        
    end
end % classdef