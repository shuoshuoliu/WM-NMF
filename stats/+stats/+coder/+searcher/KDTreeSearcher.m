classdef KDTreeSearcher < stats.coder.searcher.NeighborSearcher
    %#codegen
    %KDTreeSearcher Codegen compatible Neighbor search object using kd-tree search.
    %   A KDTreeSearcher object performs KNN (K-nearest-neighbor) search
    %   or radius search using kd-tree search.
    %
    %   Rows of X correspond to observations and columns correspond to
    %   variables. When one of the above methods creates a KDTreeSearcher
    %   object, it it creates and saves a kd-tree based on X. The KNNSEARCH
    %   method and RANGESEARCH method compute the distance values from all
    %   the points in X to the query points to find requested points.
    %
    %   KDTreeSearcher properties:
    %       X               - Data used to create the object.
    %       Distance        - The distance metric.
    %       DistParameter   - The additional parameter for the distance metric.
    %
    %   KDTreeSearcher methods:
    %       KDTreeSearcher  - Construct a KDTreeSearcher object
    %       knnsearch       - Find nearest neighbors for query points.
    %       rangesearch     - Find points within a given radius for query points.

    %   Copyright 2017 The MathWorks, Inc.
    properties
        X               % Input data
        Distance        % Distance method: 'euclidean', 'cityblock', 'chebychev', 'minkowski'
        DistParameter   % Extra parameter for a distance parametr. For minkowski, this stores the exponent.
        nx_nonan        % The number of points with no-NaN values
        cutDim          % The dimension along which each node is split, or 0 for a leaf node
        cutVal          % Cutoff value for split.
        lowerBounds     % The lower bounds of the corresponding node along each dimension
        upperBounds     % The upper bounds of the corresponding node along each dimension
        leftChild       % Left child index for each node
        rightChild      % Right child index for each node
        leafNode        % A logical vector indicating whether the node is a leaf node.
        idxAll          % The point index for all nodes
        idxDim          % Length of each node
        wasnanIdx       % The index that has NaN values
    end
    
    methods(Access = protected)
        function obj = KDTreeSearcher(cgStruct)
            %KDTreeSearcher Construct a KDTreeSearcher object.
            %   NS = KDTreeSearcher(X,'NAME1',VALUE1,...,'NAMEN',VALUEN) creates
            %   a KDTreeSearcher object, specifying the following optional
            %   argument name/value pairs:
            %
            %     Name         Value
            %     'Distance'    A string specifying the default distance metric used
            %                   when you call the KNNSEARCH method. It can be one of
            %                   the following:
            %                     'euclidean'   - Euclidean distance (default)
            %                     'cityblock'   - City Block distance
            %                     'chebychev'   - Chebychev distance (maximum
            %                                     coordinate  difference)
            %                     'minkowski'   - Minkowski distance
            %
            %     'P'           A positive scalar indicating the exponent for Minkowski
            %                   distance. This argument is only valid when 'Distance' is
            %                   'minkowski'. Default is 2.
            %            
            
            % Validate parameters
            validateKSParams(cgStruct);
            
            coder.internal.errorIf(~coder.internal.isConst(ismatrix(cgStruct.X)) || ~ismatrix(cgStruct.X),'stats:pdist2:UnsupportedND');
            
            coder.internal.errorIf( ~isreal(cgStruct.X),'stats:pdist2:ComplexData');
            classreg.learning.coderutils.checkSupportedNumeric('X',cgStruct.X,false,false,false);
            obj.Distance = cgStruct.Distance;
            obj.X = cgStruct.X;
            obj.DistParameter = cgStruct.DistParameterValue;            
            
            obj.nx_nonan = cgStruct.nx_nonan;
            obj.cutDim = cgStruct.cutDim;
            obj.cutVal = cgStruct.cutVal;
            
            obj.lowerBounds = cgStruct.lowerBounds;
            obj.upperBounds = cgStruct.upperBounds;
            if coder.target('MATLAB')
                obj.lowerBounds(cgStruct.lowerBoundsInfIdx) = -inf;
                obj.upperBounds(cgStruct.upperBoundsInfIdx) = inf;
            else
                obj.lowerBounds(cgStruct.lowerBoundsInfIdx) = -coder.internal.inf;
                obj.upperBounds(cgStruct.upperBoundsInfIdx) = coder.internal.inf;
            end
            
            obj.leftChild = cgStruct.leftChild;
            obj.rightChild = cgStruct.rightChild;
            obj.leafNode = cgStruct.leafNode;
            obj.idxAll = cgStruct.idxAll;
            obj.idxDim = cgStruct.idxDim;
            
            obj.wasnanIdx = cgStruct.wasnanIdx;
            
        end %  KDTreeSearcher constructor
    end
    methods
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
            
            narginchk(2,10);
            
            nDims = size(obj.X,2);
            
            coder.internal.errorIf(~coder.internal.isConst(obj),'stats:coder:pdist2:ExpectedConstant','KDTreeSearcher Object');
            
            [includeTies,numNN,distMetric,P,Cov,Scale,bucketsize] =  stats.coder.distutils.extractNSParams(varargin{:});
            
            [includeTies,numNN,~,distMetric,~] = stats.coder.distutils.validateNSParams('XSize',nDims,...
                'includeTies',includeTies,'Distance',distMetric,'K',numNN,'P',P,...
                'DistExtra',obj.Distance,'DistParamExtra',obj.DistParameter);
            
            % Scale and Cov for knnsearch method are not validated in validateNSParams.
            coder.internal.errorIf(~isempty(Cov) || ~isempty(Scale),'stats:coder:knnsearch:covScaleNotNeeded');
            
            % BucketSize is only valid for the knnsearch function and not
            % the method
            coder.internal.errorIf(~isempty(bucketsize),'stats:coder:knnsearch:bSizeNotNeededObj');
            
            [idx, dist] = stats.coder.distutils.kdsearchfun(obj, Y, coder.internal.indexInt(numNN), distMetric, includeTies,P);
        end
        
        
        function [idx,dist] = rangesearch(obj,Y,radius,varargin)
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
            
            narginchk(3,11);
            coder.internal.errorIf(~coder.internal.isConst(obj),'stats:coder:pdist2:ExpectedConstant','KDTreeSearcher Object');
            nDims = size(obj.X,2);
            
            coder.internal.prefer_const(obj);
            [~,~,distMetric,P,Cov,Scale,bucketsize] = stats.coder.distutils.extractNSParams(varargin{:});
            
            [~,~,~,distMetric,~] = stats.coder.distutils.validateNSParams('XSize',nDims,...
                'Distance',distMetric,'P',P,'DistExtra',obj.Distance,'DistParamExtra',obj.DistParameter);
            
            % Scale and Cov for knnsearch method are not validated in validateNSParams.
            coder.internal.errorIf(~isempty(Cov) || ~isempty(Scale),'stats:coder:knnsearch:covScaleNotNeeded');
            
            % BucketSize is only valid for rangesearch function and not the method
            coder.internal.errorIf(~isempty(bucketsize),'stats:coder:knnsearch:bSizeNotNeededObj');
            
            [idx, dist]= stats.coder.distutils.kdsearchfun(obj, Y, [], distMetric, [], P, radius);
        end
        
    end %for methods
    methods (Static)
        function obj = fromStruct(cgStruct)
            % static method to return a
            % stats.coder.searcher.KDTreeSearcher object from a given
            % struct
            coder.inline('always');
            coder.internal.prefer_const(cgStruct);
            obj =  stats.coder.searcher.KDTreeSearcher(cgStruct);
        end
    end
    methods(Static, Hidden)
        
        function outObj = matlabCodegenToRedirected(inObj)
            % static method to return equivalent the target
            % MCOS instance class i.e. this class object,
            % for the given source MCOS instance.
            cgStruct = toStruct(inObj);
            outObj = stats.coder.searcher.KDTreeSearcher(cgStruct);
        end
        
    end
end % classdef

function validateKSParams(cgStruct)

% Validate fields of Struct

coder.inline('always');

validateattributes(cgStruct.Distance,{'char'},{'nonempty','row'},mfilename,'Distance');
 
validateattributes(cgStruct.cutDim,{'numeric'},{'real'},mfilename,'cutDim');

validateattributes(cgStruct.cutVal,{'numeric'},{'real'},mfilename,'cutVal');

validateattributes(cgStruct.lowerBounds,{'numeric'},{'real'},mfilename,'lowerBounds');

validateattributes(cgStruct.upperBounds,{'numeric'},{'real'},mfilename,'upperBounds');

validateattributes(cgStruct.leftChild,{'numeric'},{'real'},mfilename,'leftChild');

validateattributes(cgStruct.rightChild,{'numeric'},{'real'},mfilename,'rightChild');

validateattributes(cgStruct.leafNode,{'logical'},{'nonnan','real','finite'},mfilename,'leafNode');

validateattributes(cgStruct.idxAll,{'numeric'},{'nonnan','integer','real','finite'},mfilename,'idxAll');

validateattributes(cgStruct.idxDim,{'numeric'},{'nonnan','integer','real','finite'},mfilename,'idxDim');

validateattributes(cgStruct.nx_nonan,{'numeric'},{'scalar','integer','real','finite'},mfilename,'nx_nonan');

validateattributes(cgStruct.wasnanIdx,{'numeric'},{'nonnan','real','finite'},mfilename,'wasnanIdx');

end