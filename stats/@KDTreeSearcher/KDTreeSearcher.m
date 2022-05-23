classdef KDTreeSearcher < NeighborSearcher
%KDTreeSearcher Neighbor search object using a kd-tree.
%   A KDTreeSearcher object performs KNN (K-nearest-neighbor) search or
%   radius search using a kd-tree. You can create a KDTreeSearcher object
%   based on X using either of the following syntaxes:
%
%   CREATENS function:
%      NS = CREATENS(X,'NSMethod','kdtree')
%   KDTreeSearcher constructor:
%      NS = KDTreeSearcher(X)
%
%   Rows of X correspond to observations and columns correspond to
%   variables. When one of the above methods creates a KDTreeSearcher
%   object, it creates and saves a kd-tree based on X. The KNNSEARCH method
%   and RANGESEARCH use the kd-tree to find requested points in X to the
%   query points.
%
%   KDTreeSearcher properties:
%       X               - Data used to create the KDTreeSearcher object
%       Distance        - The distance metric
%       DistParameter   - The additional parameter for the distance
%                         metric
%       BucketSize      - Bucket size of each leaf node in the kd-tree
%   KDTreeSearcher methods:
%       KDTreeSearcher  - Construct a KDTreeSearcher object
%       knnsearch       - Find nearest neighbors for query points
%       rangesearch     - Find points within a given radius for query
%                         points.
%
%    
%   Example:
%      % Create a KDTreeSearcher object for data X with the 'euclidean' 
%      % distance.
%      X = randn(100,5);
%      ns = createns(X,'nsmethod','kdtree');
%
%      % Find 5 nearest neighbors in X and the corresponding distance
%      % values for each point in Y.
%      Y = randn(25, 5);
%      [idx, dist] = knnsearch(ns,Y,'k',5);
%
%      % Find the points in X whose distance are not greater than 0.3 to
%      % the points in Y.
%      idx = rangesearch(ns,Y,0.3);
%
%   See also CREATENS, ExhaustiveSearcher, STATS/KNNSEARCH,
%   STATS/RANGESEARCH.

%   References:
%     Friedman, J. H., Bentely, J. and Finkel, R. A. (1977) An Algorithm
%     for Finding Best Matches in Logarithmic Expected Time, ACM
%     Transactions on Mathematical Software 3, 209. 

%   Copyright 2009-2017 The MathWorks, Inc.

    
    properties (GetAccess='private', SetAccess='private')
        
        nx_nonan = 0; %the number of points with no-NaN values
        wasnanIdx = []; %the index that has NaN values
        
        % The dimension along which each node is split, or 0 for a leaf
        % node.
        cutDim = zeros(0,1);
        
        %cutoff value for split. For non-leaf node, the points will go to left
        %child node if its values on the split dimension is not greater
        %than this cutoff value; otherwise, they will go the right child
        %node
        cutVal = zeros(0,1);
        
        %The lower bounds of the corresponding node along each dimension.
        lowerBounds = zeros(0,1);
        
        %The upper bounds of the corresponding node along each dimension.
        upperBounds = zeros(0,1);
        
        %Left child index for each node
        leftChild = zeros(0,1);
        
        %Right child index for each node
        rightChild = zeros(0,1);
        
        % A logical vector indicating whether the node is a leaf node. TRUE
        % for a leaf node.
        leafNode = false(0,1);
        
        %the point index for each node
        idx = {};
        
        %The total number of nodes in the kd-tree.
        numNodes = 1;
    end
    
    
    properties (GetAccess='public', SetAccess='private')
        %BucketSize Bucket size of each leaf node in the kd-tree.
        %   The BucketSize property is a positive integer, indicating the
        %   maximum number of data points in each leaf node of the kd-tree.
        %
        %   See also KDTreeSearcher. 
        BucketSize = 50;
    end
    
    methods (Access = protected)
        function this = setDistance(this, distMetric)
            distMetric = convertStringsToChars(distMetric);
             if ischar(distMetric)
              methods = {'euclidean';  'cityblock'; 'chebychev'; ...
                    'minkowski'; };
                i =  find(strncmpi(distMetric, methods, length(distMetric)));
                if length(i) > 1
                     error(message('stats:KDTreeSearcher:AmbigiousDistance',distMetric)); 
                elseif isempty(i)
                     error(message('stats:KDTreeSearcher:BadDistance'));
                                   
                end
                this.PrivDistance =  methods{i};   
                if i < 4
                    this.PrivDistParameter = [];
                else
                    this.PrivDistParameter = 2;
                end
             else
                 error(message('stats:KDTreeSearcher:BadDistance'));
             end
                
        end
        
          function this = setDistParameter(this, minExp)
            i = find(strncmpi(this.Distance, {'minkowski'},3),1);
            if isempty(i)
                 error(message('stats:KDTreeSearcher:InvalidDistanceParam'));
            end
           
            if ~(isscalar(minExp) && isnumeric(minExp) && minExp > 0 )
                   error(message('stats:KDTreeSearcher:BadMinExp'))
            end
            
            this.PrivDistParameter = minExp;
      
        end
    end
    
    methods %(Access='public')
        
        function obj = KDTreeSearcher(X, varargin)
%KDTreeSearcher Construct a KDTreeSearcher object.
%   NS = KDTreeSearcher(X,'NAME1',VALUE1,...,'NAMEN',VALUEN) creates
%   a KDTreeSearcher object, specifying the following optional argument
%   name/value pairs:
%
%     Name          Value
%     'Distance'    A string specifying the default distance metric used
%                   when you call the KNNSEARCH method. It can be one of
%                   the following:
%                     'euclidean'   - Euclidean distance (default)
%                     'cityblock'   - City Block distance
%                     'chebychev'   - Chebychev distance (maximum
%                                     coordinate  difference)
%                     'minkowski'   - Minkowski distance
%
%    'P'            A positive scalar, indicating the exponent of
%                   Minkowski distance. This argument is only
%                   valid when 'Distance' is 'minkowski'. Default
%                   is 2.
%
%    'BucketSize'   A positive integer, indicating the maximum number of
%                   data points in each leaf node of the kd-tree.
%                   Default is 50.
%
%   See also  KDTreeSearcher, CREATENS.

            if nargin == 0
                error(message('stats:KDTreeSearcher:NoDataInput'));
            end
            
            pnames = {'distance' 'bucketsize' 'p'};
            dflts =  { 'euclidean'        50   []};
            [varargin{:}] = convertStringsToChars(varargin{:});
            [dist,bSize,minExp] = ...
                 internal.stats.parseArgs(pnames, dflts, varargin{:});
                        
            methods = {'euclidean'; 'cityblock'; 'chebychev'; 'minkowski'};
            dist = internal.stats.getParamVal(dist,methods,'Distance');
            obj.Distance = dist;
            
            if isempty(bSize)
                %accept empty bucketSize value and assign the default number
                %to it, because createns may pass empty value when calling
                %KDTreeSearcher constructor
                bSize = 50;
            elseif ~(isscalar(bSize) && isnumeric(bSize) && bSize >= 1 && ...
                    round(bSize) == bSize)
                error(message('stats:KDTreeSearcher:BadBucketSize'))
            end
            obj.BucketSize = bSize;
            
            if strcmp(dist,'minkowski')
                if  ~isempty(minExp)
                    if ~(isscalar(minExp) && isnumeric(minExp) && minExp > 0 )
                        error(message('stats:KDTreeSearcher:BadMinExp'))
                    end
                    obj.DistParameter = minExp;
                else
                    obj.DistParameter = 2;
                end
            elseif ~isempty(minExp)
                error(message('stats:KDTreeSearcher:InvalidMinExp'));
            end
            
            % Non-float data, such as Integer/logical/char/anything data
            % will be converted double.
            if ~isfloat(X)
                warning(message('stats:KDTreeSearcher:DataConversion', class( X ), 'double'));
                X = double(X);
            end
            %  Complex floating point data can't be handled
            if ~isreal(X)
                error(message('stats:KDTreeSearcher:ComplexData'));
            end
            
            %convert sparse data to full, KDTreeSearcher can't handle
            %sparse data directly
            if issparse(X)
                warning(message('stats:KDTreeSearcher:DataConversion', 'sparse', 'full'));
                X = full(X);
            end
            
            obj.X = X;
            
            [nx, nDims] = size(X);
            
            wasnan = any(isnan(X),2)';
            
            hasNaNs = any(wasnan);
            notnan = 1:nx;
            if hasNaNs
                notnan = find(~wasnan); %index of points with no missing values.
                obj.wasnanIdx = find(wasnan);
                nx = numel(notnan); %the number points with no missing values
            end
            
            obj.nx_nonan =nx;
            
            %M is the number of maximal nodes if we choose to cut at the
            %median at each cutting dimension. If the tree is not split at
            %the median in each non-leaf node, the tree may contain more
            %than M nodes.
            M = 2^(ceil(log2(max(nx/obj.BucketSize,1)))+1) - 1;
            
            % Keeping accessing the properties of KDTreeSearcher inside a loop
            % seems affect the performance, therefore We created some temporary
            % variables.
            % The dimension along which each node is split, or 0 for a leaf node.
            cutDimTemp = zeros(M,1);
            %cutoff value for split. The points go to left child node if its
            %values on the split dimension is not greater than this cutoff
            %value.
            cutValTemp = zeros(M,1);
            % each row specifies the lower bounds of the corresponding node along
            % each dimension.
            lowerBoundsTemp = -Inf(M,nDims);
            % each row specifies the upper bounds of the corresponding node along
            % each dimension.
            upperBoundsTemp = Inf(M,nDims);
            % A column vector indicating the left child index.
            leftChildTemp = zeros(M,1);
            % A column vector indicating the right child index.
            rightChildTemp= zeros(M,1);
            % A logical vector indicating whether the node is a leaf node. TRUE
            % for a leaf node.
            leafNodeTemp = false(M,1);
            %Each row is a double vector indicating the data points belong to the
            %corresponding node. Empty for a non-leaf node.
            idxTemp = cell(M,1);
            %only data with no NaNs will be used to create the kd-tree
            idxTemp{1} = notnan;
            currentNode = 1;
            nextUnusedNode = 2; %the next un-used node number
            %start to build the kd-tree
            while(currentNode < nextUnusedNode)
                % if we don't cut at the median, then the following check is
                % required
                %  if currentNode > M
                %  %   we need to expand the size of all the temporary variables
                %
                %  end
                currentIdx = idxTemp{currentNode};
                nPoints = numel(currentIdx);
                if nPoints <= bSize % obj.BucketSize
                    %become a leaf node
                    leafNodeTemp(currentNode) = true;
                else %non-leaf node
                    %find the cutting dimension with the largest spread
                    
                    [~,cuttingDim] = max(max(X(currentIdx,:),[],1) - ...
                                         min(X(currentIdx,:),[],1),[],2);
                    
                    %choose the median value as the partition value
                 
                    [sx,sidx] = sort(X(currentIdx,cuttingDim));
                    sidx= currentIdx(sidx);
                    half = ceil(size(sx,1)/2);
                    %the cutting threshold. It will not be the median when
                    %sx has odd number of element.
                    p = (sx(half)+sx(half+1))/2;
                    
                    cutDimTemp(currentNode) = cuttingDim;
                    cutValTemp(currentNode) = p;
                    
                    lChild = nextUnusedNode;
                    rChild = nextUnusedNode + 1;
                    leftChildTemp(currentNode) = lChild;
                    rightChildTemp(currentNode) = rChild;
                    
                    %Decide the upper bounds of the two children
                    temp = upperBoundsTemp(currentNode,:);
                    %right child keeps the parent's upper bounds
                    upperBoundsTemp(rChild,:) = temp ;
                    temp(cuttingDim) = p;
                    upperBoundsTemp(lChild,:) = temp;
                    %Decide the lower bounds of the two children
                    temp = lowerBoundsTemp(currentNode,:);
                    %left child keeps the parent's lower bounds
                    lowerBoundsTemp(lChild,:) = temp;
                    temp(cuttingDim)= p;
                    lowerBoundsTemp(rChild,:) = temp;
                    %add the data points of the current node to its
                    %children
                    idxTemp{currentNode} = [];
                    %   left = X(currentIdx,cuttingDim) <= p;
                    %   idxTemp{lChild} = currentIdx(left);
                    %   idxTemp{rChild} = currentIdx(~left);
                    idxTemp{lChild} = sidx(1:half);
                    idxTemp{rChild} = sidx(half+1:end);
                    nextUnusedNode = nextUnusedNode+2;
                end
                currentNode = currentNode + 1;
            end %while
            obj.numNodes = nextUnusedNode - 1 ;
            
            %Assign the temporary variables to the properties
            obj.cutDim = cutDimTemp(1:nextUnusedNode-1);
            obj.cutVal = cutValTemp(1:nextUnusedNode-1);
            obj.lowerBounds = lowerBoundsTemp(1:nextUnusedNode-1,:);
            obj.upperBounds = upperBoundsTemp(1:nextUnusedNode-1,:);
            obj.idx = idxTemp(1:nextUnusedNode-1,:);
            obj.leftChild = leftChildTemp(1:nextUnusedNode-1);
            obj.rightChild = rightChildTemp(1:nextUnusedNode-1);
            obj.leafNode = leafNodeTemp(1:nextUnusedNode-1);
            
        end % constructor
        
    end % methods block
    
    methods(Hidden)
        
        function s = toStruct(this)
            warnState  = warning('query','all');
            warning('off','MATLAB:structOnObject');
            cleanupObj = onCleanup(@() warning(warnState));
            s = struct;
            s.Distance = this.Distance;
            s.BucketSize = this.BucketSize;
            s.X = this.X;
            s.DistParameterValue = this.DistParameter;
            if strcmpi(s.Distance,'minkowski')
                s.DistParameterName = 'P';
            end
            
            s.nx_nonan = this.nx_nonan;
            s.cutDim = this.cutDim;
            s.cutVal = this.cutVal;
            
            s.lowerBounds = this.lowerBounds;
            lowerInfIdx = find(isinf(this.lowerBounds));            
            s.lowerBoundsInfIdx = lowerInfIdx;
            s.lowerBounds(lowerInfIdx) = 0;
            
            s.upperBounds = this.upperBounds;
            upperInfIdx = find(isinf(this.upperBounds));
            s.upperBoundsInfIdx = upperInfIdx;
            s.upperBounds(upperInfIdx) = 0;
            
            s.leftChild = this.leftChild;
            s.rightChild = this.rightChild;            
            s.leafNode = this.leafNode;
%             s.idx = this.idx;
            % Convert this.idx from cell array to struct
            tempIdx = uint32(zeros(0,1));
            tempDim = uint32(zeros(numel(this.idx),1));
            for c = 1:numel(this.idx)
                tempIdx = [tempIdx;this.idx{c}'];
                tempDim(c) = numel(this.idx{c});
            end
            s.idxAll = tempIdx;
            s.idxDim = tempDim;
%             s.idxStruct = cell2struct(this.idx,'idxVector',2);            
            s.wasnanIdx = this.wasnanIdx;
            
            s.FromStructFcn = 'KDTreeSearcher.fromStruct';
        end
    end
    methods(Hidden, Static)
        function name = matlabCodegenRedirect(~)
            name = 'stats.coder.searcher.KDTreeSearcher';
        end
        
        function obj = fromStruct(s)
            if ~isempty(s.DistParameterValue)
                obj = KDTreeSearcher(s.X,'Distance',s.Distance,s.DistParameterName,...
                    s.DistParameterValue,'BucketSize',s.BucketSize);
            else
                obj = KDTreeSearcher(s.X,'Distance',s.Distance,'BucketSize',s.BucketSize);
            end
            
        end
    end
end % classdef
