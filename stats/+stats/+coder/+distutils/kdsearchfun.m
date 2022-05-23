function [idx, dist] = kdsearchfun(obj,Y,numNNin,distMetric,includeTies,P,radius)
%#codegen
% kdsearchfun
% Helper function for knnsearch and rangesearch using kd-tree searcher

% Copyright 2017-2020 The MathWorks, Inc.

    coder.inline('always');

    if nargin == 6
        radius = cast([], 'like', obj.X);
    end

    nDims = size(obj.X,2);
    nDims2 = size(Y,2);

    coder.internal.errorIf(nDims2 ~= nDims,'stats:KDTreeSearcher:knnsearch:SizeMisMatch',int32(nDims));

    if coder.target('MATLAB')
        NONFINITES = true;
    else
        NONFINITES = eml_option('NonFinitesSupport');
    end

    if isempty(distMetric)
        distMetric = obj.Distance; % use the default distance saved in obj.dist
    end

    % Extra flag for indicating if nonfinite numbers are not supported
    minExpIsInf = false;

    if strncmpi(distMetric,'min',3) % 'minkowski'
        if  ~isempty(P)
            isNotOk = ~(isscalar(P) && isnumeric(P) && P(1) > 0 );
            coder.internal.errorIf(isNotOk,'stats:KDTreeSearcher:knnsearch:BadMinExp')
            minExp = cast(P(1), 'like', obj.X);
        elseif strncmpi(obj.Distance,'minkowski',3) && ~isempty(obj.DistParameter)
            minExp = cast(obj.DistParameter, 'like', obj.X);
        else
            minExp = cast(2, 'like', obj.X);
        end
    else% 'euclidean', 'cityblock' or 'chebychev' distance
        coder.internal.errorIf(~isempty(P),'stats:KDTreeSearcher:knnsearch:InvalidMinExp');
        switch distMetric(1:3)
          case 'euc'
            minExp = cast(2, 'like', obj.X);
          case 'cit'
            minExp = cast(1, 'like', obj.X);
          case 'che'
            if coder.target('MATLAB')
                minExp = inf;
            else
                if NONFINITES
                    minExp = coder.internal.inf('like', obj.X);
                else
                    minExp = coder.internal.inf('like', obj.X);
                    minExpIsInf = true;
                end
            end
          otherwise
            coder.internal.errorIf(true,'stats:internal:getParamVal:BadValueListChoices',obj.Distance,'Distance','''euclidean'', ''cityblock'', ''chebychev'', ''minkowski''');
        end
    end

    coder.internal.prefer_const(distMetric,numNNin,minExp,includeTies,minExpIsInf);

    if coder.target('MEX') && eml_option('MXArrayCodegen')
        [idx, dist] = kdsearchfunMEX(obj,Y,numNNin,includeTies,minExp,radius);

    else
        [idx, dist] = kdsearchfunCG(obj,Y,numNNin,includeTies,minExp,radius,minExpIsInf);
    end

end


function [idx, dist] = kdsearchfunCG(obj,Y,numNNin,includeTies,minExp,radius,minExpIsInf)

    if coder.target('MATLAB')
        NONFINITES = true;
        inMATLAB = true;
    else
        NONFINITES = eml_option('NonFinitesSupport');
        inMATLAB = false;
    end

    nX = size(obj.X,1);
    nY = size(Y,1);

    % Conversion to float.
    % Complex floating point data can't be handled.

    outClass = superiorfloat(obj.X,Y);

    Y = cast(Y,outClass);
    if ~strcmpi(outClass, class(obj.X))
        %only happens when X is double and Y is single
        X = cast(obj.X, outClass);
    else
        X = obj.X;
    end

    coder.internal.errorIf(~isreal(Y),'stats:KDTreeSearcher:knnsearch:ComplexData');

    if NONFINITES
        wasNaNY = any(isnan(Y),2);
    else
        wasNaNY = false(nY,1);
    end

    if ~isempty(numNNin)
        isNotOk = ~isscalar(numNNin) || ~isnumeric(numNNin) ||  ...
                  numNNin <1 || numNNin~=round(numNNin);

        coder.internal.errorIf(isNotOk,'stats:KDTreeSearcher:knnsearch:BadK');

        numNN1 = coder.internal.indexInt(min(numNNin,nX));
        assert(numNN1 <= nX);
        % Degenerate case, just return an empty matrix or cell of the proper size.
        if (nY == 0 || numNN1 ==0)
            if ~includeTies
                idx = zeros(nY, numNN1, indexType);
                dist = zeros(nY, numNN1, 'like', Y);
            else
                idx = repmat({zeros(1,0, indexType)},nY,1);
                dist = repmat({zeros(1,0,'like', Y)},nY,1);
            end
            return;
        end

        numNN = coder.internal.indexInt(min(cast(numNN1,'like', obj.nx_nonan), obj.nx_nonan));

        if numNN > 0        %ask for nearest neighbors
            if ~includeTies   % Doesn't need to check ties
                              % Initialize the index matrix
                idx = zeros(nY,numNN1, indexType);

                % Initialize the distance matrix
                dist = zeros(nY,numNN1,'like',Y);

                noNanCol = 1:numNN;
                parfor j = 1:coder.internal.indexInt(nY)
                    if NONFINITES && wasNaNY(j) %The jth point in Y has NaN
                        if inMATLAB
                            dist(j,noNanCol) = NaN;
                        else
                            dist(j,noNanCol) = coder.internal.nan;
                        end
                        idx(j,noNanCol) = cast(noNanCol, indexType);
                    else
                        pq = search_kdtree(obj,X,Y(j,:),numNN,minExp,numNNin,minExpIsInf);
                        distTemp = coorRoot(pq.D,minExp,minExpIsInf);
                        dist(j,noNanCol) = distTemp(1:numel(noNanCol));
                        idx(j,noNanCol) = cast(pq.I(1:numel(noNanCol)), indexType);
                    end

                end
            else    % Need to check ties
                    % Initialize the index matrix
                idx = coder.nullcopy(cell(coder.internal.indexInt(nY),1));
                % Initialize the distance matrix
                dist = coder.nullcopy(cell(coder.internal.indexInt(nY),1));

                extra = coder.internal.indexInt(5);
                if numNN < 3
                    extra = coder.internal.indexInt(2);
                elseif numNN < 5
                    extra = numNN;
                end
                thisK = numNN + extra;
                if thisK > obj.nx_nonan
                    thisK = coder.internal.indexInt(obj.nx_nonan);
                end

                parfor j = 1: coder.internal.indexInt(nY)
                    if ~wasNaNY(j)
                        [dist1,idx1] = search_kdtree_withties(obj, X, Y(j,:), double(numNN), minExp, thisK,...
                                                              coder.internal.indexInt(numNNin)+extra, minExpIsInf);
                        dist{j} = cast([dist1',zeros(1,numNN1 - numNN)], 'like', Y(j,:));
                        idx{j} = cast([(idx1'),zeros(1,numNN1 - numNN)], indexType);
                    end
                end
            end
        else %numNN2 ==0
            if  ~includeTies
                idx = zeros(nY,0, indexType);
                if nargout > 1
                    dist = zeros(nY,0,'like', Y);
                end
            else
                idx = coder.nullcopy(cell(coder.internal.indexInt(nY),1));
                if nargout > 1
                    dist = coder.nullcopy(cell(coder.internal.indexInt(nY),1));
                end
            end
        end
        if includeTies && any(wasNaNY)
            %Process the point in Y having NaN values
            NaNYIdx = find(wasNaNY);
            lenNaNY = numel(NaNYIdx);

            assert(numNN <= nX);
            for c = 1:lenNaNY
                idx{NaNYIdx(c)} = cast(1:numNN, indexType);
            end

            if nargout > 1  && NONFINITES
                if inMATLAB
                    for c = 1:lenNaNY
                        dist{NaNYIdx(c)} = NaN(1,numNN, outClass);
                    end
                else
                    for c = 1:lenNaNY
                        dist{NaNYIdx(c)} = coder.internal.nan(1,numNN, outClass);
                    end
                end
            end
        end

        numDiff = numNN1 - numNN;
        if numDiff > 0 && NONFINITES % The number of X points without NaN is smaller than numNN
            assert(numDiff <= nX);
            nanIdx = obj.wasnanIdx(1:numDiff);
            if ~includeTies % No need to check ties
                idx(:,numNN+1:numNN1) = cast(repmat(nanIdx,[nY,1]), indexType);
                nanUpLim = min(sum(wasNaNY),nY);
                if any(wasNaNY)
                    allInds = 1:nX;
                    vecToRep = cast(allInds(numNN+1:numNN1),'like', Y);
                    idx(wasNaNY,numNN+1:numNN1) = cast(repmat(vecToRep,nanUpLim,cast(1,'like', obj.X)), indexType);
                end
                if nargout > 1
                    nanTemp = NaN(size(dist));
                    dist(:,numNN+1:numNN1) = cast(nanTemp(:,numNN+1:numNN1), outClass);
                end
            else
                if nargout < 2
                    for i=1:coder.internal.indexInt(nY)
                        if ~wasNaNY(i)
                            idx{i}(numNN+1:numNN1) = cast(nanIdx, indexType);
                        end
                    end
                else
                    for i=1:coder.internal.indexInt(nY)
                        if ~wasNaNY(i)
                            idx{i}(numNN+1:numNN1) = cast(nanIdx, indexType);
                            if inMATLAB
                                dist{i}(numNN+1:numNN1) = cast(NaN(1,numDiff), outClass);
                            else
                                nanTemp = NaN(size(dist{i}));
                                dist{i}(numNN+1:numNN1) = cast(nanTemp(1,numNN+1:numNN1), outClass);
                            end
                        end
                    end
                end
            end
        end % numDiff > 0

    else % Radius search
        isOK = isscalar(radius) && isnumeric(radius) && radius >= 0 && ~isnan(radius) && ~isinf(radius);
        coder.internal.errorIf(~isOK,'stats:KDTreeSearcher:rangesearch:BadRadius');

        % Degenerate case, create cell with proper size.
        if(nX == 0|| nY== 0)
            idx = repmat({zeros(1,0, indexType)},nY,1);
            if nargout > 1
                dist = repmat({zeros(1,0,outClass)},nY,1);
            end
            return;
        end

        powRadius = coorPow(radius,minExp,minExpIsInf);
        % Initialize the index matrix
        idx = coder.nullcopy(cell(coder.internal.indexInt(nY),1));
        % Initialize the distance matrix
        dist = coder.nullcopy(cell(coder.internal.indexInt(nY),1));
        parfor j = 1:nY
            if wasNaNY(j) %The jth point in Y has NaN
                dist{j} = zeros(1,0,'like',X);
                idx{j} = zeros(1,0, indexType);
            else
                [dist1,idx1] = search_kdtree_radius(obj, X, Y(j,:), minExp, powRadius, minExpIsInf);
                dist{j} = dist1';
                % The returned idx should be double to match MATLAB
                idx{j} = cast(idx1', indexType);
            end
        end
    end

end


function pq = search_kdtree(obj, X, queryPt, numNN, expVal,numNNmax,expValIsInf)

    [nX,nDims] = size(X);
    % Find the node containing the query point
    start_node = get_starting_node(queryPt, obj.cutDim , obj.cutVal, obj.leafNode, obj.leftChild, obj.rightChild);

    pq = struct('D',zeros(0,1,'like',X),'I',uint32(zeros(0,1)));
    if coder.internal.isConst(numNNmax)
        coder.varsize('pq.D',[numNNmax,1],[1,0]);
        coder.varsize('pq.I',[numNNmax,1],[1,0]);
    elseif coder.internal.isConst(nX)
        coder.varsize('pq.D',[nX,1],[1,0]);
        coder.varsize('pq.I',[nX,1],[1,0]);
    else
        coder.varsize('pq.D');
        coder.varsize('pq.I');
        assert(size(pq.D,1) <= nX);
        assert(size(pq.I,1) <= nX);
    end

    % Search the starting node
    node_idx_this = getNodeFromArray(obj.idxAll, obj.idxDim, start_node);

    pq = search_node(X, queryPt, node_idx_this, numNN, expVal, pq, expValIsInf);

    % If found enough nearest points and the ball is within the bounds of
    % the starting node, the search is done
    if ~isempty(pq.D)
        ballIsWithinBounds = ball_within_bounds(queryPt, obj.lowerBounds(start_node,:),...
                                                obj.upperBounds(start_node,:), pq.D(end), expVal, expValIsInf);
    else
        ballIsWithinBounds = false;
    end

    if numel(pq.D) == numNN && ballIsWithinBounds
        return
    end
    nNodes = cast(size(obj.cutDim,1),'like',X);
    if coder.internal.isConst(nNodes)
        coder.varsize('nodeStack',[nNodes,1],[1,0]);
    elseif coder.internal.isConst(nX)
        coder.varsize('nodeStack',[ceil(nX/10),1],[1,0]);
    else
        coder.varsize('nodeStack',[],[1,0]);
    end

    nodeStack = cast(1, 'like', X);  % Start from the root node

    while ~isempty(nodeStack)
        assert(size(nodeStack,1) <= nNodes);
        currentNode = nodeStack(1);     % Get the next node to be visited
        nodeStack(1) = [];
        %     nodeStack = nodeStack(2:end);   % Remove the node that is being visited
        boundsOverlapBall = bounds_overlap_ball(queryPt, obj.lowerBounds(currentNode,:),...
                                                obj.upperBounds(currentNode,:), pq.D(end), nDims, expVal, expValIsInf);

        if numel(pq.D) < numNN || boundsOverlapBall
            % we haven't found enough neighbors or the current node overlaps the ball
            if ~obj.leafNode(currentNode)
                if (queryPt(obj.cutDim(currentNode)) <= obj.cutVal(currentNode))
                    % The point is on the left side of the current node
                    % push the right child and then left child, so that
                    % the left child will be visited first
                    nodeStack = [obj.leftChild(currentNode);obj.rightChild(currentNode);nodeStack]; %#ok<AGROW>
                else
                    % The point is on the right side
                    % visit the right child first and then the left child
                    nodeStack = [obj.rightChild(currentNode);obj.leftChild(currentNode);nodeStack]; %#ok<AGROW>
                end
            elseif currentNode ~= start_node
                % current node is a leaf node that we haven't visited
                node_idx_this = getNodeFromArray(obj.idxAll,obj.idxDim,currentNode);
                pq = search_node(X, queryPt, node_idx_this, numNN, expVal, pq, expValIsInf);
            end
        end
    end

end

% Search the nearest points within a fixed radius.
% This function is similar to search_kdtree
function [dist,idx] = search_kdtree_radius(obj, X, queryPt, expVal, powRadius, expValIsInf)

    [nX,nDims] = size(X);
    % Find the node containing the query point
    start_node = get_starting_node(queryPt, obj.cutDim , obj.cutVal, obj.leafNode, obj.leftChild, obj.rightChild);

    pq = struct('D',zeros(0,1,'like',X),'I',uint32(zeros(0,1)));
    if coder.internal.isConst(nX)
        coder.varsize('pq.D',[nX,1],[1,0]);
        coder.varsize('pq.I',[nX,1],[1,0]);
    else
        coder.varsize('pq.D');
        coder.varsize('pq.I');
    end

    node_idx_this = getNodeFromArray(obj.idxAll, obj.idxDim, start_node);
    pq = search_node_radius(X, queryPt, node_idx_this, powRadius, expVal, pq,expValIsInf);

    % If found enough nearest points and the ball is within the bounds of
    % the starting node, the search is done; otherwise, continue the search

    if ~isempty(pq.D)
        ballIsWithinBounds = ball_within_bounds(queryPt, obj.lowerBounds(start_node,:),...
                                                obj.upperBounds(start_node,:), powRadius, expVal,expValIsInf);
    else
        ballIsWithinBounds = false;
    end

    if ~ballIsWithinBounds
        nNodes = cast(size(obj.cutDim,1), 'like',X);
        if coder.internal.isConst(nNodes)
            coder.varsize('nodeStack',[nNodes,1],[1,0]);
        elseif coder.internal.isConst(nX)
            coder.varsize('nodeStack',[ceil(nX/10),1],[1,0]);
        else
            coder.varsize('nodeStack',[],[1,0]);
        end
        nodeStack = cast(1, 'like', X);  % Start from the root node
        while ~isempty(nodeStack)
            assert(size(nodeStack,1) <= nNodes);
            currentNode = nodeStack(1);     % Get the next node to be visited
            nodeStack(1) = [];

            boundsOverlapBall = bounds_overlap_ball(queryPt, obj.lowerBounds(currentNode,:),...
                                                    obj.upperBounds(currentNode,:), powRadius, nDims, expVal,expValIsInf);
            if boundsOverlapBall
                if  ~obj.leafNode(currentNode)
                    if (queryPt(obj.cutDim(currentNode)) <= obj.cutVal(currentNode))
                        % The point is on the left side of the current node
                        % push the right child and then left child, so that
                        % the left child will be visited first
                        nodeStack = [obj.leftChild(currentNode);obj.rightChild(currentNode);nodeStack]; %#ok<AGROW>
                    else
                        % The point is on the right side
                        % visit the right child first and then the left child
                        nodeStack = [obj.rightChild(currentNode);obj.leftChild(currentNode);nodeStack]; %#ok<AGROW>
                    end
                elseif currentNode ~= start_node
                    node_idx_this = getNodeFromArray(obj.idxAll,obj.idxDim,currentNode);
                    pq = search_node_radius(X, queryPt, node_idx_this, powRadius, expVal, pq,expValIsInf);
                end
            end
        end
    end

    [distSorted, sortedIDx] = sort(pq.D);
    pq.D = distSorted;
    pq.I = pq.I(sortedIDx);

    dist = coorRoot(pq.D,expVal,expValIsInf);
    idx = pq.I;
end

% Search the K nearest points in X. The output will include the ties if the
% Kth neighbor has any ties
function [dist,idx] = search_kdtree_withties(obj, X, queryPt, numNN, expVal, thisK, numNNmax, expValIsInf)
    pq = search_kdtree(obj, X, queryPt, thisK, expVal,numNNmax, expValIsInf);

    len = coder.internal.indexInt(0);
    for c = numNN+1:thisK
        if pq.D(c) > pq.D(numNN)
            len = c;
            break
        end
    end

    if len ~= 0
        dist = coorRoot(pq.D(1:len-1),expVal, expValIsInf);
        idx = pq.I(1:len-1);

    else
        [dist,idx] = search_kdtree_radius(obj, X, queryPt, expVal, pq.D(numNN), expValIsInf);
    end

end


function boundsOverlapBall = bounds_overlap_ball(queryPt, lowBounds,...
                                                 upBounds, radius, nDims, expVal, expValIsInf)

    boundsOverlapBall = true;
    sumDist = zeros(1,'like',queryPt);
    for c = 1:nDims
        if queryPt(c) < lowBounds(c)
            distToAdd = coorPow(queryPt(c) - lowBounds(c), expVal, expValIsInf);
            sumDist = accuDist([sumDist,distToAdd],expVal, expValIsInf);
        elseif queryPt(c) > upBounds(c)
            distToAdd = coorPow(queryPt(c) - upBounds(c), expVal, expValIsInf);
            sumDist = accuDist([sumDist,distToAdd],expVal, expValIsInf);
        end
        if sumDist > radius
            boundsOverlapBall = false;
            return
        end
    end

end


function ballIsWithinBounds = ball_within_bounds(queryPt, lowBounds,...
                                                 upBounds, poweredRadius, expVal,expValIsInf)

    lowDist = coorPow(queryPt - lowBounds, expVal,expValIsInf);
    upDist = coorPow(queryPt - upBounds, expVal,expValIsInf);

    if isempty(queryPt) || min(lowDist) <= poweredRadius || min(upDist) <= poweredRadius
        ballIsWithinBounds = false;
        return
    end
    ballIsWithinBounds = true;

end

function pq = search_node(X, queryPt, node_idx_start, numNN, expVal,pq,expValIsInf)

    diffAllDim = bsxfun(@minus, X(node_idx_start,:), queryPt);

    % distInP is a column vector with the same size as node_idx_start
    distInPAll = coorPow(diffAllDim,expVal,expValIsInf);
    distInP = accuDist(distInPAll,expVal,expValIsInf);

    [distInPSorted,sortedInd] = sort(distInP);
    % pqHeight = numel(pq.D);
    nodeHeight = numel(distInPSorted);

    if isempty(pq.D)
        if nodeHeight <= numNN
            pq.D = distInPSorted;
            pq.I = node_idx_start(sortedInd);
        else
            pq.D = distInPSorted(1:numNN);
            pq.I = node_idx_start(sortedInd(1:numNN));
        end
    else
        [pq.D,pq.I] = mergeSort(pq.D,distInPSorted,pq.I,node_idx_start(sortedInd),numNN);
    end
    assert(size(pq.D,1) <= size(X,1));
    assert(size(pq.I,1) <= size(X,1));
end

function pq = search_node_radius(X, queryPt, node_idx_this, powRadius, expVal, pq, expValIsInf)
% Search all the points in the given node within or equal to the given radius
% Saves the nearest points in the vector.
% The nearest points in pq won't be sorted

    diffAllDim = bsxfun(@minus, X(node_idx_this,:), queryPt);

    % distInP is a column vector with the same size as node_idx_start
    distInPAll = coorPow(diffAllDim,expVal,expValIsInf);
    distInP = accuDist(distInPAll,expVal,expValIsInf);

    filteredInd = distInP <= powRadius;

    nD1 = numel(pq.D);
    nD2 = numel(distInP(filteredInd));

    pqTempD = zeros(size(X,1),1,'like',pq.D);
    pqTempD(1:nD1) = pq.D;

    pqTempI = zeros(size(X,1),1,'like',pq.I);
    pqTempI(1:nD1) = pq.I;
    if nD2 > 0
        pqTempD(nD1+1:nD1+nD2) = distInP(filteredInd);
        pqTempI(nD1+1:nD1+nD2) = node_idx_this(filteredInd);
    end
    pq.D = pqTempD(1:nD1+nD2,1);
    pq.I = pqTempI(1:nD1+nD2,1);
    % pq.D = [pq.D;distInP(filteredInd)];
    % pq.I = [pq.I;node_idx_this(filteredInd)];

end


function node = get_starting_node(queryPt, cutDim , cutVal, leafNode, leftChild, rightChild)
    node = cast(1,'like',leftChild);
    while ~leafNode(node)
        if queryPt(cutDim(node)) <= cutVal(node)
            node = leftChild(node);
        else
            node = rightChild(node);
        end
    end
end


function powDist = coorPow(pRadIn,expVal,expValIsInf)
    if expVal == 2
        powDist = pRadIn .* pRadIn;
    elseif expVal == 1 || isinf(expVal) || expValIsInf
        powDist = abs(pRadIn);
    else
        powDist = abs(pRadIn) .^ expVal;
    end

end

function distOut = coorRoot(powDist,expVal,expValIsInf)
    if expVal == 2
        distOut = sqrt(powDist);
    elseif expVal == 1 || isinf(expVal) || expValIsInf
        distOut = powDist;
    else
        distOut = powDist .^ (1/expVal);
    end

end

function distOut = accuDist(distIn,expVal,expValIsInf)
    if isinf(expVal) || expValIsInf
        aDistOut = max(distIn,[],2);
    else
        aDistOut = sum(distIn,2);
    end
    distOut = aDistOut(:,1);
end

function node_idx_this = getNodeFromArray(idxAll,idxDim,this_node)
    coder.varsize('node_idx_this');
    if idxDim(this_node) == 0
        node_idx_this = zeros(0,1,'like',idxAll);
    elseif this_node == 1
        node_idx_this = idxAll(1:idxDim(this_node));
    else
        nIdx = coder.internal.indexInt(sum(idxDim(1:this_node-1)));
        node_idx_this = idxAll(coder.internal.indexPlus(nIdx,1): ...
                               coder.internal.indexPlus(nIdx,coder.internal.indexInt(idxDim(this_node))));
    end

end

function [dOut,iOut] = mergeSort(D1,D2,I1,I2,N)
% Merging two sorted arrays in ascending order
    nD1 = numel(D1);
    nD2 = numel(D2);
    cD1 = coder.internal.indexInt(1);
    cD2 = coder.internal.indexInt(1);

    if nargin < 5
        uBound = nD1+nD2;
    else
        uBound = min(nD1+nD2,N);
    end

    dOut = zeros(uBound,1,'like',D1);
    iOut = zeros(uBound,1,'like',I1);

    for c = 1:coder.internal.indexInt(uBound)
        if D1(cD1) <=D2(cD2)
            dOut(c) = D1(cD1);
            iOut(c) = I1(cD1);
            cD1 = cD1 + 1;
            if cD1 > nD1    % Remaining candidates are in D2
                for cc = c+1:uBound
                    dOut(cc) = D2(cD2);
                    iOut(cc) = I2(cD2);
                    cD2 = cD2 + 1;
                end
                break
            end
        else
            dOut(c) = D2(cD2);
            iOut(c) = I2(cD2);
            cD2 = cD2 + 1;
            if cD2 > nD2    % Remaining candidates are in D1
                for cc = c+1:uBound
                    dOut(cc) = D1(cD1);
                    iOut(cc) = I1(cD1);
                    cD1 = cD1 + 1;
                end
                break
            end
        end
    end
    assert(size(dOut,1) <= N);
    assert(size(iOut,1) <= N);
end

function [idx, dist] = kdsearchfunMEX(obj,Y,numNNin,includeTies,minExp,radius)

    if coder.target('MATLAB')
        inMATLAB = true;
    else
        inMATLAB = false;
    end

    nX = size(obj.X,1);
    nY = size(Y,1);

    outClass = superiorfloat(obj.X,Y);

    Y = cast(Y,outClass);
    if ~strcmpi(outClass, class(obj.X))
        %only happens when X is double and Y is single
        X = cast(obj.X, outClass);
    else
        X = obj.X;
    end

    numNN1 = coder.internal.indexInt(min(numNNin,nX));
    
    % If obj.nx_nonan is a single, this call fails, so cast to double
    numNN = coder.internal.indexInt(min(numNN1, double(obj.nx_nonan)));

    wasNaNY = any(isnan(Y),2);
    c = cast(1,'like',obj.idxDim);
    numCell = numel(obj.idxDim);
    idxCell = cell(numCell,1);
    idxArray = double(obj.idxAll);
    for i = 1:numCell
        idxCellOne = idxArray(c:c+obj.idxDim(i)-1);
        if~isempty(idxCellOne)
            idxCell{i} = idxCellOne';
        else
            idxCell{i} = [];
        end
        c = c + obj.idxDim(i);
    end

    doSort = true; % This is a new PV pair in 19a (on the MATLAB-side only)
    coder.extrinsic('internal.stats.KDTreeSearcher.knnsearchmex');

    if ~isempty(numNNin) %knnsearch
                         % Degenerate case, just return an empty matrix or cell of the proper size.
        if (nY == 0 || numNN1 ==0)
            if ~includeTies
                idx = zeros(nY, numNN1);
                dist = zeros(nY, numNN1, 'like', Y);
            else
                idx = repmat({zeros(1,0)},nY,1);
                dist = repmat({zeros(1,0,'like', Y)},nY,1);
            end
            return;
        end

        if numNN > 0
            if nargout < 2
                % If some of the inputs are single, this errors, so cast to double
                idxt = internal.stats.KDTreeSearcher.knnsearchmex(X', Y', numNN, minExp, double(obj.cutDim), double(obj.cutVal), ...
                                                                  double(obj.lowerBounds'), double(obj.upperBounds'),double(obj.leftChild), double(obj.rightChild), ...
                                                                  obj.leafNode,idxCell,double(obj.nx_nonan),wasNaNY, includeTies,[],doSort);
            else

                [idxt, distt]= internal.stats.KDTreeSearcher.knnsearchmex(X', Y',numNN,minExp, double(obj.cutDim), double(obj.cutVal), ...
                                                                  double(obj.lowerBounds'), double(obj.upperBounds'),double(obj.leftChild), double(obj.rightChild), ...
                                                                  obj.leafNode, idxCell,double(obj.nx_nonan),wasNaNY,includeTies,[],doSort);

            end
            if  ~includeTies % Numeric array: transpose
                             % Initialize outputs to native types
                idxTemp = zeros(nY,numNN1);
                idxTemp(:) = idxt';
                if nargout > 1
                    distTemp = zeros(nY,numNN1,'like', Y);
                    distTemp(:) = distt';
                end
            else % Cell array
                idxTemp = coder.nullcopy(cell(coder.internal.indexInt(nY),1)); %#ok<NASGU>
                idxTemp = idxt;
                if nargout > 1
                    distTemp = coder.nullcopy(cell(coder.internal.indexInt(nY),1)); %#ok<NASGU>
                    distTemp = distt;
                end
            end
        else %numNN ==0
            if  ~includeTies
                idxTemp = zeros(nY,0);
                if nargout > 1
                    distTemp = zeros(nY,0,'like', Y);
                end
            else
                idxTemp = coder.nullcopy(cell(coder.internal.indexInt(nY),1));
                if nargout > 1
                    distTemp = coder.nullcopy(cell(coder.internal.indexInt(nY),1));
                end
            end
        end

        if includeTies && any(wasNaNY)
            %Process the point in Y having NaN values
            for c = 1:numel(wasNaNY)
                if wasNaNY(c)
                    idxTemp{c} = double(1:numNN);
                end
            end
            if coder.internal.isConst(nX)
                coder.varsize('idxTemp{:}',[1,nX],[1,1]);
            else
                coder.varsize('idxTemp{:}',[1,Inf],[1,1]);
            end
            if nargout > 1
                for c = 1:numel(wasNaNY)
                    if wasNaNY(c)
                        distTemp{c} = coder.internal.nan(1,numNN, outClass);
                    end
                end
                if coder.internal.isConst(nX)
                    coder.varsize('distTemp{:}',[1,nX],[1,1]);
                else
                    coder.varsize('distTemp{:}',[1,Inf],[1,1]);
                end
            end

        end

        numDiff = numNN1 - numNN;
        if numDiff > 0 % The number of X points without NaN is smaller than numNN
            nanIdx = obj.wasnanIdx(1:numDiff);
            if ~includeTies % No need to check ties
                idxTemp(:,numNN+1:numNN1) = repmat(nanIdx,[nY,1]);
                nanUpLim = min(sum(wasNaNY),nY);
                if any(wasNaNY)
                    allInds = 1:nX;
                    vecToRep = cast(allInds(numNN+1:numNN1),'like', Y);
                    idxTemp(wasNaNY,numNN+1:numNN1) = repmat(vecToRep,nanUpLim,1);
                end
                if nargout > 1
                    nanTemp = NaN(size(distTemp));
                    distTemp(:,numNN+1:numNN1) = nanTemp(:,numNN+1:numNN1);
                end
            else % Need to check ties
                if nargout < 2
                    for i=1:coder.internal.indexInt(nY)
                        if ~wasNaNY(i)
                            idxTemp{i}(numNN+1:numNN1) = nanIdx;
                        end
                    end
                else
                    for i=1:coder.internal.indexInt(nY)
                        if ~wasNaNY(i)
                            idxTemp{i}(numNN+1:numNN1) = nanIdx;
                            if inMATLAB
                                distTemp{i}(numNN+1:numNN1) = cast(NaN(1,numDiff), outClass);
                            else
                                nanTemp = NaN(size(distTemp{i}));
                                distTemp{i}(numNN+1:numNN1) = cast(nanTemp(1,numNN+1:numNN1), outClass);
                            end
                        end
                    end
                end
            end
        end % numDiff > 0
        if ~includeTies
            idx = idxTemp;
            if nargout > 1
                dist = distTemp;
            end
        else
            idx = coder.nullcopy(cell(coder.internal.indexInt(nY),1));
            for i=1:nY
                if isempty(idxTemp{i})
                    idx{i} = zeros(1,0);
                    if coder.internal.isConst(nX)
                        coder.varsize('idx{:}',[1,nX],[0,1]);
                    else
                        coder.varsize('idx{:}',[1,Inf],[0,1]);
                    end
                else
                    idx{i} = idxTemp{i};
                end

            end
            if nargout > 1
                dist = coder.nullcopy(cell(coder.internal.indexInt(nY),1));
                for i=1:nY
                    if isempty(idxTemp{i})
                        dist{i} = zeros(1,0,'like',Y);
                        if coder.internal.isConst(nX)
                            coder.varsize('dist{:}',[1,nX],[0,1]);
                        else
                            coder.varsize('dist{:}',[1,Inf],[0,1]);
                        end
                    else
                        dist{i} = distTemp{i};
                    end

                end
            end
        end
    else % Radius search (rangesearch)
        isOK = isscalar(radius) && isnumeric(radius) && radius >= 0 && ~isnan(radius) && ~isinf(radius);
        coder.internal.errorIf(~isOK,'stats:KDTreeSearcher:rangesearch:BadRadius');

        if nargout < 2
            if nX ~= 0 && nY ~= 0
                % If some of the inputs are single, this errors, so cast to double
                idx = internal.stats.KDTreeSearcher.knnsearchmex(X', Y', [], minExp, double(obj.cutDim), double(obj.cutVal), ...
                                                                 double(obj.lowerBounds'), double(obj.upperBounds'),double(obj.leftChild), double(obj.rightChild), ...
                                                                 obj.leafNode,idxCell,double(obj.nx_nonan),wasNaNY,includeTies,radius,doSort);
            else
                % Degenerate case, create cell with proper size.
                idx = repmat({zeros(1,0)},nY,1);
            end
        else
            if nX ~= 0 && nY ~= 0
                [idx, dist]= internal.stats.KDTreeSearcher.knnsearchmex(X', Y',[],minExp, double(obj.cutDim), double(obj.cutVal), ...
                                                                  double(obj.lowerBounds'), double(obj.upperBounds'),double(obj.leftChild), double(obj.rightChild), ...
                                                                  obj.leafNode,idxCell,double(obj.nx_nonan),wasNaNY,includeTies,radius,doSort);
            else
                % Degenerate case, create cell with proper size.
                idx = repmat({zeros(1,0)},nY,1);
                dist = repmat({zeros(1,0,outClass)},nY,1);
            end

        end
    end

end

function outType = indexType
    if coder.target('MATLAB') || coder.target('MEX')
        outType = 'double';
    else
        outType = 'int32'; % Using int32 directly because Simulink cannot simulate with coder.internal.indexIntClass returned as outputs
    end
end
