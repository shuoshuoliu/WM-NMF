function [idx,dist]=knnsearch(obj,Y,varargin)
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
%                   smallest distances in D.
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
%    'SortIndices'  A flag to indicate if output distances and the
%                   corresponding indices should be sorted in the order of 
%                   distances ranging from the smallest to the largest 
%                   distance.  Default is true.
% 
%   Example:
%      % Create an ExhaustiveSearcher object for data X with the cosine 
%      % distance. The kd-tree method does not support nearest neighbors 
%      % search for the cosine distance, therefore CREATENS creates an
%      % ExhaustiveSearcher object.
%      X = randn(100,5);
%      Y = randn(25,5);
%      ns = createns(X, 'distance', 'cosine');
%
%      % Find 5 nearest neighbors in X for each point in Y
%      idx = knnsearch(ns,Y,'k',5);
%
%      % Choose a different distance metric from the one used when creating
%      % this ExhaustiveSearcher object
%      idx2 = knnsearch(ns,Y,'k',5, 'distance','correlation');
%
%   See also KNNSEARCH, ExhaustiveSearcher, CREATENS, RANGESEARCH,
%   KDTreeSearcher, PDIST2.

%   Copyright 2011-2018 The MathWorks, Inc.


[NX,nDims] = size(obj.X);
[NY,nDims2]= size(Y);

if nDims2 ~= nDims
    error(message('stats:ExhaustiveSearcher:knnsearch:SizeMisMatch', nDims));
end

[varargin{:}] = convertStringsToChars(varargin{:});
pnames = { 'k'  'distance'  'p', 'cov', 'scale','includeties' 'sortindices'};
dflts =  { 1    []           []   []     []    false            true};
[numNN,distMetric, minExp, mahaCov, seucInvWgt, includeTies, doSort] = ...
     internal.stats.parseArgs(pnames, dflts, varargin{:});

validateattributes(doSort,{'logical','numeric'},{'scalar'},'','sortindices');
doSort = logical(doSort);
 
if  ~islogical(includeTies)
     error(message('stats:ExhaustiveSearcher:knnsearch:InvalidIncludeTies'));
 end
 
if isempty(distMetric)
      distMetric = obj.Distance;
elseif ischar(distMetric)
     methods = {'euclidean'; 'seuclidean'; 'cityblock'; 'chebychev'; ...
        'mahalanobis'; 'minkowski'; 'cosine'; 'correlation'; ...
        'spearman'; 'hamming'; 'jaccard'};
    i =  find(strncmpi(distMetric, methods, length(distMetric)));
    if length(i) > 1
        error(message('stats:ExhaustiveSearcher:knnsearch:BadDistance', distMetric));
    elseif isempty(i)
        error(message('stats:ExhaustiveSearcher:knnsearch:UnrecognizedDistance', distMetric));
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

if ~logical(includeTies)

       [dist,idx] = pdist2(obj.X,Y, distMetric, arg{:}, 'smallest',numNN, 'sortindices', doSort);
       
       if obj.checkNegativeDistance  && isa(distMetric,'function_handle') && any(dist(:)<0)
            error(message('stats:ExhaustiveSearcher:knnsearch:NegativeDistanceValue'));
       end  
       dist = dist';
    
    idx = idx';
else %need to include ties if the Kth neighbor has any ties
    
    if (numNN < 3)
        extra = 2;     
    elseif numNN < 5
        extra = numNN;
    else
        extra = 5;
    end
    K=numNN+extra;
    
    [dist2, idx2]=pdist2(obj.X,Y, distMetric, arg{:}, 'smallest',K, 'sortindices', doSort);
    if obj.checkNegativeDistance && isa(distMetric,'function_handle') && any(dist2(:)<0)
         error(message('stats:ExhaustiveSearcher:knnsearch:NegativeDistanceValue'));
    end
    
    dist2=dist2';
    idx2=idx2';
    outClass = class(dist2);
    %  Degenerate case, just return an empty cell of the proper size.
    if (NY == 0 || NX ==0)
        idx =repmat({zeros(1,0)},NY,1);
        dist=repmat({zeros(1,0,outClass)},NY,1);
        return;
    end
    
    dist=cell(NY,1);
    idx = cell(NY,1);
    
    %mat2cell and num2cell are slower than using the loops when NY is
    %much larger than numNN
    if numNN >= NX 
        for i =1 : NY
            idx{i}= idx2(i, :);
            dist{i} = dist2(i,:);
        end
        return;
    end
    

    %Expect that most query points don't have extra tie neighbors
    done=(dist2(:,numNN+1)>dist2(:,numNN)) | isnan(dist2(:,numNN+1));
           
    doneIdx=find(done);
   
    %mat2cell and num2cell are slower than the loops when the number of query
    %points is much larger than the value of K
    doneLen=numel(doneIdx);
    for i=1:doneLen
        tempIdx=doneIdx(i);
        idx{tempIdx}=idx2(tempIdx,1:numNN);
        dist{tempIdx}=dist2(tempIdx,1:numNN);
    end
      
    notDone=~done;
    notDoneIdx=find(notDone);
    dist2= dist2(notDoneIdx,:);
    idx2 = idx2(notDoneIdx,:);
      
    if K>=NX %we have scanned all the points in X
        for i=1:numel(notDoneIdx)
            tempT=find(dist2(i,numNN+1:end)>dist2(i,numNN),1);
            if ~isempty(tempT)
                N=numNN+tempT-1;
            else
                tempT = find(isnan(dist2(i,numNN+1:end)));
                if ~isempty(tempT)
                    N= numNN+tempT-1;
                else
                    N = size(dist2,2);
                end
            end
            dist{notDoneIdx(i)}=dist2(i,1:N);
            idx{notDoneIdx(i)}= idx2(i,1:N);  
        end
               
    else
        for i=1:numel(notDoneIdx)
            tempT=find(dist2(i,numNN+1:end)>dist2(i,numNN),1);
            tempIdx = notDoneIdx(i);
            if ~isempty(tempT)
                N=numNN+tempT-1;
                dist{tempIdx}=dist2(i,1:N);
                idx{tempIdx}= idx2(i,1:N);            
            else
                if isinf(dist2(i,numNN))
                   r =  dist2(i,numNN);
                else
                   r = dist2(i,numNN)+10*eps(dist2(i,numNN));
                end
                [dTemp, iTemp] = pdist2(obj.X,Y(notDoneIdx(i),:),...
                    distMetric, arg{:}, 'radius',r, 'sortindices', doSort);
                dist(tempIdx) = dTemp;
                idx(tempIdx) = iTemp;
                
            end
        end 
       
    end
   
end

end %method knnsearch

%Give an error if an extra input is provided
function checkExtraArg(distMetric,minExp,seucInvWgt,mahaCov)
if ~isempty(minExp) && ~strncmp(distMetric,'min',3)
    error(message('stats:ExhaustiveSearcher:knnsearch:InvalidMinExp'));
end
if ~isempty(seucInvWgt) && ~strncmp(distMetric,'seu',3)
    error(message('stats:ExhaustiveSearcher:knnsearch:InvalidSeucInvWgt'));
end
if ~isempty(mahaCov) && ~strncmp(distMetric,'mah',3)
    error(message('stats:ExhaustiveSearcher:knnsearch:InvalidMahaCov'));
end

end


