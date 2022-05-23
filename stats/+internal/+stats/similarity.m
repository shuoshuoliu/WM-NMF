function S = similarity(X,varargin)
%SIMILARITY Construct a similarity matrix.
%   S = SIMILARITY(X, K) constructs an N-by-N symmetric similarity matrix S.
%   Rows of X correspond to observations, columns correspond to variables.
%
%   S = SIMILARITY(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional
%   parameter name/value pairs to control the algorithm used by SPECTRAL.
%   Parameters are:
%
%   'Distance'      -   A distance metric which can be any of the distance 
%                       measures accepted by the KNNSEARCH function. The 
%                       default is 'euclidean'. For more information on 
%                       KNNSEARCH and available distances, type HELP
%                       KNNSEARCH.
%
%   'P'             -   A positive scalar indicating the exponent of Minkowski
%                       distance. This argument is only valid when 'Distance'
%                       is 'minkowski'. Default is 2.
%  
%   'Cov'           -   A positive definite matrix indicating the covariance
%                       matrix when computing the Mahalanobis distance. This
%                       argument is only valid when 'Distance' is
%                       'mahalanobis'. Default is NANCOV(X).
%  
%   'Scale'         -   A vector containing nonnegative values, with
%                       length equal to the number of columns in X. This
%                       argument is only valid when 'Distance' is
%                       'seuclidean'. Default is NANSTD(X).
%
%   'SimilarityGraph' - Type of graph similarity graph. Choices are:
%          'knn'        - Construct graph using nearest neighbors 
%                       (the default).
%          'epsilon'     - Construct graph using radius search. You must
%                       specify a value for Radius if you use this option.
%
%   'NumNeighbors' -    A positive integer, specifying the number of
%                       nearest neighbors used to construct the similarity
%                       graph. SimilarityGraph must be equal to knn (Default
%                       is log(size(X,1)))
%
%   'KNNGraphType' -    Type of nearest neighbor graph. Choices are:
%          'complete'   - Two points are said to be connected when either 1
%                       is a Nearest-Neighbor of 2 OR vice versa (the
%                       default)
%          'mutual'     - Two points are said to be connected when 1
%                       is a Nearest-Neighbor of 2 AND vice versa
%
%   'Radius'        -   A nonnegative scalar, specifying the radius of the
%                       nearest neighbors search query used to construct
%                       the similarity graph. SimilarityGraph must be equal
%                       to epsilon. You must specify a value for Radius if
%                       SimilarityGraph is equal to epsilon.
%
%   'KernelScale'   -   Either string 'auto' or positive scalar specifying
%                       the scale factor. If you pass 'auto', SPECTRALCLUSTER
%                       selects an appropriate scale factor using a
%                       heuristic procedure. The 'auto' option is supported
%                       only for 'euclidean' and 'seuclidean' distances.
%                       Default: 1
%                       NOTE: The heuristic procedure for estimation of
%                           the scale factor uses subsampling. Estimates
%                           obtained by this procedure can vary from one
%                           application of SPECTRALCLUSTER to another. Set
%                           the random number generator seed prior to
%                           calling SPECTRALCLUSTER for reproducibility.
%
%   Example:
%      % Compute the similarity matrix using the default distance metric 
%      % 'euclidean'.
%      X = [rand(20,2)+2; rand(20,2)];
%      S = internal.stats.similarity(X);
%
%   See also SPECTRALCLUSTER, KNNSEARCH.

%   Copyright 2019 The MathWorks, Inc.


[varargin{:}] = convertStringsToChars(varargin{:});

% Size of X
[N,D] = size(X);

funcName = mfilename;
% Parse inputs using inputParser
p = inputParser;
addRequired(p,'X',@(x) validateattributes(x,{'single','double'},...
    {'2d','real','nonempty'},funcName,'X'));
addParameter(p,'Distance','euclidean',@(x) validateattributes(x,{'char',...
    'string','function_handle'},{'nonempty'},funcName,'Distance'));
addParameter(p,'P',[],@(x) validateattributes(x,{'single','double'},...
    {'scalar','real','positive'},funcName,'P'));
addParameter(p,'Cov',[],@(x) validateattributes(x, {'single','double'},...
    {'2d','real','square'},funcName,'Cov'));
addParameter(p,'Scale',[],@(x) validateattributes(x, {'single','double'},...
    {'row','real','nonnegative','ncols',D},funcName,'Scale'));
addParameter(p,'KernelScale',1,@(x) validateattributes(x, {'single','double',...
    'char'},{},funcName,'KernelScale'));
addParameter(p,'NumNeighbors',[],@(x) validateattributes(x,...
    {'single','double'},{'scalar','integer','real','positive'},funcName,'NumNeighbors'));
addParameter(p,'SimilarityGraph','knn',@(x) validateattributes(x,{'char','string'},...
    {'scalartext','nonempty'},funcName,'SimilarityGraph'));
addParameter(p,'KNNGraphType','complete',@(x) validateattributes(x,...
    {'char','string'},{'scalartext','nonempty'},funcName,'KNNGraphType'));
addParameter(p,'Radius',[],@(x) validateattributes(x, {'single','double'},...
    {'real','nonsparse','nonnegative','scalar'},funcName,'Radius'));
parse(p,X,varargin{:});
p = p.Results;

% Obtain parsed outputs
dist = p.Distance;
minExp = p.P;
seucInvWgt = p.Scale;
mahaCov = p.Cov;
kscale = p.KernelScale;
numN = p.NumNeighbors;
simGraph = p.SimilarityGraph;
knnType = p.KNNGraphType;
radius = p.Radius;

if ischar(dist)
    % Allowed distances
    validDists = {'euclidean'; 'seuclidean'; 'cityblock'; 'chebychev'; ...
        'mahalanobis'; 'minkowski'; 'cosine'; 'correlation'; ...
        'spearman'; 'hamming'; 'jaccard'};
    
    % Validate strings
    dist = validatestring(dist,validDists,funcName,'Distance');
end
simGraph = validatestring(simGraph,{'knn','epsilon'},funcName,'SimilarityGraph');
knnType = validatestring(knnType,{'complete','mutual'},funcName,'KNNGraphType');

if strcmp(simGraph,'epsilon') 
    % NumN cannot be specified for epsilon graph
    if ~isempty(numN)
        error(message('stats:similarity:InvalidNumNContext'));
    end
    
    % Empty radius is not allowed for epsilon graph
    if isempty(radius)
        error(message('stats:similarity:EmptyRadius'));
    end
end

if strcmp(simGraph,'knn')
    % Radius cannot be specified for knn graph
    if ~isempty(radius)
        error(message('stats:similarity:RadiusForKNN'));
    end
    
    % NumNeighbors is log(N) by default
    if isempty(numN)
        numN = max(ceil(log(N)),1);
    end
end

% Validate optional distance arguments
checkExtraArg(dist,minExp,seucInvWgt,mahaCov);

% If KernelScale is set to 'auto', compute value based on heuristic
if ischar(kscale)
    validatestring(kscale,{'auto'},funcName,'KernelScale');
    if ~(strcmp(dist,'euclidean') || strcmp(dist,'seuclidean'))
        error(message('stats:similarity:AutoInvalidDist'));
    end
    
    if strcmp(dist,'seuclidean')
        % If the distance is seculidean standardize X before calculating
        % optimalKernelScale
        mu = mean(X,1,'omitnan');
        
        % Standardize using user specified Scale. If Scale is empty, use
        % nanstd
        if isempty(seucInvWgt)
            seucInvWgt = std(X,[],1,'omitnan');
        end
        
        Xstd = X - mu;
        nonZeroSigma = seucInvWgt>0;
        seucInvWgt(~nonZeroSigma) = 0;
        if any(nonZeroSigma)
            Xstd(:,nonZeroSigma) = Xstd(:,nonZeroSigma)./seucInvWgt(nonZeroSigma);
        end
        
        % Optimal Kernel Scale for seuclidean distance
        kscale = classreg.learning.svmutils.optimalKernelScale(Xstd,[],1);
    else
        % Optimal Kernel Scale for euclidean distance
        kscale = classreg.learning.svmutils.optimalKernelScale(X,[],1);
    end
else
    validateattributes(kscale,{'single','double'},{'scalar','real',...
        'nonempty','positive','nonnan'},funcName,'KernelScale');
end

% Combine extra distance arguments. Validation shall be done in knnsearch
% or rangesearch
extraArgs = {};
if ~isempty(minExp)
    extraArgs = [extraArgs(:)',{'P'},{minExp}];
end
if ~isempty(mahaCov)
    extraArgs = [extraArgs(:)',{'Cov'},{mahaCov}];
end
if ~isempty(seucInvWgt)
    extraArgs = [extraArgs(:)',{'Scale'},{seucInvWgt}];
end

% Construct either a KNN or Epsilon graph
if strcmp(simGraph,'knn')
    % KNN graph generated using knnsearch
    [rows,weights] = knnsearch(X,X,'K',numN+1,'Sort',false,'Distance',dist,...
        'IncludeTies',true,extraArgs{:});
    S = makeSimilarityMatrixFromIndices(rows,weights,kscale,N);
    
    % The KNN graph is not guaranteed to be symmetric. Make the matrix
    % symmetric using either 'complete' or 'mutual'
    if strcmp(knnType,'complete')
        S = max(S,S');
    else % mutual
        S = min(S,S');
    end
else
    % Epsilon graph generated using rangesearch
    [rows,weights] = rangesearch(X,X,radius,'Sort',false,'Distance',dist,...
        extraArgs{:});
    S = makeSimilarityMatrixFromIndices(rows,weights,kscale,N);
end
end

function S = makeSimilarityMatrixFromIndices(rows,weights,kscale,N)
% rows and weights are always cells (rangesearch and knnsearch-with-ties)
weights = exp(-([weights{:}]/kscale).^2);

% Get column indices
cols = ones(size(weights));
colind = 1;
for idx = 1:N
    num = length(rows{idx});
    cols(colind:colind+num-1) = idx*cols(colind:colind+num-1);
    colind = colind+num;
end
rows = [rows{:}];

% Return a sparse similarity matrix. 'sparse' does not accept weights with
% single precision. Convert weights to double.
S = sparse(rows,cols,double(weights),N,N);
end

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