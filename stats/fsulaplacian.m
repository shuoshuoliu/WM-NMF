function [ranked,score] = fsulaplacian(X,varargin)
%FSULAPLACIAN Importance of features (variables) for unsupervised learning using Laplacian scores
%   IDX=FSULAPLACIAN(X) ranks features for an N-by-P matrix X with one row
%   per observation and one column per feature.
%
%   IDX is a 1-by-P vector for P variables. IDX are indices of columns in X
%   ordered by importance, meaning IDX(1) is the index of the most
%   important feature.
%
%   [IDX,SCORES]=FSULAPLACIAN(...) also returns feature scores SCORES, a
%   1-by-P array for P features. SCORES have the same order as features in
%   the input data, meaning SCORES(1) is the score for the first feature.
%   Large score indicates important feature.
%
%   [...]=FSULAPLACIAN(X,'PARAM1',val1,'PARAM2',val2,...) specifies
%   optional parameter name/value pairs:
%       'Similarity'   - Either empty or N-by-N symmetric matrix with
%                        nonnegative elements specifying pairwise
%                        similarities between observations in X. If you
%                        pass a nonempty similarity matrix, FSULAPLACIAN
%                        does not accept any other optional parameters.
%                        Default: []
%       'Distance'     - Distance metric which can be any of the distance
%                        measures accepted by the KNNSEARCH function. For
%                        more information on KNNSEARCH and available
%                        distances, type HELP KNNSEARCH. Default:
%                        'euclidean'
%       'P'            - Positive integer scalar indicating the exponent
%                        of Minkowski distance. This argument is only
%                        valid when 'Distance' is 'minkowski'. Default: 2
%       'Cov'          - Positive-definite matrix specifying the
%                        covariance matrix for computing the Mahalanobis
%                        distance. This argument is only valid when
%                        'Distance' is 'mahalanobis'. Default: nancov(X)
%       'Scale'        - Vector of nonnegative values, with length equal
%                        to the number of columns in X. This argument is
%                        only valid when 'Distance' is 'seuclidean'.
%                        Default: nanstd(X)
%       'NumNeighbors' - Positive integer specifying the number of
%                        nearest neighbors used to construct the similarity
%                        graph. Default: log(size(X,1)))
%       'KernelScale'  - Either string 'auto' or positive scalar specifying
%                        the scale factor for the Gaussian kernel. If you
%                        pass 'auto', FSULAPLACIAN selects an appropriate
%                        scale factor using a heuristic procedure. The
%                        'auto' option is supported only for 'euclidean'
%                        and 'seuclidean' distances. Default: 1
%                          NOTE: The heuristic procedure for estimation of
%                            the scale factor uses subsampling. Estimates
%                            obtained by this procedure can vary from one
%                            application of FSULAPLACIAN to another. Set
%                            the random number generator seed prior to
%                            calling FSULAPLACIAN for reproducibility.
%
%   Example:
%       % Identify important variables in the ionosphere dataset:
%       load ionosphere
%       rng(1) % for reproducibility with the 'auto' option
%       [idx,scores] = fsulaplacian(X,'KernelScale','auto');
%       bar(scores(idx))
%       xlabel('Feature rank')
%       ylabel('Feature importance score')

%   Copyright 2019 The MathWorks, Inc.

[varargin{:}] = convertStringsToChars(varargin{:});

% Check X: Support sparse, do not support integer and complex
internal.stats.checkSupportedNumeric('X',X,false,true,false);
[N,J] = size(X);

% Prepare output
ranked = zeros(1,J);
score = zeros(1,J,'like',X);

% Get similarity matrix if passed
args = {'similarity'};
defs = {          []};
[S,~,extraArgs] = internal.stats.parseArgs(args,defs,varargin{:});

% Good features
goodColIdx = 1:J;
Jgood = J;

% Handle nans
if any(isnan(X(:)))
    if ~isempty(S)
        error(message('stats:fsulaplacian:NaNsInXWithSimilarityMatrix'));
    end
    
    % Remove features that have nothing but NaN's
    badcols = all(isnan(X),1);
    X(:,badcols) = [];
    Jgood = J - sum(badcols);
    ranked(Jgood+1:end) = find(badcols);
    goodColIdx = find(~badcols);
    
    % Remove rows with NaN's
    badrows = any(isnan(X),2);
    X(badrows,:) = [];
    if isempty(X)
        error(message('stats:fsulaplacian:NoDataAfterNaNsRemoved'));
    end
end

if isempty(S)
    S = internal.stats.similarity(X,extraArgs{:},...
        'SimilarityGraph','knn','KNNGraphType','complete');
else
    validateattributes(S,{'single','double'},...
        {'real','square','nonnegative','nonempty','size',[N N]},...
        'fsulaplacian','Similarity');

    if ~issymmetric(S)
        error(message('stats:spectraleigs:InvalidSimMat'));
    end
    
    if ~isempty(extraArgs)
        error(message('stats:fsulaplacian:ExtraParamsWithSimilarityMatrix'));
    end
end

%
% Scores for good features. (For bad features, scores are zero.)
%

d = full(sum(S,1))';
sumd = sum(d);
maxd = max(d);

X = X - (d'*X)/sumd;

for j=1:numel(goodColIdx)
    xj = double(X(:,j));
    denom = xj'*(d.*xj);
    if denom < eps(maxd)
        score(goodColIdx(j)) = 0;
    else
        score(goodColIdx(j)) = min( 1, (xj'*S*xj)/denom );
    end
end

% Ranks for good features

[~,goodRanked] = sort(score(goodColIdx),'descend');
ranked(1:Jgood) = goodColIdx(goodRanked);

end
