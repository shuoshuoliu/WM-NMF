function [idx,corepts] = dbscan(X,epsilon,minpts,varargin)
%DBSCAN Density-Based algorithm for clustering
%   IDX = DBSCAN(X, EPSILON, MINPTS) partitions the points in the N-by-P
%   data matrix X into clusters based on parameters EPSILON and MINPTS.
%   EPSILON is a threshold for a neighborhood search query. MINPTS is a
%   positive integer used as a threshold to determine whether a point is a
%   core point. IDX is an N-by-1 vector containing cluster indices. An
%   index equal to '-1' implies a noise point.
%
%   IDX = DBSCAN(D, EPSILON, MINPTS, 'DISTANCE', 'PRECOMPUTED') is an
%   alternative syntax that accepts distances D between pairs of
%   observations instead of raw data. D may be a vector or matrix as
%   computed by PDIST or PDIST2, or a more general dissimilarity vector or
%   matrix conforming to the output format of PDIST or PDIST2.
%
%   [IDX, COREPTS] = DBSCAN(...) returns a logical vector COREPTS
%   indicating indices of core-points as identified by DBSCAN.
%
%   IDX = DBSCAN(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional
%   parameter name/value pairs to control the algorithm used by DBSCAN.
%   Parameters are:
%
%   'Distance'      -   a distance metric which can be any of the distance 
%                       measures accepted by the PDIST2 function. The 
%                       default is 'euclidean'. For more information on 
%                       PDIST2 and available distances, type HELP PDIST2. 
%                       An additional choice is:
%     'precomputed' -   Needs to be specified when a custom distance matrix
%                       is passed in
%
%    'P'            -   A positive scalar indicating the exponent of Minkowski
%                       distance. This argument is only valid when 'Distance'
%                       is 'minkowski'. Default is 2.
%  
%    'Cov'          -   A positive definite matrix indicating the covariance
%                       matrix when computing the Mahalanobis distance. This
%                       argument is only valid when 'Distance' is
%                       'mahalanobis'. Default is NANCOV(X).
%  
%    'Scale'        -   A vector S containing non-negative values, with length
%                       equal to the number of columns in X. Each coordinate
%                       difference between X and a query point is scaled by the
%                       corresponding element of S. This argument is only valid
%                       when 'Distance' is 'seuclidean'. Default is NANSTD(X).
%
%   Example:
%      % Find clusters in data X, using the default distance metric 
%      % 'euclidean'.
%      X = [rand(20,2)+2; rand(20,2)];
%      idx = dbscan(X,0.5,2);
%
%   See also KMEANS, KMEDOIDS, PDIST2, PDIST.

%   Copyright 2018 The MathWorks, Inc.

% Convert strings to chars
[varargin{:}] = convertStringsToChars(varargin{:});

% Parse X
validateattributes(X,{'single','double','logical'},{'2d','nonsparse','real',...
    'nonempty'},'','X');

% Parse epsilon
validateattributes(epsilon,{'single','double'},{'real',...
    'nonsparse','nonnegative'},'','epsilon');

% Epsilon must be [] when X is logical and must be non-empty when X is
% numeric
if islogical(X) && ~isempty(epsilon)
    error(message('stats:dbscan:EpsilonEmpty')); % Epsilon must be empty
elseif isnumeric(X) && isempty(epsilon)
    error(message('stats:dbscan:EpsilonNonEmpty')); % Epsilon must not be empty
elseif ~isscalar(epsilon) && ~isempty(epsilon)
    error(message('stats:dbscan:EpsilonScalar'));
end

% Parse minpts
validateattributes(minpts,{'numeric'},{'scalar','real','nonempty',...
    'nonsparse','integer','positive'},'','minpts');

% Store size of X
[nobs,ndims] = size(X);

% Parse Name-Value pairs
pnames = {'distance','p','scale','cov'};
dflts =  {'euc',[],[],[]};
[dist,p,scale,cov] = internal.stats.parseArgs(pnames, dflts, varargin{:});

% ------ The following lines parse parameters that shall be passed to
% pdist2 from the C++ routine. We parse once here instead of parsing
% multiple times within pdist2. These lines are similar but not identical
% to the lines in pdist2. ------ %
additionalArg = [];
isFcnHndl = false;
if ischar(dist)
    methods = {'euclidean'; 'seuclidean'; 'cityblock'; 'chebychev'; ...
        'mahalanobis'; 'minkowski'; 'cosine'; 'correlation'; ...
        'spearman'; 'hamming'; 'jaccard'; 'squaredeuclidean';'precomputed'};
    i = find(strncmpi(dist,methods,length(dist)));
    if length(i) > 1
        error(message('stats:pdist2:AmbiguousDistance', dist));
    elseif isempty(i)
        error(message('stats:pdist2:UnrecognizedDistance', dist));
    else
        if i == 12 %'squaredeuclidean'
            dist = 'sqe'; % represent squared Euclidean
        else
            dist = methods{i}(1:3);
        end
        
        if ~isempty(p) || ~isempty(cov) || ~isempty(scale)
            
            % Can only specify one of Cov, P, or Scale
            if sum(~isempty(p) + ~isempty(cov) + ~isempty(scale)) > 1
                error(message('stats:dbscan:SingleOptArg'));
            end
            
            if ~isempty(p)
                if ~strcmp(dist,'min')
                    error(message('stats:dbscan:InvalidOptArg',...
                        'P', 'Minkowski'));
                end
                validateattributes(p,{'numeric'},{'scalar','real','positive',...
                    'nonempty','nonsparse'},'','P');
                arg = p;
            elseif ~isempty(cov)
                if ~strcmp(dist,'mah')
                    error(message('stats:dbscan:InvalidOptArg',...
                        'Cov', 'Mahalanobis'));
                end
                validateattributes(cov,{'numeric'},{'2d','real',...
                    'nonempty','nonsparse'},'','Cov');
                arg = cov;
            else
                if ~strcmp(dist,'seu')
                    error(message('stats:dbscan:InvalidOptArg',...
                        'Scale','Seuclidean'));
                end
                validateattributes(scale,{'numeric'},{'vector','real',...
                    'positive','nonempty','nonsparse'},'','Scale');                
                arg = scale;
            end
            
            % Get the additional distance argument from the inputs
            if isnumeric(arg)
                switch dist
                    case {'seu' 'mah' 'min'}
                        additionalArg = arg;
                end
            end
        end
    end
elseif isa(dist, 'function_handle')
    distfun = dist;
    dist = 'usr';
    isFcnHndl = true;
else
    error(message('stats:pdist2:BadDistance'));
end

switch dist
    case 'seu' % Standardized Euclidean weights by coordinate variance
        if isempty(additionalArg)
            additionalArg =  nanvar(X,[],1);
            if any(additionalArg == 0)
                warning(message('stats:pdist2:ConstantColumns'));
            end
            additionalArg = 1./ additionalArg;
        else
            if ~(isvector(additionalArg) && length(additionalArg) == ndims...
                    && all(additionalArg >= 0))
                error(message('stats:pdist2:InvalidWeights'));
            end
            if any(additionalArg == 0)
                warning(message('stats:pdist2:ZeroInverseWeights'));
            end
            additionalArg = 1./ (additionalArg .^2);
        end
        
    case 'mah' % Mahalanobis
        if isempty(additionalArg)
            if nobs == 1
                error(message('stats:pdist2:tooFewXRowsForMah'));
            end
            additionalArg = nancov(X);
            [T,flag] = chol(additionalArg);
        else %provide the covariance for mahalanobis
            if ~isequal(size(additionalArg),[ndims,ndims])
                error(message('stats:pdist2:InvalidCov'));
            end
            %cholcov will check whether the covariance is symmetric
            [T,flag] = cholcov(additionalArg,0);
        end
        
        if flag ~= 0
            error(message('stats:pdist2:InvalidCov'));
        end
        
        if ~issparse(X)
            additionalArg = T \ eye(ndims); %inv(T)
        end
        
    case 'min' % Minkowski
        if isempty(additionalArg)
            additionalArg = 2;
        elseif ~(isscalar(additionalArg) && additionalArg > 0)
            error(message('stats:pdist2:BadMinExp'));
        end
    case 'cos' % Cosine
        [X,flag] = normalizeX(X);
        if flag
            warning(message('stats:pdist2:ZeroPoints'));
        end
        
    case 'cor' % Correlation
        X = bsxfun(@minus,X,mean(X,2));
        [X,flag] = normalizeX(X);
        if flag
            warning(message('stats:pdist2:ConstantPoints'));
        end
        
    case 'spe'
        X = tiedrank(X')'; % treat rows as a series
        X = X - (ndims+1)/2; % subtract off the (constant) mean
        [X,flag] = normalizeX(X);
        if flag
            warning(message('stats:pdist2:TiedPoints'));
        end
        
    otherwise
        
end

% Note that if the above switch statement is modified to include the
% 'che', 'euc', or 'cit' distances, that code may need to be repeated
% in the corresponding block below.
if strcmp(dist,'min') % Minkowski distance
    if isinf(additionalArg) %the exponent is inf
        dist = 'che';
        additionalArg = [];
    elseif additionalArg == 2 %the exponent is 2
        dist = 'euc';
        additionalArg = [];
    elseif additionalArg == 1 %the exponent is 1
        dist = 'cit';
        additionalArg = [];
    end
end
% ------ End of parsing for pdist2 inputs ------ %

% If X is a logical matrix, assume that it is a distance matrix and issue a
% warning
if islogical(X) && ~strcmp(dist,'pre')
    dist = 'pre';
    warning(message('stats:dbscan:WarnPrecomputed'));
end

% The C++ routine finds neighbors based on 'method'. It has 3 options:
%   0   - User passed in NxN distance matrix Y, like the output of pdist2
%   1   - User passed in Nx(N-1)/2 dist matrix Y, like the output of pdist
%   2   - User passed in raw data X, call pdist2 to compute distances of
%         each point to find neighbors.
method = 2;

if strcmp(dist,'pre') % First input is distance matrix
    
    % DBSCAN does not need the actual distances but just an indication of
    % neighbors of a given point. Hence, if the user has passed in a 
    % numeric distance matrix, convert numeric values to logical before
    % passing to C++ routine.
    if ~islogical(X)
        X = X <= epsilon;
    end
    
    % Don't pass empty epsilon to builtin code as it expects a scalar
    % value. Epsilon is no longer needed.
    if isempty(epsilon)
        epsilon = 0;
    end
    
    % Tell the builtin routine whether we are working with pdist or pdist2
    % type of distance matrix by setting "method"
    if isvector(X)
        % Output is a lower triangular matrix similar to the output of
        % pdist
        % Estimate the number of observations. This estimation is inspired
        % by the linkage.m
        if nobs < 1
            error(message('stats:dbscan:TooFewDistances'));
        end
        
        m = ceil(sqrt(2*ndims)); % m = (1+sqrt(1+8*nobs))/2, but works for large nobs
        if nobs ~= 1 || m*(m-1)/2 ~= ndims
            error(message('stats:dbscan:BadDistanceSize'));
        end
    
        nobs = m;
        method = 1;
    else
        % Output is the full NxN distance matrix, similar to the output of
        % pdist2
        if ~isequal(nobs,ndims)
            error(message('stats:dbscan:NonSquareDistMat'));
        end
        method = 0;
    end
    
    % Call internal C++ builtin method
    if nargout > 1
        [idx,corepts] = internal.stats.dbscan([],epsilon,minpts,dist,method,X,...
            nobs,ndims,additionalArg,isFcnHndl);
    else
        idx = internal.stats.dbscan([],epsilon,minpts,dist,method,X,...
            nobs,ndims,additionalArg,isFcnHndl);
    end
elseif strcmp(dist,'usr') % dist is a user specified function handle
    if nargout > 1
        [idx,corepts] = internal.stats.dbscan(X,epsilon,minpts,distfun,...
            method,false,nobs,ndims,additionalArg,isFcnHndl);
    else
        idx = internal.stats.dbscan(X,epsilon,minpts,distfun,...
            method,false,nobs,ndims,additionalArg,isFcnHndl);
    end
else % dist is one of the in-built distance metrics
    if nargout > 1
        [idx,corepts] = internal.stats.dbscan(X',epsilon,minpts,dist,...
            method,false,nobs,ndims,additionalArg,isFcnHndl);
    else
        idx = internal.stats.dbscan(X',epsilon,minpts,dist,...
            method,false,nobs,ndims,additionalArg,isFcnHndl);
    end
end
end

%---------------------------------------------
% Normalize the data matrix X to have unit norm
function [X,flag] = normalizeX(X)
Xmax = max(abs(X),[],2);
X2 = X./Xmax;
Xnorm = sqrt(sum(X2.^2, 2));

% Find out points for which distance cannot be computed.

% The norm will be NaN for rows that are all zeros, fix that for the test
% below.
Xnorm(Xmax==0) = 0;

% The norm will be NaN for rows of X that have any +/-Inf. Those should be
% Inf, but leave them as is so those rows will not affect the test below.
% The points can't be normalized, so any distances from them will be NaN
% anyway.

% Find points that are effectively zero relative to the point with largest norm.
flag =  any(Xnorm <= eps(max(Xnorm)));
Xnorm = Xnorm .* Xmax;
X = X./Xnorm;
end