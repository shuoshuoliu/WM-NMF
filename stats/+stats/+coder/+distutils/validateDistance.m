function dist = validateDistance(distance)
%#codegen
%VALIDATEDISTANCE helper function to validate Distance input for
%   knnsearch and rangesearch.

%   Copyright 2017 The MathWorks, Inc.

coder.inline('always');
validateattributes(distance,{'char','string'},{'nonempty','row'},mfilename,'Distance');
methods = { 'minkowski','seuclidean', 'mahalanobis','euclidean','cityblock',...
    'chebychev', 'cosine', 'correlation','spearman', 'hamming', 'jaccard','squaredeuclidean'};
if strcmpi(distance,'chebyshev')
    distancetemp = 'chebychev';
else
    distancetemp = distance;
end

dist = validatestring(distancetemp,methods,mfilename, 'Distance');
end


