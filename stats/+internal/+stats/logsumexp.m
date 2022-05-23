function L = logsumexp(X, Dim)
%   L = logsumexp(X), for matrix X, returns log(sum(exp(X))) while avoiding
% over/underflow.
%   L = logsumexp(X,Dim) returns log(sum(exp(X),Dim)). logsumexp(X) is the
% same as logsumexp(X,1)

%   Copyright 2016 The MathWorks, Inc.

if nargin < 2
    Dim = 1;
end

M = max(X,[],Dim,'includenan');
L = M + log(sum(exp(bsxfun(@minus, X, M)),Dim));
% Now handle Infs properly. If max is Inf, return Inf. If max is -Inf,
% return -Inf:
InfRows = isinf(M);
L(InfRows) = M(InfRows);
end
