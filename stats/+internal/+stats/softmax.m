function Y = softmax(X, Dim)
%   Y = softmax(X), for matrix X, returns matrix Y which is the softmax of
% each column of X, i.e., bsxfun(@rdivide,exp(X),sum(exp(X))). Uses
% logsumexp to avoid over/underflow.
%   Y = softmax(X,Dim) returns the softmax along dimension Dim. softmax(X)
% is the same as softmax(X,1)

%   Copyright 2016 The MathWorks, Inc.

if nargin < 2
    Dim = 1;
end

Y = exp(bsxfun(@minus, X, internal.stats.logsumexp(X,Dim)));

% Handle Infs correctly. The result is what you would get if you had finite
% X components, and then let them go to infinity uniformly. So, if you have
% any +Inf components, each of the +Infs maps to 1/NumPositiveInfs, and all
% other components map to zero. If you have all -Infs, all components map
% to 1/NumComponents.    
M = max(X,[],Dim,'includenan');
NegInfIdx = isneginf(M);
PosInfIdx = isposinf(M);
if Dim==1
    Y(:,NegInfIdx) = 1/size(X,1);
    IsInf = isposinf(X(:,PosInfIdx));
    Y(:,PosInfIdx) = IsInf./sum(IsInf);
else
    Y(NegInfIdx,:) = 1/size(X,2);
    IsInf = isposinf(X(PosInfIdx,:));
    Y(PosInfIdx,:) = IsInf./sum(IsInf,2);
end
end

function Y = isposinf(X)
Y = isinf(X) & X>0;
end

function Y = isneginf(X)
Y = isinf(X) & X<0;
end
