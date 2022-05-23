function parmci=statexpci(parmhat,cv,alpha,x,cens,freq) %#codegen %#internal
% Code generation compatible statexpci, a helper for distribution fitting functions

%   Copyright 2019 The MathWorks, Inc.

%STATSEXPCI Exponential parameter confidence interval.



% Number of observations
if isempty(freq)
    if isvector(x)
        n = length(x);
    else
        n = size(x,1);
    end
else
    n = sum(freq);
end

% Number not censored
if ~isempty(cens)
    if ~isempty(freq)
        nunc = n - sum(freq.*cens);
    else
        nunc = n - sum(cens);
    end
else
    nunc = n;
end

% Confidence interval
if nunc > 0
    gamcrit = gaminv([alpha/2 1-alpha/2], nunc, 1);
    parmci = [nunc.*parmhat ./ gamcrit(2);
              nunc.*parmhat./ gamcrit(1)];
else
    parmci = coder.internal.nan(2,numel(parmhat),'like',x);
end
