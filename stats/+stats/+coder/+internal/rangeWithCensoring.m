function [w,maxx,minx] = rangeWithCensoring(x , censoring) %#codegen %#internal
% Code generation helper for calculating range with censored values, used
% by the distribution fitting functions which in turn are used by fitdist

%   Copyright 2019 The MathWorks, Inc.

tempMin = coder.internal.inf; % use inf/realmax
tempMax = -coder.internal.inf; % use inf/realmax
for idx = 1:numel(x)
    if isempty(censoring) || (~isempty(censoring) && ~censoring(idx)) % only use uncensored values in determining
        %             x(idx) = coder.internal.inf;
        if x(idx) > tempMax
            tempMax = x(idx);
        end
        if x(idx) < tempMin
            tempMin = x(idx);
        end
    end
end
maxx = tempMax;
minx = tempMin;
w = maxx - minx;



