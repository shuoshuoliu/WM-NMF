function [m,v]= poisstat(lambda)
%POISSTAT Mean and variance for the Poisson distribution.
%   [M,V] = POISSTAT(LAMBDA) returns the mean and variance of
%   the Poisson distribution with parameter LAMBDA.
%
%   See also POISSCDF, POISSFIT, POISSINV, POISSPDF, POISSRND.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.22.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin <  1, 
    error(message('stats:poisstat:TooFewInputs')); 
end

% Initialize mean and variance to zero.
m = zeros(size(lambda),'like',lambda);
v = zeros(size(lambda),'like',lambda);

% Lambda must be positive.
k = find(lambda <= 0);
if any(k)
    m(k) = NaN;
    v(k) = NaN;
end

k = find(lambda > 0);
if any(k)
    m(k) = lambda(k);
    v(k) = lambda(k);
end
