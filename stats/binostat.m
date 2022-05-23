function [m,v] = binostat(n,p)
% BINOSTAT Mean and variance of the binomial distribution.
%   [M, V] = BINOSTAT(N,P) returns the mean and variance of the
%   binomial distribution with parameters N and P.
%
%   See also BINOCDF, BINOFIT, BINOINV, BINOPDF, BINORND.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.20.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin < 2, 
    error(message('stats:binostat:TooFewInputs')); 
end

[errorcode, n, p] = distchck(2,n,p);

if errorcode > 0
    error(message('stats:binostat:InputSizeMismatch'));
end

if ~isfloat(n)
   n = double(n);
end

m = n .* p;
v = n .* p .* (1 - p);

% Return NaN for parameter values outside their respective limits.
k = (p<0 | p>1 | n<0 | round(n)~=n);
if any(k(:))
    m(k) = NaN;
    v(k) = NaN;
end
