function [m,v] = betastat(a,b)
%BETASTAT Mean and variance for the beta distribution.
%   [M,V] = BETASTAT(A,B) returns the mean and variance 
%   of the beta distribution with parameters A and B.
%
%   See also BETACDF, BETAFIT, BETAINV, BETALIKE, BETAPDF, BETARND.
    
%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.33.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin < 2, 
    error(message('stats:betastat:TooFewInputs'));
end

[errorcode, a, b] = distchck(2,a,b);

if errorcode > 0
    error(message('stats:betastat:InputSizeMismatch'));
end

m = zeros(size(a),'like',internal.stats.dominantType(a,b)); % Single if a, or b is
v = m;

%   Return NaN if the parameter values are outside their respective limits.
k = (a <= 0 | b <= 0);
if any(k(:)) 
    m(k) = NaN;
    v(k) = NaN;
end

k1 = ~k;
if any(k1(:))
    m(k1) = a(k1) ./ (a(k1) + b(k1));
    v(k1) = m(k1) .* b(k1) ./ ((a(k1) + b(k1)) .* (a(k1) + b(k1) + 1));
end

