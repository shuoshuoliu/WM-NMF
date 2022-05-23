function [mn,v] = hygestat(m,k,n);
%HYGESTAT Mean and variance for the hypergeometric distribution.
%   [MN,V] = HYGESTAT(M,K,N) returns the mean and variance 
%   of the hypergeometric distribution with parameters M, K, and N.
%
%   See also HYGECDF, HYGEINV, HYGEPDF, HYGERND.

%   Reference:
%      [1]  Mood, Alexander M., Graybill, Franklin A. and Boes, Duane C.,
%      "Introduction to the Theory of Statistics, Third Edition", McGraw Hill
%      1974 p. 538.

%   Copyright 1993-2019 The MathWorks, Inc. 


if nargin < 3
   error(message('stats:hygestat:TooFewInputs')); 
end

[errorcode, m, k, n] = distchck(3,m,k,n);

if errorcode > 0
   error(message('stats:hygestat:InputSizeMismatch'));
end

% Initialize the mean and variance to zero.
outType = internal.stats.dominantType(m,k,n);
mn = zeros(size(m),"like",outType);
v = mn;

%   Return NaN for values of the parameters outside their respective limits.
k1 = (m < 0 | k < 0 | n < 0 | round(m) ~= m | round(k) ~= k | ...
     round(n) ~= n | n > m | k > m);
if any(k1,"all")
   mn(k1) = NaN;
   v(k1)  = NaN;
end

kc = 1:numel(m);
kc(k1) = [];

if any(kc)
   nc = n(kc);
   mc = m(kc);
   mn(kc) = nc .* k(kc) ./ mc;
   v(kc) = nc .* k(kc) .* (mc - k(kc)) .* (mc - nc) ./ (mc .* mc .* (mc - 1));
end
