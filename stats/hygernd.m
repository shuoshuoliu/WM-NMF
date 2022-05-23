function r = hygernd(m,k,n,varargin)
%HYGERND Random arrays from the hypergeometric distribution.
%   R = HYGERND(M,K,N) returns an array of random numbers chosen from a
%   hypergeometric distribution with parameters M, K, and N.  The size of R
%   is the common size of M, K, and N if all are arrays.  If any parameters
%   are scalars, the size of R is the size of the other parameter(s).
%
%   R = HYGERND(M,K,N,MM,NN,...) or R = HYGERND(M,K,N,[MM,NN,...]) returns
%   an MM-by-NN-by-... array.
%
%   See also HYGECDF, HYGEINV, HYGEPDF, HYGESTAT, RANDOM.

%   HYGERND uses the inversion method.

%   Copyright 1993-2019 The MathWorks, Inc. 


if nargin < 3
    error(message('stats:hygernd:TooFewInputs')); 
end

[err, sizeOut, numelOut] = internal.stats.statsizechk(3,m,k,n,varargin{:});
if err > 0
    error(message('stats:hygernd:InputSizeMismatch'));
end
tyOut = internal.stats.dominantType(m,k,n);

% Return NaN for elements corresponding to illegal parameter values.
k1 = find((n <= m & k <= m) & (0 <= n & round(n) == n) ...
        & (0 <= k & round(k) == k) & (0 <= m & round(m) == m));

% Handle the simplest situation, all "nice" parameter values, efficiently
if (isscalar(m) && isscalar(k) && isscalar(n) && numel(k1) == 1) ...   & scalar params
                                              || numel(k1) == numelOut % at least one matrix params
    % Generate uniform random values, and apply the hypergeometric inverse CDF.
    r = hygeinv(rand(sizeOut,"like",tyOut), m, k, n);

% Handle calls with some illegal parameter values
else
    r = NaN(sizeOut, "like", tyOut);
    if numel(m) > 1, m = m(k1); end
    if numel(k) > 1, k = k(k1); end
    if numel(n) > 1, n = n(k1); end
    r(k1) = hygeinv(rand(size(k1),"like",tyOut), m, k, n);
end
