function [m,v]= unidstat(n)
%UNIDSTAT Mean and variance for uniform (discrete) distribution.
%   [M,V] = UNIDSTAT(N) returns the mean and variance of
%   the (discrete) uniform distribution on {1,2,...,N}
%
%   See also UNIDCDF, UNIDINV, UNIDPDF, UNIDRND.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin < 1, 
    error(message('stats:unidstat:TooFewInputs'));   
end

%   Initialize the mean and variance to zero.
m = zeros(size(n),'like',n); % Single if n is.
v = m;

k = find(n > 0 & round(n) == n);
m(k) = (n(k) + 1) / 2;
v(k) = (n(k) .^ 2 - 1) / 12;

k1 = (n <= 0 | round(n) ~= n);
if any(k1(:))
    m(k1) = NaN;
    v(k1) = NaN;
end
