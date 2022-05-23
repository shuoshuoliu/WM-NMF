function y = unidpdf(x,n)
%UNIDPDF Uniform (discrete) probability density function (pdf).
%   Y = UNIDPDF(X,N) returns the (discrete) uniform probability 
%   density function on (1,2,...,N) at the values in X.
%
%   The size of Y is the common size of X and N. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   Y is zero (or NaN) unless X is an integer between 1 and N.
%
%   See also UNIDCDF, UNIDINV, UNIDRND, UNIDSTAT, PDF.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin < 2, 
    error(message('stats:unidpdf:TooFewInputs')); 
end

[errorcode, x, n] = distchck(2,x,n);

if errorcode > 0
    error(message('stats:unidpdf:InputSizeMismatch'));
end

% Initialize Y to zero.
y = zeros(size(x),'like',internal.stats.dominantType(x,n)); % Single if x or n is.

k = find(round(x) == x & x >= 1 & x <= n & n ~= 0);
if any(k),
    y(k) = 1 ./ n(k);
end

y(n < 1 | round(n) ~= n) = NaN;
