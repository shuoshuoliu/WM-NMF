function x = unidinv(p,n)
%UNIDINV Inverse of uniform (discrete) distribution function.
%   X = UNIDINV(P,N) returns the inverse of the uniform
%   (discrete) distribution function at the values in P.
%   X takes the values (1,2,...,N).
%   
%   The size of X is the common size of P and N. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   See also UNIDCDF, UNIDPDF, UNIDRND, UNIDSTAT, ICDF.

%   Copyright 1993-2015 The MathWorks, Inc.


if nargin < 2, 
    error(message('stats:unidinv:TooFewInputs')); 
end

[errorcode,p,n] = distchck(2,p,n);

if errorcode > 0
    error(message('stats:unidinv:InputSizeMismatch'));
end

% Initialize X to zero.
x = zeros(size(p),'like',internal.stats.dominantType(p,n)); % Single if p or n is.

k = p > 0 & p <= 1 & n >= 1 & round(n) == n;

% Return NaN if the arguments are outside their respective limits.
x(~k) = NaN;

if any(k(:))
    x(k) = ceil(n(k) .* p(k));
end
