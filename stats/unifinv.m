function x = unifinv(p,a,b)
%UNIFINV Inverse of uniform (continuous) distribution function.
%   X = UNIFINV(P,A,B) returns the inverse of the uniform
%   (continuous) distribution function on the interval [A,B],
%   at the values in P. By default A = 0 and B = 1.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also UNIFCDF, UNIFIT, UNIFPDF, UNIFRND, UNIFSTAT, ICDF.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin < 1, 
    error(message('stats:unifinv:TooFewInputs')); 
end

if nargin == 1
    a = 0;
    b = 1;
end

[errorcode, p, a, b] = distchck(3,p,a,b);

if errorcode > 0
    error(message('stats:unifinv:InputSizeMismatch'));
end

% Initialize X to zero.
x = zeros(size(p),'like',internal.stats.dominantType(p,a,b)); % single if p, a, or b is single

% Return NaN if the arguments are outside their respective limits.
x(a >= b | p < 0 | p > 1) = NaN;

k = find(~(a >= b | p < 0 | p > 1));
if any(k)
    x(k) = a(k) + p(k) .* (b(k) - a(k));
end
