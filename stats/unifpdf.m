function y = unifpdf(x,a,b)
%UNIFPDF Uniform (continuous) probability density function (pdf).
%   Y = UNIFPDF(X,A,B) returns the continuous uniform pdf on the
%   interval [A,B] at the values in X. By default A = 0 and B = 1.
%   
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also UNIFCDF, UNIFINV, UNIFIT, UNIFRND, UNIFSTAT, PDF.

%   Reference:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.34.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin < 1
    error(message('stats:unifpdf:TooFewInputs')); 
end

if nargin == 1
    a = 0;
    b = 1;
end

[errorcode, x, a, b] = distchck(3,x,a,b);

if errorcode > 0
    error(message('stats:unifpdf:InputSizeMismatch'));
end

% Initialize Y to zero.
y = zeros(size(x),'like',internal.stats.dominantType(x,a,b)); % single if any of x, a, b is single

y(a >= b) = NaN;

k = find(x >= a & x <= b & a < b);
if any(k),
    y(k) = 1 ./ (b(k) - a(k));
end

% Pass NaN inputs through to outputs
y(isnan(x)) = NaN;
