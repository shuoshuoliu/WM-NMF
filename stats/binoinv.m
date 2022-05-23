function x = binoinv(y,n,p)
%BINOINV Inverse of the binomial cumulative distribution function (cdf).
%   X = BINOINV(Y,N,P) returns the inverse of the binomial cdf with 
%   parameters N and P. Since the binomial distribution is
%   discrete, BINOINV returns the least integer X such that 
%   the binomial cdf evaluated at X, equals or exceeds Y.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   Note that X takes the values 0,1,2,...,N.
%
%   See also BINOCDF, BINOFIT, BINOPDF, BINORND, BINOSTAT, ICDF.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin < 3, 
    error(message('stats:binoinv:TooFewInputs')); 
end 

[errorcode, y, n, p] = distchck(3,y,n,p);

if errorcode > 0
    error(message('stats:binoinv:InputSizeMismatch'));
end

k = 1:numel(y);
k1 = find(n < 0 | p < 0 | p > 1 | round(n) ~= n | y < 0 | y > 1 | isnan(y)); 
k2 = find(y == 1);

% Initialize X to 0.
x = zeros(size(y),'like',internal.stats.dominantType(y,n,p)); % Single if y, n or p is.

cumdist = x;
x(k2) = n(k2);
x(k1) = NaN;
k([k1(:); k2(:)]) = [];
if isempty(k), return; end

cumdist(k) = binopdf(0,n(k),p(k));

count = 0;

k = k(cumdist(k) < y(k));
while ~isempty(k)
   x(k) = x(k) + 1;
   count = count + 1;
   cumdist(k) = cumdist(k) + binopdf(count,n(k),p(k));
   k(cumdist(k) >= y(k)) = [];
   k(x(k) >= n(k)) = [];
end

