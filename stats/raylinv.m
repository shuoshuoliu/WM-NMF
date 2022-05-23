function x = raylinv(p,b)
%RAYLINV  Inverse of the Rayleigh cumulative distribution function (cdf).
%   X = RAYLINV(P,B) returns the Rayleigh cumulative distribution 
%   function with parameter B at the probabilities in P.
%
%   The size of X is the common size of P and B. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   See also RAYLCDF, RAYLFIT, RAYLPDF, RAYLRND, RAYLSTAT, ICDF.

%   Reference:
%      [1]  Evans, Merran, Hastings, Nicholas and Peacock, Brian,
%      "Statistical Distributions, Second Edition", Wiley
%      1993 p. 134-136.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin < 1
    error(message('stats:raylinv:TooFewInputs')); 
end
if nargin<2
    b = 1;
end

[errorcode, p, b] = distchck(2,p,b);

if errorcode > 0
    error(message('stats:raylinv:InputSizeMismatch'));
end

% Initialize x to zero.
x = zeros(size(p),'like',internal.stats.dominantType(p,b)); % single if p or b is


% Return NaN if the arguments are outside their respective limits.
x(b <= 0| p < 0 | p > 1) = NaN;

% Put in the correct values when P is 1.
x(p == 1) = Inf;

k=find(b > 0 & p > 0  &  p < 1);
if any(k),
    pk = p(k);
    bk = b(k);
    x(k) = sqrt((-2*bk .^ 2) .* log(1 - pk));
end

% Return NaN if the input probability values are NaN.
x(isnan(p)) = NaN;

