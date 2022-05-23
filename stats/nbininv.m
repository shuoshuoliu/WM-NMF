function x = nbininv(y,r,p)
%NBININV Inverse of the negative binomial cumulative distribution function (cdf).
%   X = NBININV(Y,R,P) returns the inverse of the negative binomial cdf with
%   parameters R and P. Since the negative binomial distribution is discrete,
%   NBININV returns the least integer X such that the negative
%   binomial cdf evaluated at X, equals or exceeds Y.
%
%   The size of X is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also NBINCDF, NBINFIT, NBINLIKE, NBINPDF, NBINRND, NBINSTAT, ICDF.

%   Copyright 1993-2017 The MathWorks, Inc.


if nargin < 3,
    error(message('stats:nbininv:TooFewInputs'));
end

[errorcode y r p] = distchck(3,y,r,p);

if errorcode > 0
    error(message('stats:nbininv:InputSizeMismatch'));
end

% Initialize X to 0.
if isa(y,'single') || isa(r,'single') || isa(p,'single')
    x = zeros(size(y),'single');
else
    x = zeros(size(y));
end

% Out of range or missing parameters or probabilities return NaN.
% Infinite values for R correspond to a Poisson, but its mean cannot
% be determined from the (R,P) parametrization.
nans = ~(0 < r & isfinite(r) & 0 < p & p <= 1) | ~(0 <= y & y <= 1);
x(nans) = NaN;

% Compute X when 0 <= Y <= 1.
k = find(~nans);

% Return Inf if the probability is 1.
k1 = find(y(k) == 1);
if any(k1)
    x(k(k1)) = Inf;
    k(k1) = [];
end

% Accumulate probabilities to satisfy the values in Y.
if any(k)
    cumdist = zeros(size(y));
    fudge = 40*eps(y);
    count = 0;
    cumdist(k) = nbinpdf(0,r(k),p(k));
    k = k(cumdist(k) < y(k));
    while ~isempty(k)
        x(k) = x(k) + 1;
        count = count + 1;
        % We are essentially doing this:
        % cumdist(k) = cumdist(k) + nbinpdf(count,r(k),p(k));
        cumdist(k) = iUpdateCDF(cumdist(k),count,r(k),p(k));
        k = k(cumdist(k) + fudge(k) < y(k));
    end
end
end

function cumdist = iUpdateCDF(cumdist,count,r,p)
t = cumdist > 0.5;
if any(t)
    cumdist(t) = 1 - nbincdf(count,r(t),p(t),'upper');
end
if any(~t)
    cumdist(~t) = cumdist(~t) + nbinpdf(count,r(~t),p(~t));
end
end