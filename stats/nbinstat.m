function [m, v] = nbinstat(r,p)
%NBINSTAT Mean and variance of the negative binomial distribution.
%   [M, V] = NBINSTAT(R,P) returns the mean and variance of the
%   negative binomial distribution with parameters R and P.
%
%   See also NBINCDF, NBINFIT, NBININV, NBINLIKE, NBINPDF, NBINRND.

%   Copyright 1993-2019 The MathWorks, Inc. 


if nargin < 2
    error(message('stats:nbinstat:TooFewInputs')); 
end

[errorcode, r, p] = distchck(2,r,p);

if errorcode > 0
    error(message('stats:nbinstat:InputSizeMismatch'));
end

if isa(r,'single') || isa(p,'single')
   m = zeros(size(r),'single');
else
   m = zeros(size(r));
end
v = m;

% Out of range or missing parameters return NaN.  Infinite values
% for R correspond to a Poisson, but its stats cannot be determined
% from the (R,P) parametrization.
nans = ~(0 < r & isfinite(r) & 0 < p & p <= 1);
m(nans) = NaN;
v(nans) = NaN;

k = find(~nans);
if any(k)
    q = 1 - p(k);
    m(k) = r(k) .* q ./ p(k);
    v(k) = r(k) .* q ./ (p(k) .* p(k));
end
