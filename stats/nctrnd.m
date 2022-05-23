function r = nctrnd(v,delta,varargin)
%NCTRND Random arrays from the noncentral t distribution.
%   R = NCTRND(V,DELTA) returns an array of random numbers chosen from the
%   noncentral t distribution with parameters V and DELTA.  The size of R is
%   the common size of V and DELTA if both are arrays.  If either parameter
%   is a scalar, the size of R is the size of the other parameter.
%
%   R = NCTRND(V,DELTA,M,N,...) or R = NCTRND(V,DELTA,[M,N,...]) returns an
%   M-by-N-by-... array.
%
%   See also NCTCDF, NCTINV, NCTPDF, NCTSTAT, TRND, RANDOM.

%   NCTRND generates values using the definition of a noncentral t random
%   variable, as the ratio of a normal with non-zero mean and the sqrt of a
%   chi-square.

%   Reference:
%      [1] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2018 The MathWorks, Inc. 


if nargin < 2
    error(message('stats:nctrnd:TooFewInputs'));
end

[err, sizeOut] = internal.stats.statsizechk(2,v,delta,varargin{:});
if err > 0
    error(message('stats:nctrnd:InputSizeMismatch'));
end

% Return NaN for elements corresponding to illegal parameter values.
% Non-integer degrees of freedom are allowed.
v(v <= 0) = NaN;

% Infinite degrees of freedom is equivalent to a normal, and valid.
% Prevent Inf/Inf==NaN for the standardized chi-square in the denom.
v(isinf(v)) = realmax;

outType = real(full(internal.stats.dominantType(v,delta)));
r = (randn(sizeOut, "like", outType)+delta) ...
    ./ sqrt(2.*randg(v./2,sizeOut, "like", outType) ./ v);
