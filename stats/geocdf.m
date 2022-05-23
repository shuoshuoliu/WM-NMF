function y = geocdf(x,p,uflag)
%GEOCDF Geometric cumulative distribution function.
%   Y = GEOCDF(X,P) returns the cdf of the geometric distribution with
%   probability parameter P, evaluated at the values in X.  The size of Y
%   is the common size of the input arguments.  A scalar input functions as
%   a constant matrix of the same size as the other input.
%
%   Y = GEOCDF(X,P,'upper') returns the upper tail probability of the 
%   geometric distribution with probability parameter P, evaluated at the 
%   values in X.  
%
%   See also CDF, GEOINV, GEOPDF, GEORND, GEOSTAT.

%   References:
%      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
%          Functions, Dover, New York, 1046pp., sec. 26.1.24.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2017 The MathWorks, Inc. 


if nargin > 2
    uflag = convertStringsToChars(uflag);
end

if nargin < 2
    error(message('stats:geocdf:TooFewInputs'));
end

% Return NaN for out of range parameters.
p(p <= 0 | 1 < p) = NaN;

% Force a zero for negative x.
x(x < 0) = -1;

if ~isfloat(x)
   x = double(x);
end

if nargin>2
    % Compute upper tail
    if ~strcmpi(uflag,'upper')
        error(message('stats:cdf:UpperTailProblem'));
    end
    [errorcode,x,p] = distchck(2,x,p);
    if errorcode > 0
        error(message('stats:geocdf:InputSizeMismatch'));
    end    
    y = zeros(size(x),'like',internal.stats.dominantType(x,p)); % single if x or p is
    k1 = (isnan(x)|isnan(p));
    y(k1) = NaN;
    k2 = (x < 0 & 0 < p & p <= 1);
    y(k2) = 1;
    kc = (~k1 & ~k2);
    if any(kc(:))
        y(kc)= betainc(p(kc),1,floor(x(kc)) + 1,'upper');
    end
else 
    try
        y = 1 - (1 - p) .^ (floor(x) + 1);
    catch
        error(message('stats:geocdf:InputSizeMismatch'));
    end
end
