function y = gampdf(x,a,b)
%GAMPDF Gamma probability density function.
%   Y = GAMPDF(X,A,B) returns the gamma probability density function with
%   shape and scale parameters A and B, respectively, at the values in X.
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Some references refer to the gamma distribution with a single
%   parameter.  This corresponds to the default of B = 1.
%
%   See also GAMCDF, GAMFIT, GAMINV, GAMLIKE, GAMRND, GAMSTAT, GAMMA,
%            GAMMALN.

%   References:
%      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
%          Functions, Dover, New York, section 26.1.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley.

%   Copyright 1993-2019 The MathWorks, Inc.


if nargin < 2
    error(message('stats:gampdf:TooFewInputs'));
elseif nargin < 3
    b = 1;
end

[errorcode, x, a, b] = distchck(3,x,a,b);

if errorcode > 0
    error(message('stats:gampdf:InputSizeMismatch'));
end

% Initialize y to zero.
y = zeros(size(x),'like',internal.stats.dominantType(x,a,b));

% Return NaN for out of range parameters.
y(a < 0) = NaN;
y(b <= 0) = NaN;
y(isnan(a) | isnan(b) | isnan(x)) = NaN;

% Scale
z = x./b;

% Special cases
i = z==0 & a==1 & b>0;
y(i) = 1./b(i);
y(z==0 & a<1 & b>0) = Inf;

% Normal cases
k = find(z>0 & z<Inf & a>0 & a<Inf & b>0);

if ~isempty(k)
    z = z(k);
    a = a(k)-1;
    b = b(k);
    
    i = a<0;
    noti = ~i;
    y(k(i)) = f(z(i),a(i)+1) .* exp(log(a(i)+1)-log(z(i))) ./ b(i);
    y(k(noti)) = f(z(noti),a(noti)) ./ b(noti);
end

end


function y = f(z,a)
% Compute gampdf without error checking for z>0 and a>0.
y = zeros(size(z),'like',internal.stats.dominantType(z,a));

% z term dominates
i1 = a<=realmin*z;
y(i1) = exp(-z(i1));

% Normal expansion through logs
i2 = z<realmin*a;
y(i2) = exp( a(i2).*log(z(i2)) -z(i2) -gammaln(a(i2)+1) );

% Loader's saddle point expansion
i3 = ~i1 & ~i2;
lnsr2pi = 0.9189385332046727; % log(sqrt(2*pi))
ai3 = a(i3);
y(i3) = exp(-lnsr2pi -0.5*log(ai3) - stirlerr(ai3) ...
    - binodeviance(ai3,z(i3)));
end
