function r = betarnd(a,b,varargin);
%BETARND Random arrays from beta distribution.
%   R = BETARND(A,B) returns an array of random numbers chosen from the
%   beta distribution with parameters A and B.  The size of R is the common
%   size of A and B if both are arrays.  If either parameter is a scalar,
%   the size of R is the size of the other parameter.
%
%   R = BETARND(A,B,M,N,...) or R = BETARND(A,B,[M,N,...]) returns an
%   M-by-N-by-... array.
%
%   See also BETACDF, BETAINV, BETALIKE, BETAPDF, BETASTAT, RANDOM.

%   BETARND uses a transformation method, expressing a beta random variable
%   in terms of gamma random variables (Devroye, page 430).

%   Reference:
%      [1]  Devroye, L. (1986) Non-Uniform Random Variate Generation, 
%           Springer-Verlag.

%   Copyright 1993-2018 The MathWorks, Inc.


if nargin < 2
    error(message('stats:betarnd:TooFewInputs'));
end

[err, sizeOut] = internal.stats.statsizechk(2,a,b,varargin{:});
if err > 0
    error(message('stats:betarnd:InputSizeMismatch'));
end

% Generate gamma random values and take ratio of the first to the sum.
g1 = randg(a,sizeOut); % could be Infs or NaNs
g2 = randg(b,sizeOut); % could be Infs or NaNs
r = g1 ./ (g1 + g2);

% For a and b both very small, we often get 0/0.  Since the distribution is
% essentially a Bernoulli(a/(a+b)), we can replace those NaNs.
t = (g1==0 & g2==0);
if any(t(:))
    p = a ./ (a+b);
    if ~isscalar(p), p = p(t); end
    r(t) = binornd(1,p(:),sum(t(:)),1);
end
