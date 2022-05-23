function r = geornd(p,varargin)
%GEORND Random arrays from the geometric distribution.
%   R = GEORND(P) returns an array of random numbers chosen from the
%   geometric distribution with probability parameter P.  The size of R is
%   the size of P.
%
%   R = GEORND(P,M,N,...) or R = GEORND(P,[M,N,...]) returns an
%   M-by-N-by-... array.
%
%   See also GEOCDF, GEOINV, GEOPDF, GEOSTAT, RANDOM.

%   GEORND uses the inversion method.

%   References:
%      [1]  Devroye, L. (1986) Non-Uniform Random Variate Generation, 
%           Springer-Verlag.

%   Copyright 1993-2018 The MathWorks, Inc. 


if nargin < 1
    error(message('stats:geornd:TooFewInputs'));
end

[err, sizeOut] = internal.stats.statsizechk(1,p,varargin{:});
if err > 0
    error(message('stats:geornd:InputSizeMismatch'));
end

% Return NaN for elements corresponding to illegal parameter values.
p(p <= 0 | p > 1) = NaN;

% log(1-u) and log(1-p) "really" should both be negative, and
% the abs() here keeps the correct sign for their ratio when
% roundoff makes one of them exactly zero.
r = ceil(abs(log(rand(sizeOut,'like',p)) ./ log(1 - p)) - 1); % == geoinv(u,p)

% Force a zero when p==1, instead of -1.
r(r < 0) = 0;
