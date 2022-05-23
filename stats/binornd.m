function r = binornd(n,p,varargin)
% BINORND Random arrays from the binomial distribution.
%   R = BINORND(N,P,MM,NN) returns an array of random numbers chosen from a
%   binomial distribution with parameters N and P.  The size of R is the
%   common size of N and P if both are arrays.  If either parameter is a
%   scalar, the size of R is the size of the other parameter.
%   
%   R = BINORND(N,P,MM,NN,...) or R = BINORND(N,P,[MM,NN,...]) returns an
%   MM-by-NN-by-... array.
%
%   See also BINOCDF, BINOINV, BINOPDF, BINOSTAT, RANDOM.

%   BINORND generates values using the definition of the binomial
%   distribution, as a sum of Bernoulli random variables.  See Devroye,
%   Lemma 4.1 on page 428, method on page 524.

%   References:
%      [1]  Devroye, L. (1986) Non-Uniform Random Variate Generation, 
%           Springer-Verlag.

%   Copyright 1993-2018 The MathWorks, Inc.


if nargin < 2
    error(message('stats:binornd:TooFewInputs')); 
end

[err, sizeOut] = internal.stats.statsizechk(2,n,p,varargin{:});
ndimsOut = numel(sizeOut);
if err > 0
    error(message('stats:binornd:InputSizeMismatch'));
end

outType = internal.stats.dominantType(n,p); % Single if n or p is.

% Handle the scalar params case efficiently
if isscalar(n) && isscalar(p)
    if (0 <= p && p <= 1) && (0 <= n && round(n) == n)
        r = cast(...
            sum(rand([sizeOut,n],'like',outType) < p, ndimsOut+1), ...
            'like', outType);
    else
        r = NaN(sizeOut,'like',outType);
    end

% Handle the scalar n case efficiently
elseif isscalar(n)
    if 0 <= n && round(n) == n
        r = cast(...
            sum(bsxfun(@lt, rand([sizeOut,n],'like',outType), p), ndimsOut+1), ...
            'like', outType);
        r(p < 0 | 1 < p) = NaN;
    else
        r = NaN(sizeOut,'like',outType);
    end

% Handle the non-scalar params case
else 
    if isscalar(p), p = p * ones(sizeOut,'like',outType); end
    r = zeros(sizeOut,'like',outType);
    if ~isempty(n)
        for i = 1:max(n(:))
            k = find(n >= i);
            r(k) = r(k) + (rand(size(k),'like',outType) < p(k));
        end
        r(p < 0 | 1 < p | n < 0 | round(n) ~= n) = NaN;
    end
end
