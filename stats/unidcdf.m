function p = unidcdf(x,n,uflag)
%UNIDCDF Uniform (discrete) cumulative distribution function.
%   P = UNIDCDF(X,N) returns the cumulative distribution function
%   for a random variable uniform on (1,2,...,N), at the values in X.
%
%   The size of P is the common size of X and N. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   P = UNIDCDF(X,N,'upper') returns the upper tail probability
%   for a random variable uniform on (1,2,...,N), at the values in X.
%
%   See also UNIDINV, UNIDPDF, UNIDRND, UNIDSTAT, CDF.

%   Copyright 1993-2017 The MathWorks, Inc. 


if nargin > 2
    uflag = convertStringsToChars(uflag);
end

if nargin < 2, 
    error(message('stats:unidcdf:TooFewInputs')); 
end

[errorcode, x, n] = distchck(2,x,n);

if errorcode > 0
    error(message('stats:unidcdf:InputSizeMismatch'));
end

% Initialize P to zero.
p = zeros(size(x),'like',internal.stats.dominantType(x,n)); % Single if x or n is.

% P = 1 when X >= N
p(x >= n) = 1;

xx=floor(x);

k = find(xx >= 1 & xx <= n);
if any(k),
    p(k) = xx(k) ./ n(k);
end

p(n < 1 | round(n) ~= n) = NaN;

if nargin>2
    % Compute upper tail
    if ~strcmpi(uflag,'upper')
        error(message('stats:cdf:UpperTailProblem'));
    end
    p = 1 - unidcdf(x,n);
end