function p = poisscdf(x,lambda,uflag)
%POISSCDF Poisson cumulative distribution function.
%   P = POISSCDF(X,LAMBDA) computes the Poisson cumulative
%   distribution function with parameter LAMBDA at the values in X.
%
%   The size of P is the common size of X and LAMBDA. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   P = POISSCDF(X,LAMBDA,'upper') computes the upper tail probability of
%   the Poisson distribution with parameter LAMBDA at the values in X.
%
%   See also POISSFIT, POISSINV, POISSPDF, POISSRND, POISSTAT.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.22.

%   Copyright 1993-2017 The MathWorks, Inc.

 
if nargin > 2
    uflag = convertStringsToChars(uflag);
end

if nargin < 2, 
    error(message('stats:poisscdf:TooFewInputs')); 
end

[errorcode, x, lambda] = distchck(2,x,lambda);

if errorcode > 0
    error(message('stats:poisscdf:InputSizeMismatch'));
end

% Initialize P to zero.
p = zeros(size(x),'like',internal.stats.dominantType(x,lambda)); % Single if x or lambda is
p(isnan(x) | isnan(lambda)) = NaN;
if ~isfloat(x)
   x = double(x);
end
x = floor(x);

% Return NaN for invalid or indeterminate inputs
t = (lambda<0) | (x==Inf & lambda==Inf);
if any(t(:))
    p(t) = NaN;
end
todo = (x>=0) & ~t & isfinite(lambda);

% Return 1 for x=Inf as long as lambda is valid
t = (x==Inf & lambda>0 & isfinite(lambda));
if any(t(:))
    todo(t) = false;
    if nargin>2 && strcmpi(uflag,'upper')
        p(t) = 0;
    else
        p(t) = 1;
    end
end

t = (x<0 & lambda>0 & isfinite(lambda));
if any(t(:))
    todo(t) = false;
    if nargin>2 && strcmpi(uflag,'upper')
        p(t) = 1;
    else
        p(t) = 0;
    end
end

% Compute P when X for the remaining cases
x = x(todo);
lambda = lambda(todo);

if nargin>2
    if ~strcmpi(uflag,'upper')
        error(message('stats:cdf:UpperTailProblem'));
    end
    p(todo) = gammainc(lambda,floor(x+1));
else
    p(todo) = gammainc(lambda,x+1,'upper');
end
