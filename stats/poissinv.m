function x = poissinv(p,lambda)
%POISSINV Inverse of the Poisson cumulative distribution function (cdf).
%   X = POISSINV(P,LAMBDA) returns the inverse of the Poisson cdf 
%   with parameter lambda. Since the Poisson distribution is discrete,
%   POISSINV returns the smallest value of X, such that the Poisson 
%   cdf evaluated, at X, equals or exceeds P.
%
%   The size of X is the common size of P and LAMBDA. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   See also POISSCDF, POISSFIT, POISSPDF, POISSRND, POISSTAT, ICDF.

%   Copyright 1993-2019 The MathWorks, Inc.


if nargin < 2
    error(message('stats:poissinv:TooFewInputs')); 
end

[errorcode, p, lambda] = distchck(2,p,lambda);

if errorcode > 0
    error(message('stats:poissinv:InputSizeMismatch'));
end

x = zeros(size(p),'like',internal.stats.dominantType(p,lambda)); % Single if p or lambda is
x(isnan(p) | isnan(lambda)) = NaN;

% Return NaN if the arguments are outside their respective limits.
x(lambda < 0 | p < 0 | p > 1) = NaN;

% Return Inf if p = 1 and lambda is positive.
x(lambda > 0 & p == 1) = Inf;

% For large lambda, use normal approximation.
k = lambda>10 & p>0 & p<1;
if any(k(:))
    lambda1 = lambda(k);
    p1 = p(k);
    x1 = -sqrt(2)*erfcinv(2*p1);
    count = max(ceil(sqrt(lambda1).*x1 + lambda1),0);
    y = gammainc(lambda1,count+1,'upper');
    
    % Check upward
    ynew = zeros(size(y));
    xnew = count;
    under = y<p1;
    while any(under(:))
        ynew(under) = gammainc(lambda1(under),xnew(under)+2,'upper');
        xnew(under) = xnew(under)+1;
        under = under & ynew<p1;
    end
    count = xnew;
    
    % Check downward
    ynew = zeros(size(y),'like',y);
    xnew = count;
    over = y>p1 & ~under;
    while any(over(:))
        ynew(over) = gammainc(lambda1(over),xnew(over),'upper');
        over = over & ynew>=p1;
        xnew(over) = xnew(over)-1;
    end
    count = xnew;

    % Get the final count
    x(k) = count;
end

% For small lambda, start at x=0.
k = lambda > 0 & lambda <= 10 & p > 0 & p < 1;
if any(k(:))
    p1 = p(k);
    lambda1 = lambda(k);
    count = zeros(size(p1),'like',x); % will be assigned into X so must match X
    cumdist = zeros(size(p1),'like',x);
    ok = true(size(p1));
    while any(ok(:))
        cumdist(ok) = gammainc(lambda1(ok),count(ok)+1,'upper');
        ok = ok & cumdist<p1;
        count(ok) = count(ok) + 1;
    end
    x(k) = count;
end
