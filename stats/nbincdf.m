function y = nbincdf(x,r,p,uflag)
%NBINCDF Negative binomial cumulative distribution function.
%   Y=NBINCDF(X,R,P) returns the negative binomial cumulative distribution
%   function with parameters R and P at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   The algorithm uses the cumulative sums of the binomial masses.
%
%   Y=NBINCDF(X,R,P,'upper') returns the upper tail probability of 
%   the negative binomial distribution with parameters R and P 
%   at the values in X.
%
%   See also NBINFIT, NBININV, NBINLIKE, NBINPDF, NBINRND, NBINSTAT, CDF.

%   Copyright 1993-2013 The MathWorks, Inc. 


if nargin > 3
    uflag = convertStringsToChars(uflag);
end

if nargin < 3, 
    error(message('stats:nbincdf:TooFewInputs')); 
end 

scalarrp = (isscalar(r) & isscalar(p));

[errorcode, x, r, p] = distchck(3,x,r,p);

if errorcode > 0
    error(message('stats:nbincdf:InputSizeMismatch'));
end

% Initialize Y to 0.
if isa(x,'single') || isa(r,'single') || isa(p,'single')
   y = zeros(size(x),'single');
else
   y = zeros(size(x));
end
if ~isfloat(x)
   x = double(x);
end

% Out of range or missing parameters and missing data return NaN.
% Infinite values for R correspond to a Poisson, but its mean cannot
% be determined from the (R,P) parametrization.
nans = ~(0 < r & isfinite(r) & 0 < p & p <= 1) | isnan(x);
y(nans) = NaN;

% Compute Y when 0 <= X.  Data outside this range return 0.
xx = floor(x);
k = find(0 <= xx  &  ~nans);

% Positive infinite values of X return 1.
k1 = find(isinf(xx(k)));
if any(k1) 
    if nargin>3 && strcmpi(uflag,'upper')
        y(k(k1)) = 0;
    else
        y(k(k1)) = 1;
    end
    k(k1) = [];
end

% Uppertail when X<=0, return 1.
k1 = (x<0 & ~nans);
if any(k1)
    if nargin>3 && strcmpi(uflag,'upper')
        y(k1) = 1;
    end
end

% Accumulate probabilities up to the maximum value in X.
if any(k)
    if nargin>3
        if ~strcmpi(uflag,'upper')
            error(message('stats:cdf:UpperTailProblem'));
        end
        y(k)=betainc(p(k),r(k),xx(k)+1,'upper');
    else
        val = max(xx(k));
        if scalarrp
            tmp = cumsum(nbinpdf(0:val,r(1),p(1)));
            y(k) = tmp(xx(k) + 1);
        else
            idx = (0:val)';
            compare = idx(:,ones(size(k)));
            index = xx(k);
            index = index(:);
            index = index(:,ones(size(idx)))';
            rbig = r(k);
            rbig = rbig(:);
            rbig = rbig(:,ones(size(idx)))';
            pbig = p(k);
            pbig = pbig(:);
            pbig = pbig(:,ones(size(idx)))';
            y0 = nbinpdf(compare,rbig,pbig);
            indicator = find(compare > index);
            y0(indicator) = zeros(size(indicator));
            y(k) = sum(y0,1);
        end
    end
end

% Make sure that round-off errors never make Y greater than 1.
y(y > 1) = 1;
