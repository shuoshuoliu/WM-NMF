function y = tpdf(x,v)
% TPDF  Probability density function (pdf) for Student's T distribution
%   Y = TPDF(X,V) returns the pdf of Student's T distribution with
%   V degrees of freedom, at the values in X.
%   
%   The size of Y is the common size of X and V. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   See also TCDF, TINV, TRND, TSTAT, PDF.

%   References:
%      [1]  E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, New York, 1970, Section 10.3, pages 144-146.

%   Copyright 1993-2014 The MathWorks, Inc. 


if nargin < 2, 
    error(message('stats:tpdf:TooFewInputs')); 
end

[errorcode,x,v] = distchck(2,x,v);

if errorcode > 0
    error(message('stats:tpdf:InputSizeMismatch'));
end

% Initialize Y to zero, or NaN for invalid d.f.
y = zeros(size(x),'like',internal.stats.dominantType(x,v));
y(v <= 0) = NaN;

% Use gammaln function to avoid overflows.
k = find(v>0 & v<Inf);
if any(k)
    vk = v(k);
    term = exp(gammaln((vk + 1) / 2) - gammaln(vk/2));
    y(k) = term ./ (sqrt(vk*pi) .* (1 + (x(k) .^ 2) ./ vk) .^ ((vk + 1)/2));
end

k = find(v==Inf); % approaches normal as v->Inf
if any(k)
    y(k) = normpdf(x,0,1);
end
