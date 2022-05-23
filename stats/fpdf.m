function y = fpdf(x,v1,v2)
%FPDF   F probability density function.
%   Y = FPDF(X,V1,V2) returns the F distribution probability density
%   function with V1 and V2 degrees of freedom at the values in X.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also FCDF, FINV, FRND, FSTAT, PDF.

%   References:
%      [1] J. K. Patel, C. H. Kapadia, and D. B. Owen, "Handbook
%      of Statistical Distributions", Marcel-Dekker, 1976.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin < 3, 
    error(message('stats:fpdf:TooFewInputs')); 
end

[errorcode, x, v1, v2] = distchck(3,x,v1,v2);

if errorcode > 0
    error(message('stats:fpdf:InputSizeMismatch'));
end

% Initialize y to zero; it stays zero for x<0
y = zeros(size(x),'like',internal.stats.dominantType(x,v1,v2)); % Single if x, v1, or v2 is.

% Special case for invalid parameters
k1 = ( (v1 <= 0) | (v2 <= 0) | isnan(x) | isnan(v1) | isnan(v2) );
if any(k1(:))
    y(k1) = NaN;
end

% Regular case with valid parameters and positive x
k = ( (x > 0 & x < Inf) & (v1 > 0 & v1 < Inf) & (v2 > 0 & v2 < Inf) );
if any(k(:))
    xk = x(k);    
    ratio = v1(k) .* xk ./ v2(k);    
    logtemp = -log(xk) - 0.5 * v1(k) .* log1p(1./ratio) - 0.5 * v2(k) .* log1p(ratio) - betaln(v1(k)/2, v2(k)/2);    
    y(k) = exp(logtemp);
end

% Special case - if either v1 or v2 is too large then return NaN
const = 1e8;
k = ( (v1>const) & (v1<Inf) ) | ( (v2>const) & (v2<Inf) );
k = k & (x > 0 & x < Inf);
if any(k(:))
    y(k) = NaN;
end

% Special case - at x==Inf, result should be 0, regardless of any valid v1 or v2
k = (x==Inf) & (v1>0) & (v2>0);
if any(k(:))
   y(k) = 0; 
end

% Special case - at x==0, result depends only on the value of v1
k = (x==0) & (v1>0) & (v2>0);
if any(k(:))
    y(k & v1==2) = 1;
    y(k & v1<2) = Inf;
    y(k & v1>2) = 0;
end

% Special cases - v1 and/or v2 = Inf

% If v1 = Inf AND v2 is not Inf AND both v1, v2 > 0 AND x is not NaN
% then fpdf( t, v1, v2 ) = (v2/t^2) * chi2pdf( v2/t, v2 )
kinf1 = (v2 > 0 & v2 < Inf) & (v1 == Inf) & (x > 0 & x < Inf);
if any(kinf1(:))
    y(kinf1) = ( v2(kinf1) ./ (x(kinf1).^2) ) .* chi2pdf( v2(kinf1) ./ x(kinf1) , v2(kinf1) );
end
clear kinf1;

% If v1 is not Inf AND v2 = Inf AND both v1, v2 > 0 AND x is not NaN
% then fpdf( t, v1, v2 ) = v1 * chi2pdf( v1*t, v1 )
kinf2 = (v1 > 0 & v1 < Inf) & (v2 == Inf) & (x > 0 & x < Inf);
if any(kinf2(:))
    y(kinf2) = v1(kinf2) .* chi2pdf( v1(kinf2) .* x(kinf2), v1(kinf2) );
end
clear kinf2;

% If v1 is Inf AND v2 = Inf AND x is not NaN
% then set fpdf( t, v1, v2 ) = 0 if t not equal to 1 and Inf if t = 1
kinf3 = (v1 == Inf) & (v2 == Inf) & (x > 0 & x < Inf);
if any(kinf3(:))
    % set all kinf3 values to 0
    y(kinf3) = 0;
    % if x = 1 then set values to Inf
    y(kinf3 & (x==1)) = Inf;
end
clear kinf3;

end
