function [m,v]= fstat(v1,v2)
%FSTAT  Mean and variance for the F distribution.
%   [M,V] = FSTAT(V1,V2) returns the mean (M) and variance (V)
%   of the F distribution with V1 and V2 degrees of freedom.
%   Note that the mean of the F distribution is undefined if V1 
%   is less than 3. The variance is undefined for V2 less than 5.
%
%   See also FCDF, FINV, FPDF, FRND.

%   References:
%      [1]  W. H. Beyer, "CRC Standard Probability and Statistics",
%      CRC Press, Boston, 1991, page 23.

%   Copyright 1993-2015 The MathWorks, Inc. 


if nargin < 2, 
    error(message('stats:fstat:TooFewInputs')); 
end

[errorcode, v1, v2] = distchck(2,v1,v2);

if errorcode > 0
    error(message('stats:fstat:InputSizeMismatch'));
end

%   Initialize the M and V to zero.
m = zeros(size(v1),'like',internal.stats.dominantType(v1,v2)); % Single if v1, or v2 is.
v = m;

k = find(v1 <= 0 | v2 <= 0 | isnan(v1) | isnan(v2));
if any(k(:))
    m(k) = NaN;
    v(k) = NaN;
end

k = find(v2 > 2 & v1 > 0);
if any(k(:))
    m(k) = v2(k) ./ (v2(k) - 2);
end

% The mean is undefined for V2 less than or equal to 2.
k1 = find(v2 <= 2);
if any(k1(:))
    m(k1) = NaN;
end

k = find(v2 > 4 & v1 > 0);
if any(k(:))
    v(k) = m(k) .^ 2 * 2 .* (v1(k) + v2(k) - 2) ./ (v1(k) .* (v2(k) - 4));
end

% The variance is undefined for V2 less than or equal to 4.
k2 = find(v2 <= 4);
if any(k2(:))
    v(k2) = NaN;
end

% Special cases - v1 and/or v2 = Inf

% If v1 = Inf AND v2 is not Inf, then F(v1,v2) = v2/chi2(v2)
kinf1 = (v1 == Inf) & (v2 > 0 & v2 < Inf); 
% Mean defined only for v2 > 2
kinf1A = kinf1 & (v2 > 2); 
if any(kinf1A(:))
    m(kinf1A) = v2(kinf1A) ./ (v2(kinf1A) - 2);
end
clear kinf1A;
% Variance defined only for v2 > 4
kinf1B = kinf1 & (v2 > 4);
if any(kinf1B(:))
    v(kinf1B) = 2*(v2(kinf1B).^2) ./ ( (v2(kinf1B) - 2).^2 .* (v2(kinf1B) - 4) );
end
clear kinf1B;
clear kinf1;

% If v1 is not Inf AND v2 = Inf, then F(v1,v2) = chi2(v1)/v1
kinf2 = (v2 == Inf) & (v1 > 0 & v1 < Inf);
if any(kinf2(:))
    m(kinf2) = 1;
    v(kinf2) = 2 ./ v1(kinf2);
end
clear kinf2;

% If v1 is Inf and v2 = Inf, then F(v1,v2) concentrates around 1 with variance 0
kinf3 = (v2 == Inf) & (v1 == Inf);
if any(kinf3(:))
    m(kinf3) = 1;
    v(kinf3) = 0;
end
clear kinf3;
