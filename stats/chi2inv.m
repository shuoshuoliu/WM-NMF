function x = chi2inv(p,v)
%CHI2INV Inverse of the chi-square cumulative distribution function (cdf).
%   X = CHI2INV(P,V)  returns the inverse of the chi-square cdf with V  
%   degrees of freedom at the values in P. The chi-square cdf with V 
%   degrees of freedom, is the gamma cdf with parameters V/2 and 2.   
%
%   The size of X is the common size of P and V. A scalar input
%   functions as a constant matrix of the same size as the other input.   
%
%   See also CHI2CDF, CHI2PDF, CHI2RND, CHI2STAT, ICDF.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.4.
%      [2] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, section 10.2 (page 144)

%   Copyright 1993-2019 The MathWorks, Inc.


if nargin < 2
    error(message('stats:chi2inv:TooFewInputs'));
end

errorcode = distchck(2,p,v);

if errorcode > 0
    error(message('stats:chi2inv:InputSizeMismatch'));
end

% Call the gamma inverse function. 
x = gaminv(p,v/2,2);
