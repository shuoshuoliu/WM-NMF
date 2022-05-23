function y = chi2pdf(x,v)
%CHI2PDF Chi-square probability density function (pdf).
%   Y = CHI2PDF(X,V) returns the chi-square pdf with V degrees  
%   of freedom at the values in X. The chi-square pdf with V 
%   degrees of freedom, is the gamma pdf with parameters V/2 and 2.
%
%   The size of Y is the common size of X and V. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   See also CHI2CDF, CHI2INV, CHI2RND, CHI2STAT, PDF.

%   Notice that we do not check if the degree of freedom parameter is integer
%   or not. In most cases, it should be an integer. Numerically, non-integer
%   values still gives a numerical answer, thus, we keep them.

%   References:
%      [1]  L. Devroye, "Non-Uniform Random Variate Generation", 
%      Springer-Verlag, 1986, pages 402-403.

%   Copyright 1993-2019 The MathWorks, Inc. 


if nargin < 2
    error(message('stats:chi2pdf:TooFewInputs')); 
end

[errorcode, x, v] = distchck(2,x,v);

if errorcode > 0
    error(message('stats:chi2pdf:InputSizeMismatch'));
end

y = gampdf(x,v/2,2);
