function p = chi2cdf(x,v,uflag)
%CHI2CDF Chi-square cumulative distribution function.
%   P = CHI2CDF(X,V) returns the chi-square cumulative distribution
%   function with V degrees of freedom at the values in X.
%   The chi-square density function with V degrees of freedom,
%   is the same as a gamma density function with parameters V/2 and 2.
%
%   The size of P is the common size of X and V. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   P = CHI2CDF(X,V,'upper') returns the upper tail probability of 
%   the chi-square distribution. 
%
%   See also CHI2INV, CHI2PDF, CHI2RND, CHI2STAT, CDF.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.4.
%
%   Notice that we do not check if the degree of freedom parameter is integer
%   or not. In most cases, it should be an integer. Numerically, non-integer 
%   values still gives a numerical answer, thus, we keep them.

%   Copyright 1993-2013 The MathWorks, Inc. 


if nargin > 2
    uflag = convertStringsToChars(uflag);
end

if nargin < 2, 
    error(message('stats:chi2cdf:TooFewInputs'));
end

[errorcode x v] = distchck(2,x,v);

if errorcode > 0
    error(message('stats:chi2cdf:InputSizeMismatch'));
end

if nargin>2
    if ~strcmpi(uflag,'upper')
        error(message('stats:cdf:UpperTailProblem'));
    end
else
    uflag=[];
end

% Call the gamma distribution function. 
p = gamcdf(x,v/2,2,uflag);
