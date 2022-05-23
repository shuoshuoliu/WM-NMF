function [m,v]= unifstat(a,b)
%UNIFSTAT Mean and variance of the continuous uniform distribution.
%   [M,V] = UNIFSTAT(A,B) returns the mean and variance of
%   the uniform distribution on the interval [A,B].
%
%   See also UNIFCDF, UNIFINV, UNIFIT, UNIFPDF, UNIFRND.

%   Reference:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.34.

%   Copyright 1993-2004 The MathWorks, Inc. 


if nargin < 2,
 error(message('stats:unifstat:TooFewInputs'));
end

[errorcode,a,b] = distchck(2,a,b);

if errorcode > 0
    error(message('stats:unifstat:InputSizeMismatch'));
end

m = (a + b) / 2;
v = (b - a) .^ 2 / 12;


% Return NaN if the lower limit is greater than the upper limit.
k1 = (a >= b);
if any(k1(:))
    m(k1) = NaN;
    v(k1) = NaN;
end

