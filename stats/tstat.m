function [mn,v]= tstat(nu)
%TSTAT  Mean and variance for the student's t distribution.
%   [MN,V] = TSTAT(NU) returns the mean and variance of
%   Student's T distribution with NU degrees of freedom.
%
%   See also TCDF, TINV, TPDF, TRND.

%   References:
%      [1]  E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, New York, 1970, Section 10.3, pages 144-146.

%   Copyright 1993-2014 The MathWorks, Inc. 


if nargin < 1,   
    error(message('stats:tstat:TooFewInputs'));       
end

% Initialize the mean and variance to zero.
mn = zeros(size(nu),'like',nu);

k = find(nu <= 0);
if any(k)
    mn(k) = NaN;
end
v = mn;

% The mean of the t distribution is zero unless there is only one 
% degree of freedom. In that case the mean does not exist.
mn(nu <= 1) = NaN;

% The variance of the t distribution is undefined for one and two
% degrees of freedom.
v(nu <= 2) = NaN;

k = (nu > 2);
v(k) = nu(k) ./ (nu(k) - 2);

