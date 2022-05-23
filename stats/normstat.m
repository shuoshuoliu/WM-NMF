function [m,v]= normstat(mu,sigma)
%NORMSTAT Mean and variance for the normal distribution.
%   [M,V] = NORMSTAT(MU,SIGMA) returns the mean and variance of the normal
%   distribution with mean MU and standard deviation SIGMA.  The sizes of M
%   and V are the common size of the input arguments.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also NORMCDF, NORMFIT, NORMINV, NORMLIKE, NORMPDF, NORMRND.

%   References:
%      [1] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2019 The MathWorks, Inc. 


if nargin < 2
    error(message('stats:normstat:TooFewInputs'));
end

try
    prototype = mu+sigma;
catch ME
    if strcmp(ME.identifier,'MATLAB:dimagree')
        error(message('stats:normstat:InputSizeMismatch'));
    else
        rethrow(ME)
    end
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

m = prototype;
v = prototype;

m(:) = mu + zeros(size(sigma)); % expand m's size to match sigma if necessary
v(:) = sigma.^2 + zeros(size(mu)); % expand v's size to match mu if necessary

m(isnan(v)) = NaN;
v(isnan(m)) = NaN;
