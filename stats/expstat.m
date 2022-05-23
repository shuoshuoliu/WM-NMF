function [m,v] = expstat(mu)
%EXPSTAT Mean and variance of the exponential distribution.
%   [M,V] = EXPSTAT(MU) returns the mean and variance of the exponential
%   distribution with mean parameter MU.  M and V are the size of the
%   input argument.
%
%   See also EXPCDF, EXPFIT, EXPINV, EXPLIKE, EXPPDF, EXPRND.

%   References:
%     [1] Lawless, J.F. (1982) Statistical Models and Methods for Lifetime Data, Wiley,
%         New York.
%     [2} Meeker, W.Q. and L.A. Escobar (1998) Statistical Methods for Reliability Data,
%         Wiley, New York.
%     [3] Crowder, M.J., A.C. Kimber, R.L. Smith, and T.J. Sweeting (1991) Statistical
%         Analysis of Reliability Data, Chapman and Hall, London.

%     Copyright 1993-2009 The MathWorks, Inc. 


if nargin < 1
    error(message('stats:expstat:TooFewInputs'));
end

% Return NaN for out of range parameters.
mu(mu <= 0) = NaN;

m = mu;
v = mu.^2;
