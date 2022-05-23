function [lambdahat, lambdaci] = poissfit(x,alpha)
%POISSFIT Parameter estimates and confidence intervals for Poisson data.
%   POISSFIT(X) Returns the estimate of the parameter of the Poisson
%   distribution give the data X. 
%
%   [LAMBDAHAT, LAMBDACI] = POISSFIT(X,ALPHA) gives MLEs and 100(1-ALPHA) 
%   percent confidence intervals given the data. By default, the
%   optional parameter ALPHA = 0.05 corresponding to 95% confidence intervals.
%
%   See also POISSCDF, POISSINV, POISSPDF, POISSRND, POISSTAT, MLE. 

%   Copyright 1993-2019 The MathWorks, Inc. 


if nargin < 2 
    alpha = 0.05;
else
    alpha = double(alpha);
end

% Initialize params to zero.
[m, n] = size(x);
if min(m,n) == 1
   x = x(:);
   m = max(m,n);
   n = 1;
end

if any(x<0)
    error(message('stats:poissfit:InvalidX'))
end

lambdahat = mean(x);
if ~isfloat(lambdahat)
   lambdahat = double(lambdahat);
end

if nargout > 1
   lambdaci = statpoisci(m,lambdahat,alpha);
end

