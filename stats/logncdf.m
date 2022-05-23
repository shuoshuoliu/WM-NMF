function [varargout] = logncdf(x,varargin)
%LOGNCDF Lognormal cumulative distribution function (cdf).
%   P = LOGNCDF(X,MU,SIGMA) returns values at X of the lognormal cdf with 
%   distribution parameters MU and SIGMA.  MU and SIGMA are the mean and 
%   standard deviation, respectively, of the associated normal distribution.  
%   The size of P is the common size of X, MU and SIGMA.  A scalar input 
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   [P,PLO,PUP] = LOGNCDF(X,MU,SIGMA,PCOV,ALPHA) returns confidence bounds
%   for P when the input parameters MU and SIGMA are estimates.  PCOV is a
%   2-by-2 matrix containing the covariance matrix of the estimated parameters.
%   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)% confidence
%   bounds.  PLO and PUP are arrays of the same size as P containing the lower
%   and upper confidence bounds.
%
%   [...] = LOGNCDF(...,'upper') returns the upper tail probability of 
%   the lognormal distribution. 
%
%   See also ERF, ERFC, LOGNFIT, LOGNINV, LOGNLIKE, LOGNPDF, LOGNRND, LOGNSTAT.

%   References:
%      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
%          Functions, Dover, New York, 1046pp., sections 7.1, 26.2.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2019 The MathWorks, Inc.


if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin<1
   error(message('stats:logncdf:TooFewInputsX'));
end
if nargin>1 && strcmpi(varargin{end},'upper')
    uflag=true;
    varargin(end)= [];
elseif nargin>1 && ischar(varargin{end})&& ~strcmpi(varargin{end},'upper')
    error(message('stats:cdf:UpperTailProblem'));
else
    uflag=false;
end
[varargout{1:max(1,nargout)}] = locallogncdf(uflag,x,varargin{:});

function [p,plo,pup]=locallogncdf(uflag,x,mu,sigma,pcov,alpha)
if nargin < 3
    mu = 0;
end
if nargin < 4
    sigma = 1;
end

% More checking if we need to compute confidence bounds.
if nargout>1
   if nargin<5
      error(message('stats:logncdf:TooFewInputsCovariance'));
   end
   if ~isequal(size(pcov),[2 2])
      error(message('stats:logncdf:BadCovarianceSize'));
   end
   if nargin<6
      alpha = 0.05;
   elseif ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
      error(message('stats:logncdf:BadAlpha'));
   end
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

% Negative data would create complex values, which erfc cannot handle.
x(x < 0) = 0;

try
    z = (log(x)-mu) ./ sigma;
    if uflag == true
        z = -z;
    end
catch ME
    error(message('stats:logncdf:InputSizeMismatch'));
end

% Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
% to produce accurate near-zero results for large negative x.
p = 0.5 * erfc(-z ./ sqrt(2));

% Compute confidence bounds if requested.
if nargout>=2
   pcov = double(pcov);
   alpha = double(alpha);

   zvar = (pcov(1,1) + 2*pcov(1,2)*z + pcov(2,2)*z.^2) ./ (sigma.^2);
   if any(zvar(:)<0)
      error(message('stats:logncdf:BadCovarianceSymPos'));
   end
   normz = -norminv(alpha/2);
   halfwidth = normz * sqrt(zvar);
   zlo = z - halfwidth;
   zup = z + halfwidth;

   plo = 0.5 * erfc(-zlo./sqrt(2));
   pup = 0.5 * erfc(-zup./sqrt(2));
end
