function [varargout] = expcdf(x,varargin)
%EXPCDF Exponential cumulative distribution function.
%   P = EXPCDF(X,MU) returns the cdf of the exponential distribution with
%   mean parameter MU, evaluated at the values in X.  The size of P is
%   the common size of the input arguments.  A scalar input functions as a
%   constant matrix of the same size as the other input.
%
%   The default value for MU is 1.
%
%   [P,PLO,PUP] = EXPCDF(X,MU,PCOV,ALPHA) produces confidence bounds
%   for P when the input parameter MU is an estimate.  PCOV is the
%   variance of the estimated MU.  ALPHA has a default value of 0.05, and
%   specifies 100*(1-ALPHA)% confidence bounds.  PLO and PUP are arrays of
%   the same size as P containing the lower and upper confidence bounds.
%   The bounds are based on a normal approximation for the distribution of
%   the log of the estimate of MU.  You can get a more accurate set of
%   bounds simply by using EXPFIT to get a confidence interval for MU,
%   and evaluating EXPCDF at the lower and upper end points of that interval.
%
%   [...] = EXPCDF(...,'upper') returns the upper tail probability of 
%   the exponential distribution. 
%
%   See also CDF, EXPFIT, EXPINV, EXPLIKE, EXPPDF, EXPRND, EXPSTAT.

%   References:
%     [1] Lawless, J.F. (1982) Statistical Models and Methods for Lifetime Data, Wiley,
%         New York.
%     [2} Meeker, W.Q. and L.A. Escobar (1998) Statistical Methods for Reliability Data,
%         Wiley, New York.
%     [3] Crowder, M.J., A.C. Kimber, R.L. Smith, and T.J. Sweeting (1991) Statistical
%         Analysis of Reliability Data, Chapman and Hall, London.

%     Copyright 1993-2019 The MathWorks, Inc. 


if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin<1
    error(message('stats:expcdf:TooFewInputsX'));
end

if nargin>1 && strcmpi(varargin{end},'upper')
    %Compute upper tail and remove 'upper' flag
    uflag=true;
    varargin(end) = [];
elseif nargin>1 && ischar(varargin{end})&& ~strcmpi(varargin{end},'upper')
    error(message('stats:cdf:UpperTailProblem'));
else
    uflag=false;
end
[varargout{1:max(1,nargout)}] = localexpcdf(uflag,x,varargin{:});


function [p,plo,pup] = localexpcdf(uflag,x,mu,pcov,alpha) 
if nargin < 3
    mu = 1;
end

% More checking if we need to compute confidence bounds.
if nargout>1
   if nargin<4
      error(message('stats:expcdf:TooFewInputsVariance'));
   end
   if numel(pcov)~=1
      error(message('stats:expcdf:BadVarianceScalar'));
   end
   if nargin<5
      alpha = 0.05;
   elseif ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
      error(message('stats:expcdf:BadAlpha'));
   end
end

% Return NaN for out of range parameters.
mu(mu <= 0) = NaN;

try
    z = x./mu;
catch
    error(message('stats:expcdf:InputSizeMismatch'));
end

z(z < 0) = 0; % Force a zero for negative x.
if uflag == true
    p = exp(-z);
else    
    p = -expm1(-z);
end

% Compute confidence bounds if requested.
if nargout>=2
   pcov = double(pcov);
   alpha = double(alpha);

   % Work on log scale.
   logz = log(z);
   if pcov<0
      error(message('stats:expcdf:BadVarianceNonNeg'));
   end
   normz = -norminv(alpha/2);
   halfwidth = normz * sqrt(pcov ./ (mu.^2));
   zlo = logz - halfwidth;
   zup = logz + halfwidth;
   
   % Convert back to original scale
   if uflag == true
       plo = exp(-exp(zup));
       pup = exp(-exp(zlo));
   else
       plo = - expm1(-exp(zlo));
       pup = - expm1(-exp(zup));
   end
end
