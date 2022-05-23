function [varargout] = evcdf(x,varargin)
%EVCDF Extreme value cumulative distribution function (cdf).
%   P = EVCDF(X,MU,SIGMA) returns the cdf of the type 1 extreme value
%   distribution with location parameter MU and scale parameter SIGMA,
%   evaluated at the values in X.  The size of P is the common size of the
%   input arguments.  A scalar input functions as a constant matrix of the
%   same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   [P,PLO,PUP] = EVCDF(X,MU,SIGMA,PCOV,ALPHA) produces confidence bounds
%   for P when the input parameters MU and SIGMA are estimates.  PCOV is a
%   2-by-2 matrix containing the covariance matrix of the estimated parameters.
%   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)% confidence
%   bounds.  PLO and PUP are arrays of the same size as P containing the lower
%   and upper confidence bounds.
%
%   The type 1 extreme value distribution is also known as the Gumbel
%   distribution.  The version used here is suitable for modeling minima; the
%   mirror image of this distribution can be used to model maxima by negating
%   X.  If Y has a Weibull distribution, then X=log(Y) has the type 1 extreme
%   value distribution.
%
%   [...] = EVCDF(...,'upper') returns the upper tail probability of 
%   the extreme value distribution. 
%
%   See also CDF, EVFIT, EVINV, EVLIKE, EVPDF, EVRND, EVSTAT.

%   References:
%     [1] Lawless, J.F. (1982) Statistical Models and Methods for Lifetime Data, Wiley,
%         New York.
%     [2} Meeker, W.Q. and L.A. Escobar (1998) Statistical Methods for Reliability Data,
%         Wiley, New York.
%     [3] Crowder, M.J., A.C. Kimber, R.L. Smith, and T.J. Sweeting (1991) Statistical
%         Analysis of Reliability Data, Chapman and Hall, London

%   Copyright 1993-2019 The MathWorks, Inc.


if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin<1
    error(message('stats:evcdf:TooFewInputs'));
end

if nargin>1 && strcmpi(varargin{end},'upper')
    %Compute upper tail
    uflag=true;
    varargin(end) = [];
elseif nargin>1 && ischar(varargin{end})&& ~strcmpi(varargin{end},'upper')
    error(message('stats:cdf:UpperTailProblem'));
else
    uflag=false;
end
[varargout{1:max(1,nargout)}] = localevcdf(uflag,x,varargin{:});


function [p,plo,pup] = localevcdf(uflag,x,mu,sigma,pcov,alpha)
if nargin < 3
    mu = 0;
end
if nargin < 4
    sigma = 1;
end

% More checking if we need to compute confidence bounds.
if nargout>1
   if nargin<5
      error(message('stats:evcdf:MissedCov'));
   end
   if ~isequal(size(pcov),[2 2])
      error(message('stats:evcdf:BadShapedCovariance'));
   end
   if nargin<6
      alpha = 0.05;
   elseif ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
      error(message('stats:evcdf:BadAlpha'));
   end
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    z = (x-mu)./sigma;
    if uflag == true
        p = exp( -exp(z) );
    else
        p = -expm1( -exp(z) );
    end
catch
    error(message('stats:evcdf:InputSizeMismatch'));
end

% Compute confidence bounds if requested.
if nargout>=2
   pcov = double(pcov);
   alpha = double(alpha);

   zvar = (pcov(1,1) + 2*pcov(1,2)*z + pcov(2,2)*z.^2) ./ (sigma.^2);
   if any(zvar<0)
      error(message('stats:evcdf:BadCovariance'));
   end
   normz = -norminv(alpha/2);
   halfwidth = normz * sqrt(zvar);
   zlo = z - halfwidth;
   zup = z + halfwidth;
   
   if uflag == true
       plo = exp(-exp(zup));
       pup = exp(-exp(zlo));
   else
       plo = -expm1(-exp(zlo));
       pup = -expm1(-exp(zup));
   end
end
