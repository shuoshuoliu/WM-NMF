function [varargout] = wblcdf(x,varargin)
%WBLCDF Weibull cumulative distribution function (cdf).
%   P = WBLCDF(X,A,B) returns the cdf of the Weibull distribution
%   with scale parameter A and shape parameter B, evaluated at the
%   values in X.  The size of P is the common size of the input arguments.
%   A scalar input functions as a constant matrix of the same size as the
%   other inputs.
%
%   Default values for A and B are 1 and 1, respectively.
%
%   [P,PLO,PUP] = WBLCDF(X,A,B,PCOV,ALPHA) produces confidence
%   bounds for P when the input parameters A and B are estimates.
%   PCOV is a 2-by-2 matrix containing the covariance matrix of the estimated
%   parameters.  ALPHA has a default value of 0.05, and specifies
%   100*(1-ALPHA)% confidence bounds.  PLO and PUP are arrays of the same
%   size as P containing the lower and upper confidence bounds.
%
%   [...] = WBLCDF(...,'upper') returns the upper tail probability of 
%   the Weibull distribution. 
%
%   See also CDF, WBLFIT, WBLINV, WBLLIKE, WBLPDF, WBLRND, WBLSTAT.

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
    error(message('MATLAB:narginchk:notEnoughInputs'))
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
[varargout{1:max(1,nargout)}] = localwblcdf(uflag,x,varargin{:});


function [p,plo,pup] = localwblcdf(uflag,x,A,B,pcov,alpha)

if nargin < 3
    A = 1;
end
if nargin < 4
    B = 1;
end

% More checking if we need to compute confidence bounds.
if nargout>1
   if nargin<5
      error(message('stats:wblcdf:MustProvideCovMatrix'));
   end
   if ~isequal(size(pcov),[2 2])
      error(message('stats:wblcdf:BadCovariance'));
   end
   if nargin<6
      alpha = 0.05;
   elseif ~isnumeric(alpha) || numel(alpha)~=1 || alpha<=0 || alpha>=1
      error(message('stats:wblcdf:BadAlpha'));
   end
end

% Return NaN for out of range parameters.
A(A <= 0) = NaN;
B(B <= 0) = NaN;

% Force a zero for negative x.
x(x < 0) = 0;

try
    z = (x./A).^B;
catch me 
    error(message('stats:wblcdf:InputSizeMismatch'));
end

if uflag==true
    p = exp(-z);
else
    p = -expm1(-z);
end

% Compute confidence bounds if requested.
if nargout>=2
   pcov = double(pcov);
   alpha = double(alpha);

   % Work on extreme value scale (log scale).
   logz = log(z);
   dA = 1./A;
   dB = -1./(B.^2);
   logzvar = (pcov(1,1).*dA.^2 + 2*pcov(1,2).*dA.*dB.*logz + pcov(2,2).*(dB.*logz).^2) .* (B.^2);
%    deriv = [1./A, -1./(B.^2)];
%    pcov = pcov .* (deriv' * deriv);
%    logzvar = (pcov(1,1) + 2*pcov(1,2)*logz + pcov(2,2)*logz.^2) .* (B.^2);
   if any(logzvar<0)
      error(message('stats:wblcdf:BadCovariance'));
   end
   normz = -norminv(alpha/2);
   halfwidth = normz * sqrt(logzvar);
   zlo = logz - halfwidth;
   zup = logz + halfwidth;

   % Convert back to Weibull scale
   if uflag == true
       plo = exp(-exp(zup));
       pup = exp(-exp(zlo));
   else
       plo = -expm1(-exp(zlo));
       pup = -expm1(-exp(zup));
   end
end
