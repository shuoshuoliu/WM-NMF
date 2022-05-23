function [varargout] = gevcdf(x,varargin)
%GEVCDF Generalized extreme value cumulative distribution function (cdf).
%   P = GEVCDF(X,K,SIGMA,MU) returns the cdf of the generalized extreme value
%   (GEV) distribution with shape parameter K, scale parameter SIGMA, and
%   location parameter MU, evaluated at the values in X.  The size of P is the
%   common size of the input arguments.  A scalar input functions as a
%   constant matrix of the same size as the other inputs.
%
%   Default values for K, SIGMA, and MU are 0, 1, and 0, respectively.
%
%   When K < 0, the GEV is the type III extreme value distribution.  When K >
%   0, the GEV distribution is the type II, or Frechet, extreme value
%   distribution.  If W has a Weibull distribution as computed by the WBLCDF
%   function, then -W has a type III extreme value distribution and 1/W has a
%   type II extreme value distribution.  In the limit as K approaches 0, the
%   GEV is the mirror image of the type I extreme value distribution as
%   computed by the EVCDF function.
%
%   The mean of the GEV distribution is not finite when K >= 1, and the
%   variance is not finite when K >= 1/2.  The GEV distribution has positive
%   density only for values of X such that K*(X-MU)/SIGMA > -1.
%
%   [...] = GEVCDF(...,'upper') returns the upper tail probability of 
%   the generalized extreme value distribution.
%
%   See also EVCDF, GEVFIT, GEVINV, GEVLIKE, GEVPDF, GEVRND, GEVSTAT, CDF.

%   References:
%      [1] Embrechts, P., C. Klüppelberg, and T. Mikosch (1997) Modelling
%          Extremal Events for Insurance and Finance, Springer.
%      [2] Kotz, S. and S. Nadarajah (2001) Extreme Value Distributions:
%          Theory and Applications, World Scientific Publishing Company.

%   Copyright 1993-2018 The MathWorks, Inc.


if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 1
    error(message('stats:gevcdf:TooFewInputs'));
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
[varargout{1:max(1,nargout)}] = localgevcdf(x,uflag,varargin{:});

function p = localgevcdf(x,uflag,k,sigma,mu)
if nargin < 3, k = 0;     end
if nargin < 4, sigma = 1; end
if nargin < 5, mu = 0;    end

[err,sizeOut] = internal.stats.statsizechk(4,x,k,sigma,mu);
if err > 0
    error(message('stats:gevcdf:InputSizeMismatch'));
end
if isscalar(k), k = repmat(k,sizeOut); end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

p = zeros(sizeOut,superiorfloat(x,k,sigma,mu));

z = (x-mu)./sigma;
if isscalar(z), z = repmat(z,sizeOut); end

% Find the k==0 cases and fill them in.
j = (abs(k) < eps);
if uflag == true
    p(j) = -expm1(-exp(-z(j)));
else
    p(j) = exp(-exp(-z(j)));
end

% Find the k~=0 cases and fill them in.
j = ~j;
t = z.*k;

% 1 + k.*(x-mu)/sigma > 0 is the support, we'll force 0 or 1 outside that.
jj = j & (t<=-1);
t(jj) = 0; % temporarily silence warnings from log1p
if uflag == true
    p(j) = -expm1(-exp(-(1./k(j)).*log1p(t(j)))); 
else
    p(j) = exp(-exp(-(1./k(j)).*log1p(t(j))));% exp(-(1 + k.*z).^(-1./k));
end
if uflag == true
    p(jj) = (k(jj) >= 0);
else
    p(jj) = (k(jj) < 0);
end
