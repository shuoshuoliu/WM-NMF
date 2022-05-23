function [varargout] = gpcdf(x,varargin)
%GPCDF Generalized Pareto cumulative distribution function (cdf).
%   P = GPCDF(X,K,SIGMA,THETA) returns the cdf of the generalized Pareto (GP)
%   distribution with tail index (shape) parameter K, scale parameter SIGMA,
%   and threshold (location) parameter THETA, evaluated at the values in X.
%   The size of P is the common size of the input arguments.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for K, SIGMA, and THETA are 0, 1, and 0, respectively.
%
%   When K = 0 and THETA = 0, the GP is equivalent to the exponential
%   distribution.  When K > 0 and THETA = SIGMA/K, the GP is equivalent to the
%   Pareto distribution.  The mean of the GP is not finite when K >= 1, and the
%   variance is not finite when K >= 1/2.  When K >= 0, the GP has positive
%   density for X>THETA, or, when K < 0, for 0 <= (X-THETA)/SIGMA <= -1/K.
%
%   [...] = GPCCDF(...,'upper') returns the upper tail probability of 
%   the generalized Pareto distribution. 
%
%   See also GPFIT, GPINV, GPLIKE, GPPDF, GPRND, GPSTAT, CDF.

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
    error(message('stats:gpcdf:TooFewInputs'));
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
[varargout{1:max(1,nargout)}] = localgpcdf(x,uflag,varargin{:});

function p = localgpcdf(x,uflag,k,sigma,theta)
if nargin < 3 || isempty(k), k = 0;     end
if nargin < 4 || isempty(sigma), sigma = 1; end
if nargin < 5 || isempty(theta), theta = 0; end

[err,sizeOut] = internal.stats.statsizechk(4,x,k,sigma,theta);
if err > 0
    error(message('stats:gpcdf:InputSizeMismatch'));
end
if isscalar(k), k = repmat(k,sizeOut); end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

p = zeros(sizeOut,superiorfloat(x,k,sigma,theta));

% Support is 0 <= (x-theta)/sigma, force zero below that
z = (x-theta) ./ sigma; z(z<0) = 0; % max drops NaNs
if isscalar(z), z = repmat(z,sizeOut); end

% Find the k==0 cases and fill them in.
j = (abs(k) < eps);

if uflag == true
    p(j) = exp(-z(j));
else
    p(j) = - expm1(-z(j));
end

% When k<0, the support is 0 <= x/sigma <= -1/k.
t = z.*k;
jj = (t<=-1 & k<-eps);
t(jj) = 0; % temporarily silence warnings from log1p

% Find the k~=0 cases and fill them in.
j = ~j;
if uflag==true
    p(j) = exp((-1./k(j)).*log1p(t(j)));  
else
    p(j) = -expm1((-1./k(j)).*log1p(t(j))); % 1 - (1 + z.*k).^(-1./k)
end
if uflag==true
    p(jj) = 0;
else
    p(jj) = 1;
end