function [x,xlo,xup] = gaminv(p,a,b,pcov,alpha)
%GAMINV Inverse of the gamma cumulative distribution function (cdf).
%   X = GAMINV(P,A,B) returns the inverse cdf for the gamma distribution
%   with shape A and scale B, evaluated at the values in P.  The size of X
%   is the common size of the input arguments.  A scalar input functions as
%   a constant matrix of the same size as the other inputs.
%
%   [X,XLO,XUP] = GAMINV(P,A,B,PCOV,ALPHA) produces confidence bounds for
%   X when the input parameters A and B are estimates. PCOV is a 2-by-2
%   matrix containing the covariance matrix of the estimated parameters.
%   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)%
%   confidence bounds.  XLO and XUP are arrays of the same size as X
%   containing the lower and upper confidence bounds.
%
%   See also GAMCDF, GAMFIT, GAMLIKE, GAMPDF, GAMRND, GAMSTAT.

%   GAMINV uses Newton's method to find roots of GAMCDF(X,A,B) = P.

%   References:
%      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
%          Functions, Dover, New York, section 26.1.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley.

%   Copyright 1993-2019 The MathWorks, Inc.

if nargin < 2
    error(message('stats:gaminv:TooFewInputs'));
elseif nargin < 3
    b = 1;
end
if ~isscalar(a) || ~isscalar(b)
    % skip this for the common case well handled here
    [p,a,b] = internal.stats.checkCompatibleSizes(p,a,b);
end

% More checking if we need to compute confidence bounds.
if nargout > 1
   if nargin < 4
      error(message('stats:gaminv:NeedCovariance'));
   end
   if ~isequal(size(pcov),[2 2])
      error(message('stats:gaminv:BadCovarianceSize'));
   end
   if nargin < 5
      alpha = 0.05;
   elseif ~isnumeric(alpha) || numel(alpha) ~= 1 || alpha <= 0 || alpha >= 1
      error(message('stats:gaminv:BadAlpha'));
   end
end

% Weed out any out of range parameters or edge/bad probabilities.
try
    okAB = (0 < a & a < Inf) & (0 < b);
    okP = (0 < p & p < 1);
    k = (okAB & okP); % GAMMAINCINV would handle 0 or 1, but the CI code won't
catch
    error(message('stats:gaminv:InputSizeMismatch'));
end
allOK = all(k(:));

% Fill in NaNs for out of range cases, fill in edges cases when P is 0 or 1.
if ~allOK
    x = NaN(size(k),'like',internal.stats.dominantType(p,a,b)); % single if p, a, or b is
    x(okAB & p == 0) = 0;
    x(okAB & p == 1) = Inf;

    x(a==0 & okP) = 0;
    
    if nargout > 1
        % Confidence intervals also depend on type of pcov & alpha
        ciType = internal.stats.dominantType(x,pcov,alpha);
        xlo = cast(x,"like",ciType); % NaNs or zeros or Infs
        xup = cast(x,"like",ciType); % NaNs or zeros or Infs
    end

    % Remove the bad/edge cases, leaving the easy cases.  If there's
    % nothing remaining, return.
    if any(k(:))
        if numel(p) > 1, p = p(k); end
        if numel(a) > 1, a = a(k); end
        if numel(b) > 1, b = b(k); end
    else
        return;
    end
end

% Call BETAINCINV to find a root of BETAINC(Q,A,B) = P
q = gammaincinv(p,a);

tolerance = sqrt(eps(ones('like',q)));
badcdf = ((abs(gammainc(q,a) - p)./p) > tolerance);
if any(badcdf(:))   % cdf is too far off
    didnt = find(badcdf, 1, 'first');
    if numel(a) == 1, abad = a; else, abad = a(didnt); end
    if numel(b) == 1, bbad = b; else, bbad = b(didnt); end
    if numel(p) == 1, pbad = p; else, pbad = p(didnt); end
    warning(message('stats:gaminv:NoConvergence', sprintf( '%g', abad ), sprintf( '%g', bbad ), sprintf( '%g', pbad )));
end

% Add in the scale factor, and broadcast the values to the correct place if
% need be.
if allOK
    x = q .* b;
else
    x(k) = q .* b;
end

% Compute confidence bounds if requested.
if nargout >= 2
    pcov = double(pcov);
    alpha = double(alpha);

    logq = log(q);
    dqda = -dgammainc(q,a) ./ exp((a-1).*logq - q - gammaln(a));

    % Approximate the variance of x=q*b on the log scale.
    %    dlogx/da = dlogx/dq * dqda = dqda/q
    %    dlogx/db = 1/b
    logx = logq + log(b);
    varlogx = pcov(1,1).*(dqda./q).^2 + 2.*pcov(1,2).*dqda./(b.*q) + pcov(2,2)./(b.^2);
    if any(varlogx(:) < 0)
        error(message('stats:gaminv:BadCovariancePosDef'));
    end
    z = -norminv(alpha/2);
    halfwidth = z * sqrt(varlogx);

    % Convert back to original scale
    if allOK
        xlo = exp(logx - halfwidth);
        xup = exp(logx + halfwidth);
    else
        xlo(k) = exp(logx - halfwidth);
        xup(k) = exp(logx + halfwidth);
    end
end
