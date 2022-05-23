function [nlogL,acov] = gevlike(params, data)
%GEVLIKE Negative log-likelihood for the generalized extreme value distribution.
%   NLOGL = GEVLIKE(PARAMS,DATA) returns the negative of the log-likelihood
%   for the generalized extreme value (GEV) distribution, evaluated at
%   parameters PARAMS(1) = K, PARAMS(2) = SIGMA, and PARAMS(3) = MU, given
%   DATA.  NLOGL is a scalar.
%
%   [NLOGL, ACOV] = GEVLIKE(PARAMS,DATA) returns the inverse of Fisher's
%   information matrix, ACOV.  If the input parameter values in PARAMS are the
%   maximum likelihood estimates, the diagonal elements of ACOV are their
%   asymptotic variances.  ACOV is based on the observed Fisher's information,
%   not the expected information.
%
%   When K < 0, the GEV is the type III extreme value distribution.  When K >
%   0, the GEV distribution is the type II, or Frechet, extreme value
%   distribution.  If W has a Weibull distribution as computed by the WBLLIKE
%   function, then -W has a type III extreme value distribution and 1/W has a
%   type II extreme value distribution.  In the limit as K approaches 0, the
%   GEV is the mirror image of the type I extreme value distribution as
%   computed by the EVLIKE function.
%
%   The mean of the GEV distribution is not finite when K >= 1, and the
%   variance is not finite when K >= 1/2.  The GEV distribution has positive
%   density only for values of X such that K*(X-MU)/SIGMA > -1.
%
%   See also EVLIKE, GEVCDF, GEVFIT, GEVINV, GEVPDF, GEVRND, GEVSTAT.

%   References:
%      [1] Embrechts, P., C. Klüppelberg, and T. Mikosch (1997) Modelling
%          Extremal Events for Insurance and Finance, Springer.
%      [2] Kotz, S. and S. Nadarajah (2001) Extreme Value Distributions:
%          Theory and Applications, World Scientific Publishing Company.

%   Copyright 1993-2011 The MathWorks, Inc.


if nargin < 2
    error(message('stats:gevlike:TooFewInputs'));
elseif ~isvector(data)
    error(message('stats:gevlike:VectorRequired'));
end
if numel(params)~=3
    error(message('stats:probdists:WrongParameterLength',3));
end

k     = params(1);
sigma = params(2);
lnsigma = log(sigma);
mu    = params(3);

n = numel(data);
z = (data - mu) ./ sigma;

if abs(k) > eps
    u = 1 + k.*z;
    if min(u) > 0
        lnu = log1p(k.*z); % log(1 + k.*z)
        t = exp(-(1/k)*lnu); % (1 + k.*z).^(-1/k)
        nlogL = n*lnsigma + sum(t) + (1+1/k)*sum(lnu);
        if nargout > 1
            s = expm1(-(1/k)*lnu); % (1 + k.*z).^(-1/k) - 1
            r = (s - k)./u;
            q = (t - k)./u.^2;
            d2k2 = sum(2.*(z./u - lnu./k).*s + (k+1).*z.^2.*q - 2.*lnu.*z.*t./(u.*k) + lnu.^2.*t./k.^2) ./ k.^2;
            d2ksigma = sum(z.*(z - 1 - (k+1).*z.*t./k + (k.*z+1).*lnu.*t./k.^2)./u.^2) ./ sigma;
            d2kmu = sum((z - 1 - (k+1).*z.*t./k + (k.*z+1).*lnu.*t./k.^2)./u.^2) ./ sigma;
            d2sigma2 = sum((k+1).*z.^2.*q - 2.*z.*r - 1) ./ sigma.^2;
            d2sigmamu = sum(z.*t./u.^2 - r./u) ./ sigma.^2;
            d2mu2 = sum((k+1).*q) ./ sigma.^2;
            nhess = [d2k2 d2ksigma d2kmu; d2ksigma d2sigma2 d2sigmamu; d2kmu d2sigmamu d2mu2];
        end
    else
        % The support of the GEV is 1+k*z > 0, or x > mu - sigma/k.
        nlogL = NaN;
        if nargout > 1
            nhess = NaN(3,3);
        end
    end
else % limiting extreme value dist'n as k->0
    t = exp(-z);
    nlogL = n*lnsigma + sum(t + z);
    if nargout > 1
        s = expm1(-z); % exp(-z) - 1
        d2k2 = sum(t.*z.^4/4 - 2*s.*z.^3 - z.^2);
        d2ksigma = sum((z.^2.*t./2 - z.*s - 1) .* z) ./ sigma;
        d2kmu = sum(z.^2.*t/2 - z.*s) ./ sigma;
        d2sigma2 = sum((z.*(z+2).*t - ((2*z+1)))) ./ sigma.^2;
        d2sigmamu = sum(((z+1).*t - 1)) ./ sigma.^2;
        d2mu2 = sum(t) ./ sigma.^2;
        nhess = [d2k2 d2ksigma d2kmu; d2ksigma d2sigma2 d2sigmamu; d2kmu d2sigmamu d2mu2];
    end
end

if nargout > 1
    % The asymptotic cov matrix approximation is the negative inverse
    % of the hessian.
    acov = inv(nhess);
end
