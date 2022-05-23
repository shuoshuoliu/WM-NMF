function p = copulacdf(family,u,varargin)
%COPULACDF Cumulative probability function for a copula.
%   Y = COPULACDF('Gaussian',U,RHO) returns the cumulative probability of the
%   Gaussian copula with linear correlation parameters RHO, evaluated at the
%   points in U. U is an N-by-P matrix of values in [0,1], representing N
%   points in the P-dimensional unit hypercube.  RHO is a P-by-P correlation
%   matrix.  If U is an N-by-2 matrix, RHO may also be a scalar correlation
%   coefficient.
%
%   Y = COPULACDF('t',U,RHO,NU) returns the cumulative probability of the t
%   copula with linear correlation parameters RHO and degrees of freedom
%   parameter NU, evaluated at the points in U.  U is an N-by-P matrix of
%   values in [0,1]. RHO is a P-by-P correlation matrix.  If U is an N-by-2
%   matrix, RHO may also be a scalar correlation coefficient.
%   
%   Y = COPULACDF(FAMILY,U,ALPHA) returns the cumulative probability of the
%   bivariate Archimedean copula determined by FAMILY, with scalar parameter
%   ALPHA, evaluated at the points in U.  FAMILY is 'Clayton', 'Frank', or
%   'Gumbel'.  U is an N-by-2 matrix of values in [0,1].
%
%   Example:
%      u = linspace(0,1,10);
%      [U1,U2] = meshgrid(u,u);
%      F = copulacdf('Clayton',[U1(:) U2(:)],1);
%      surf(U1,U2,reshape(F,10,10));
%      xlabel('u1'); ylabel('u2');
%
%   See also COPULAPDF, COPULARND, COPULASTAT, COPULAPARAM.

%   Copyright 2005-2011 The MathWorks, Inc.


if nargin > 0
    family = convertStringsToChars(family);
end

if nargin < 3
    error(message('stats:copulacdf:WrongNumberOfInputs'));
end

[~,d] = size(u);
if d < 2
    error(message('stats:copulacdf:TooFewDimensions'));
end

% Map values outside of unit hypercube to the edges so that they get
% appropriate CDF values.
u(u<0) = 0; % doesn't overwrite NaNs
u(u>1) = 1;

families = {'gaussian','t','clayton','frank','gumbel'};
family = internal.stats.getParamVal(family,families,'FAMILY');

switch family

% Elliptical copulas -- these require a quadrature.
case {'gaussian' 't'}
    Rho = varargin{1};
    if (d == 2) && isscalar(Rho)
        if ~(-1 < Rho && Rho < 1)
            error(message('stats:copulacdf:BadScalarCorrelation'));
        end
        Rho = [1 Rho; Rho 1];
    elseif ~isequal(size(Rho), [d d]) || any(diag(Rho) ~= 1)
        error(message('stats:copulacdf:BadCorrelationMatrix'));
    end
    if isequal(family,'gaussian')
        % MVNCDF will check that Rho is symmetric and positive definite.
        p = mvncdf(norminv(u),zeros(1,d),Rho);
    else
        if nargin < 4
            error(message('stats:copulacdf:MissingDF'));
        end
        nu = varargin{2};
        if ~(isscalar(nu) && (0 < nu))
            error(message('stats:copulacdf:BadDegreesOfFreedom'));
        end
        % MVTCDF will check that Rho is symmetric and positive definite.
        p = mvtcdf(tinv(u,nu),Rho,nu);
    end

% one-parameter Archimedean copulas -- the CDF exists in closed form.
case {'clayton' 'frank' 'gumbel'}
    if d > 2
        error(message('stats:copulacdf:TooManyDimensions'));
    end
    alpha = varargin{1};
    if ~isscalar(alpha)
        error(message('stats:copulacdf:BadArchimedeanParameter'));
    end

    switch family
    case 'clayton' % a.k.a. Cook-Johnson
        % C(u1,u2) = (u1^(-alpha) + u2^(-alpha) - 1)^(-1/alpha)
        if alpha < 0
            error(message('stats:copulacdf:BadClaytonParameter'));
        elseif alpha == 0
            p = prod(u,2);
        else
            p = (sum(u.^(-alpha), 2) - 1) .^ (-1./alpha);
        end

    case 'frank'
        % C(u1,u2) = -(1/alpha)*log(1 + (exp(-alpha*u1)-1)*(exp(-alpha*u1)-1)/(exp(-alpha)-1))
        if alpha == 0
            p = prod(u,2);
        else
            p = -log((exp(-alpha) + (exp(-alpha.*sum(u,2)) - sum(exp(-alpha.*u),2))) ./ expm1(-alpha)) ./ alpha;
        end

    case 'gumbel' % a.k.a. Gumbel-Hougaard
        % C(u1,u2) = exp(-((-log(u1))^alpha + (-log(u2))^alpha)^(1/alpha))
        if alpha < 1
            error(message('stats:copulacdf:BadGumbelParameter'));
        elseif alpha == 1
            p = prod(u,2);
        else
            k = (u~=0); v = Inf(size(u)); v(k) = -log(u(k)); % avoid log(0) warnings
            v = sort(v,2); vmin = v(:,1); vmax = v(:,2); % min/max, but avoid dropping NaNs
            vmax(vmax==0) = realmin; vmin(vmin==Inf) = realmax; % avoid spurious NaNs from vmin/vmax
            p = exp(-vmax.*(1+(vmin./vmax).^alpha).^(1./alpha));
        end

    end
end
