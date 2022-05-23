function y = copulapdf(family,u,varargin)
%COPULAPDF Probability density function for a copula.
%   Y = COPULAPDF('Gaussian',U,RHO) returns the probability density of the
%   Gaussian copula with linear correlation parameters RHO, evaluated at the
%   points in U. U is an N-by-P matrix of values in [0,1], representing N
%   points in the P-dimensional unit hypercube.  RHO is a P-by-P correlation
%   matrix.  If U is an N-by-2 matrix, RHO may also be a scalar correlation
%   coefficient.
%
%   Y = COPULAPDF('t',U,RHO,NU) returns the probability density of the t
%   copula with linear correlation parameters RHO and degrees of freedom
%   parameter NU, evaluated at the points in U.  U is an N-by-P matrix of
%   values in [0,1]. RHO is a P-by-P correlation matrix.  If U is an N-by-2
%   matrix, RHO may also be a scalar correlation coefficient.
%   
%   Y = COPULAPDF(FAMILY,U,ALPHA) returns the probability density of the
%   bivariate Archimedean copula determined by FAMILY, with scalar parameter
%   ALPHA, evaluated at the points in U.  FAMILY is 'Clayton', 'Frank', or
%   'Gumbel'.  U is an N-by-2 matrix of values in [0,1].
%
%   Example:
%      u = linspace(0,1,10);
%      [U1,U2] = meshgrid(u,u);
%      F = copulapdf('Clayton',[U1(:) U2(:)],1);
%      surf(U1,U2,reshape(F,10,10));
%      xlabel('u1'); ylabel('u2');
%
%   See also COPULACDF, COPULARND, COPULASTAT, COPULAPARAM.

%   Copyright 2005-2014 The MathWorks, Inc.


if nargin > 0
    family = convertStringsToChars(family);
end

if nargin < 3
    error(message('stats:copulapdf:WrongNumberOfInputs'));
end
    
[n,d] = size(u);
if d < 2
    error(message('stats:copulapdf:TooFewDimensions'));
end

% Values outside of unit hypercube will get zero density.  For now, change
% them to harmless values.
outOfRange = any(u<=0 | 1<=u, 2); % does not pick up NaNs
anyOutOfRange = any(outOfRange);
if anyOutOfRange
    % Replace entire rows that include out-of-range values with 0.5, but
    % don't overwrite NaN. This will help protect against NaN caused by
    % data close to the edge of the space. Any NaN found in the output
    % should then be the result of NaN input or bad parameter values.
    out = u(outOfRange,:);
    out(~isnan(out)) = 0.5;
    u(outOfRange,:) = out;
end

families = {'gaussian','t','clayton','frank','gumbel'};
family = internal.stats.getParamVal(family,families,'FAMILY');

switch family

% Elliptical copulas
%
% The CDF for u ~ Unif(0,1) for these copulas cannot be expressed in closed
% form, but the PDF for u is just the PDF for x ~ MVN (or MVT), times the
% determinant of the Jacobian of the transformation u = normcdf(x) (or
% tcdf(x)).
case {'gaussian' 't'}
    Rho = varargin{1};
    if (d == 2) && isscalar(Rho)
        if ~(-1 < Rho && Rho < 1)
            error(message('stats:copulapdf:BadScalarCorrelation'));
        end
        Rho = [1 Rho; Rho 1];
    elseif ~isequal(size(Rho), [d d]) || any(diag(Rho) ~= 1)
        error(message('stats:copulapdf:BadCorrelationMatrix'));
    end
    [R,err] = cholcov(Rho,0);
    if (err ~= 0)
        error(message('stats:copulapdf:SingularCorrelationMatrix'));
    end
    if isequal(family,'gaussian')
        x = norminv(u);
        % This is mvnpdf(x,zeros(1,d),Rho) ./ prod(normpdf(x), 2), but calculated
        % so that underflow is less likely.
        logSqrtDetRho = sum(log(diag(R)));
        z = x/R;
        y = exp(-0.5 .* sum(z.^2 - x.^2,2) - logSqrtDetRho);
    else
        if nargin < 4
            error(message('stats:copulapdf:MissingDF'));
        end
        nu = varargin{2};
        if ~(isscalar(nu) && (0 < nu))
            error(message('stats:copulapdf:BadDegreesOfFreedom'));
        end
        t = tinv(u,nu);
        % This is mvtpdf(t,Rho,nu) ./ prod(tpdf(t,nu), 2), but calculated so that
        % underflow is less likely.
        z = t/R;
        logSqrtDetRho = sum(log(diag(R)));
        const = gammaln((nu+d)./2) + (d-1).*gammaln(nu./2) - d.*gammaln((nu+1)./2) - logSqrtDetRho;
        numer = -((nu+d)./2) .* log(1 + sum(z.^2,2)./nu);
        denom = sum(-((nu+1)./2) .* log(1 + (t.^2)./nu), 2);
        y = exp(const + numer - denom);
    end
        
% one-parameter Archimedean copulas
%
% The CDF for u has a closed form for these copulas, and the PDF is just
% the mixed partial derivative d2C/du1du2.
case {'clayton' 'frank' 'gumbel'}
    if d > 2
        error(message('stats:copulapdf:TooManyDimensions'));
    end
    alpha = varargin{1};
    if ~isscalar(alpha)
        error(message('stats:copulapdf:BadArchimedeanParameter'));
    end

    switch family
    case 'clayton' % a.k.a. Cook-Johnson
        % C(u1,u2) = (u1^(-alpha) + u2^(-alpha) - 1)^(-1/alpha)
        if alpha < 0
            error(message('stats:copulapdf:BadClaytonParameter'));
        elseif alpha == 0
            y = ones(n,1);
        else
            logC = (-1./alpha).*log(sum(u.^(-alpha), 2) - 1);
            y = (alpha+1) .* exp((2.*alpha+1).*logC - sum((alpha+1).*log(u), 2));
        end

    case 'frank'
        % C(u1,u2) = -(1/alpha)*log(1 + (exp(-alpha*u1)-1)*(exp(-alpha*u1)-1)/(exp(-alpha)-1))
        if alpha == 0
            y = ones(n,1);
        else
             y = alpha*(1-exp(-alpha))./ ...
                 (cosh(alpha*diff(u,1,2)/2)*2 - exp(alpha*(sum(u,2)-2)/2) - exp(-alpha*sum(u,2)/2)).^2;
        end

    case 'gumbel' % a.k.a. Gumbel-Hougaard
        % C(u1,u2) = exp(-((-log(u1))^alpha + (-log(u2))^alpha)^(1/alpha))
        if alpha < 1
            error(message('stats:copulapdf:BadGumbelParameter'));
        elseif alpha == 1
            y = ones(n,1);
        else
            v = -log(u); % u is strictly in (0,1) => v strictly in (0,Inf)
            v = sort(v,2); vmin = v(:,1); vmax = v(:,2); % min/max, but avoid dropping NaNs
            nlogC = vmax.*(1+(vmin./vmax).^alpha).^(1./alpha);
            y = (alpha - 1 + nlogC) ...
                   .* exp(-nlogC + sum((alpha-1).*log(v) + v, 2) + (1-2*alpha).*log(nlogC));
        end
    end
end

% Give values outside of unit hypercube zero density
if anyOutOfRange
    y(outOfRange & ~isnan(y)) = 0; % don't overwrite NaNs
end
