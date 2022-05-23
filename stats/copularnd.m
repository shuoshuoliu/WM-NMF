function u = copularnd(family,varargin)
%COPULARND Random vectors from a copula.
%   U = COPULARND('Gaussian',RHO,N) returns N random vectors generated from a
%   Gaussian copula with linear correlation parameters RHO.  If RHO is a
%   P-by-P correlation matrix,  U is an N-by-P matrix.  If RHO is a scalar
%   correlation coefficient, COPULARND generates U from a bivariate Gaussian
%   copula.  Each column of U is a sample from a Uniform(0,1) marginal
%   distribution.
%
%   U = COPULARND('t',RHO,NU,N) returns N random vectors generated from a t
%   copula with linear correlation parameters RHO and degrees of freedom NU.
%   If RHO is a P-by-P correlation matrix, U is an N-by-P matrix. If RHO is a
%   scalar correlation coefficient, COPULARND generates U from a bivariate t
%   copula.  Each column of U is a sample from a Uniform(0,1) marginal
%   distribution.
%
%   U = COPULARND(FAMILY,ALPHA,N) returns N random vectors generated from the
%   bivariate Archimedean copula determined by FAMILY, with scalar parameter
%   ALPHA.  FAMILY is 'Clayton', 'Frank', or 'Gumbel'.  U is an N-by-2 matrix.
%   Each column of U is a sample from a Uniform(0,1) marginal distribution.
%
%   Example:
%      % Determine the linear correlation parameter corresponding to a
%      % bivariate Gaussian copula having a rank correlation of -0.5
%      tau = -0.5
%      rho = copulaparam('gaussian',tau)
%
%      % Generate dependent beta random values using that copula
%      u = copularnd('gaussian',rho,100)
%      b = betainv(u,2,2)
%
%      % Verify that the sample has a rank correlation approximately
%      % equal to tau
%      tau_sample = corr(b,'type','kendall')
%
%   See also COPULACDF, COPULAPDF, COPULASTAT, COPULAPARAM.

%   Copyright 2005-2020 The MathWorks, Inc.


if nargin > 0
    family = convertStringsToChars(family);
end

if nargin < 3
    error(message('stats:copularnd:WrongNumberOfInputs'));
end

families = {'gaussian','t','clayton','frank','gumbel'};
family = internal.stats.getParamVal(family,families,'FAMILY');

switch family

% Elliptical copulas
%
% Random vectors from these copulas can be generated by creating random
% vectors from the standard multivariate distribution, then transforming them
% to uniform marginals using the normal (or t) CDF.
case 'gaussian'
    Rho = varargin{1};
    n = varargin{2};
    d = size(Rho,1);
    if isscalar(Rho)
        if ~(-1 < Rho && Rho < 1)
            error(message('stats:copularnd:BadScalarCorrelation'));
        end
        Rho = [1 Rho; Rho 1];
        d = 2;
    elseif any(diag(Rho) ~= 1)
        error(message('stats:copularnd:BadCorrelationMatrix'));
    end
    % MVNRND will check that Rho is square, symmetric, and positive semi-definite.
    u = normcdf(mvnrnd(zeros(1,d),Rho,n));

case 't'
    if nargin < 4
        error(message('stats:copulapdf:MissingDF'));
    end
    Rho = varargin{1};
    nu = varargin{2};
    n = varargin{3};
    d = size(Rho,1);
    if isscalar(Rho)
        if ~(-1 < Rho && Rho < 1)
            error(message('stats:copularnd:BadScalarCorrelation'));
        end
        Rho = [1 Rho; Rho 1];
        d = 2;
    elseif any(diag(Rho) ~= 1)
        error(message('stats:copularnd:BadCorrelationMatrix'));
    end
    if ~(isscalar(nu) && (0 < nu))
        error(message('stats:copularnd:BadDegreesOfFreedom'));
    end
    % MVTRND will check that Rho is square, symmetric, and positive semi-definite.
    u = tcdf(mvtrnd(Rho,nu,n),nu);

% one-parameter Archimedean copulas
%
% Random pairs from these copulae can be generated sequentially: first
% generate u1 as a uniform r.v.  Then generate u2 from the conditional
% distribution F(u2 | u1; alpha) by generating uniform random values, then
% inverting the conditional CDF.
case {'clayton' 'frank' 'gumbel'}
    alpha = varargin{1};
    n = varargin{2};
    if ~isscalar(alpha)
        error(message('stats:copularnd:BadArchimedeanParameter'));
    end
    
    switch family
    case 'clayton' % a.k.a. Cook-Johnson
        if alpha < 0
            error(message('stats:copularnd:BadClaytonParameter'));
        end
        u1 = rand(n,1);
        % The inverse conditional CDF has a closed form for this copula.
        p = rand(n,1);
        if alpha < sqrt(eps)
            u2 = p;
        else
            u2 = u1.*(p.^(-alpha./(1+alpha)) - 1 + u1.^alpha).^(-1./alpha);
        end
        u = [u1 u2];

    case 'frank'
        u1 = rand(n,1);
        % The inverse conditional CDF has a closed form for this copula.
        p = rand(n,1);
        if abs(alpha) > log(realmax)
            u2 = (alpha < 0) + sign(alpha).*u1; % u1 or 1-u1
        elseif abs(alpha) > sqrt(eps)
            u2 = -log((exp(-alpha.*u1).*(1-p)./p + exp(-alpha))./(1 + exp(-alpha.*u1).*(1-p)./p))./alpha;
%             u2 = -log(1 + (1-exp(alpha))./(exp(alpha) + exp(alpha.*(1-u1)).*(1-p)./p))./alpha;
        else
            u2 = p;
        end
        u = [u1 u2];

    case 'gumbel' % a.k.a. Gumbel-Hougaard
        if alpha < 1
            error(message('stats:copularnd:BadGumbelParameter'));
        end
        if alpha < 1 + sqrt(eps)
            u = rand(n,2);
        else
            % This uses the Marshal-Olkin method, not inversion
            % Generate gamma as Stable(1/alpha,1), c.f. Devroye, Thm. IV.6.7
            u = (rand(n,1) - .5) .* pi; % unifrnd(-pi/2,pi/2,n,1)
            u2 = u + pi/2;
            e = -log(rand(n,1)); % exprnd(1,n,1)
            t = cos(u - u2./alpha) ./ e;
            gamma = (sin(u2./alpha)./t).^(1./alpha) .* t./cos(u);
            % Frees&Valdez, eqn 3.5
            s = (-log(rand(n,2))).^(1./alpha) ./ repmat(gamma,1,2);
            u = exp(-s);
        end
    end
end
