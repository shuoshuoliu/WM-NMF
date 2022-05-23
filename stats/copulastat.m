function r = copulastat(family,varargin)
%COPULASTAT Rank correlation for a copula.
%   R = COPULASTAT('Gaussian',RHO) returns the Kendall's rank correlation R
%   that corresponds to a Gaussian copula having linear correlation parameters
%   RHO.  If RHO is a scalar correlation coefficient, R is a scalar
%   correlation coefficient corresponding to a bivariate copula.  If RHO is a
%   P-by-P correlation matrix, R is a P-by-P correlation matrix.
%
%   R = COPULASTAT('t',RHO,NU) returns the Kendall's rank correlation R that
%   corresponds to a t copula having linear correlation parameters RHO and
%   degrees of freedom NU.  If RHO is a scalar correlation coefficient, R is a
%   scalar correlation coefficient corresponding to a bivariate copula.  If
%   RHO is a P-by-P correlation matrix, R is a P-by-P correlation matrix.
%   
%   R = COPULASTAT(FAMILY,ALPHA) returns the Kendall's rank correlation R that
%   corresponds to a bivariate Archimedean copula with scalar parameter ALPHA.
%   FAMILY is one of 'Clayton', 'Frank', or 'Gumbel'.
%
%   R = COPULASTAT(...,'type',TYPE) returns the specified type of rank
%   correlation.  TYPE is 'Kendall' to compute Kendall's tau, or 'Spearman' to
%   compute Spearman's rho.
%
%   COPULASTAT uses an approximation to Spearman's rank correlation for copula
%   families when no analytic formula exists.  The approximation is based on a
%   smooth fit to values computed at discrete values of the copula parameter(s).
%   For a t copula, the approximation is accurate for degrees of freedom larger
%   than 0.05.
%
%   Example:
%      % Get the theoretical rank correlation coefficient for a bivariate
%      % Gaussian copula having linear correlation parameter -0.7071
%      rho = -.7071
%      tau = copulastat('gaussian',rho)
%
%      % Generate dependent beta random values using that copula
%      u = copularnd('gaussian',rho,100);
%      b = betainv(u,2,2);
%
%      % Verify that the sample has a rank correlation approximately
%      % equal to tau
%      tau_sample = corr(b,'type','k')
%
%   See also COPULACDF, COPULAPDF, COPULARND, COPULAPARAM.

%   Copyright 2005-2018 The MathWorks, Inc.


if nargin > 0
    family = convertStringsToChars(family);
end

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 2
    error(message('stats:copulastat:WrongNumberOfInputs'));
elseif nargin > 2 && strncmpi(varargin{end-1},'type',max(1,numel(varargin{end-1})))
    type = varargin{end};
    types = {'kendall' 'spearman'};
    type = internal.stats.getParamVal(type,types,'''Type''');
    varargin = varargin(1:end-2);
else
    type = 'kendall';
end

families = {'gaussian','t','clayton','frank','gumbel'};
family = internal.stats.getParamVal(family,families,'FAMILY');

% Already stripped off 'type' args, should only be parameters left.
if isequal(family,'t')
    numParamArgs = 2; % provisionally
else
    numParamArgs = 1;
end
if length(varargin) > numParamArgs
    error(message('stats:copulastat:UnrecognizedInput'));
end

switch family
case {'gaussian' 't'}
    Rho = varargin{1};
    if strcmp(family,'t')
        if length(varargin) > 1
            if isnumeric(varargin{2})
                % optional nu was provided for the t copula
                nu = varargin{2};
                if ~(isscalar(nu) && (0 < nu))
                    error(message('stats:copulastat:BadDegreesOfFreedom'));
                end
            else
                error(message('stats:copulastat:UnrecognizedInput'));
            end
        elseif strcmp(type,'spearman')
            error(message('stats:copulastat:WrongNumberOfInputsMissingNu'));
        end
    end
    if isscalar(Rho)
        if ~(-1 <= Rho && Rho <= 1)
            error(message('stats:copulastat:BadScalarCorrelation'));
        end
    else
        if any(diag(Rho) ~= 1)
            err = 1;
        else
            [~,err] = cholcov(Rho);
        end
        if err ~= 0
            error(message('stats:copulastat:BadCorrelationMatrix'));
        end
    end
    switch type
    case 'kendall'
        r = 2.*asin(Rho)./pi;
    case 'spearman'
        if strcmp(family,'t')
            betaShape = internal.stats.tCopulaSpearmanUtil(nu);
            Rho = 2*betainc((Rho+1)/2,betaShape,betaShape) - 1;
        end
        r = 6.*asin(Rho./2)./pi;
    end
    perfectCorr = (abs(Rho) == 1);
    r(perfectCorr) = Rho(perfectCorr);
        
case {'clayton' 'frank' 'gumbel'}
    alpha = varargin{1};
    if ~isscalar(alpha)
        error(message('stats:copulastat:BadArchimedeanParameter'));
    end
    switch family
    case 'clayton'
        if alpha < 0
            error(message('stats:copulastat:BadClaytonParameter'));
        end
        switch type
        case 'kendall'
            r = alpha ./ (2 + alpha);
        case 'spearman'
            % A quintic in terms of alpha/(2+alpha), forced through y=0 and
            % y=1 at x=0 (alpha=0) and x=1 (alpha=Inf). This is a smooth fit
            % to sample rank correlations computed from MC simulations.
            a = -0.1002; b = 0.1533; c = -0.5024; d = -0.05629;
            p = [a b c d -(a+b+c+d-1) 0];
            r = polyval(p, alpha./(2+alpha));
        end
    case 'frank'
        if alpha == 0
            r = 0;
        else
            switch type
            case 'kendall'
                r = 1 + 4 .* (debye(alpha,1)-1) ./ alpha;
            case 'spearman'
                r = 1 + 12 .* (debye(alpha,2) - debye(alpha,1)) ./ alpha;
            end
        end
    case 'gumbel'
        if alpha < 1
            error(message('stats:copulastat:BadGumbelParameter'));
        end
        switch type
        case 'kendall'
            r = 1 - 1./alpha;
        case 'spearman'
            % A quintic in terms of 1/alpha, forced through y=1 and y=0 at x=0
            % (alpha=1) and x=1 (alpha=Inf). This is a smooth fit to sample
            % rank correlations computed from MC simulations.
            a = -.2015; b = .4208; c = .2429; d = -1.453;
            p = [a b c d -(a+b+c+d+1) 1];
            r = polyval(p, 1./alpha);
        end
    end
end

