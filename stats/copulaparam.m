function param = copulaparam(family,varargin)
%COPULAPARAM Copula parameters as a function of rank correlation.
%   RHO = COPULAPARAM('Gaussian',R) returns the linear correlation parameters
%   RHO corresponding to a Gaussian copula having Kendall's rank correlation
%   R.  If R is a scalar correlation coefficient, RHO is a scalar correlation
%   coefficient corresponding to a bivariate copula.  If R is a P-by-P
%   correlation matrix, RHO is a P-by-P correlation matrix.
%
%   RHO = COPULAPARAM('t',R,NU) returns the linear correlation parameters RHO
%   corresponding to a t copula having Kendall's rank correlation R and
%   degrees of freedom NU.  If R is a scalar correlation coefficient, RHO is
%   a scalar correlation coefficient corresponding to a bivariate copula.  If
%   R is a P-by-P correlation matrix, RHO is a P-by-P correlation matrix.
%   
%   ALPHA = COPULAPARAM(FAMILY,R) returns the copula parameter ALPHA
%   corresponding to a bivariate Archimedean copula having Kendall's rank
%   correlation R.  R is a scalar.  FAMILY is one of 'Clayton', 'Frank',
%   or 'Gumbel'.
%
%   [...] = COPULAPARAM(...,'type',TYPE) assumes R is the specified type of
%   rank correlation.  TYPE is 'Kendall' for Kendall's tau, or 'Spearman' for
%   Spearman's rho.
%
%   COPULAPARAM uses an approximation to Spearman's rank correlation for copula
%   families where no analytic formula exists. The approximation is based on a
%   smooth fit to values computed at discrete values of the copula parameter(s).
%   For a t copula, the approximation is accurate for degrees of freedom larger
%   than 0.05.
%
%   Example:
%      % Get the linear correlation coefficient corresponding to a bivariate
%      % Gaussian copula having a rank correlation of -0.5
%      tau = -0.5
%      rho = copulaparam('gaussian',tau)
%
%      % Generate dependent beta random values using that copula
%      u = copularnd('gaussian',rho,100);
%      b = betainv(u,2,2);
%
%      % Verify that the sample has a rank correlation approximately
%      % equal to tau
%      tau_sample = corr(b,'type','k')
%
%   See also COPULACDF, COPULAPDF, COPULARND, COPULASTAT.

%   Copyright 2005-2019 The MathWorks, Inc.


if nargin > 0
    family = convertStringsToChars(family);
end

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 2
    error(message('stats:copulaparam:WrongNumberOfInputs'));
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
    error(message('stats:copulaparam:UnrecognizedInput'));
end

switch family
case {'gaussian' 't'}
    R = varargin{1};
    if strcmp(family,'t')
        if length(varargin) > 1
            if isnumeric(varargin{2})
                % optional nu was provided for the t copula
                nu = varargin{2};
                if ~(isscalar(nu) && (0 < nu))
                    error(message('stats:copulaparam:BadDegreesOfFreedom'));
                end
            else
                error(message('stats:copulaparam:UnrecognizedInput'));
            end
        elseif strcmp(type,'spearman')
            error(message('stats:copulaparam:WrongNumberOfInputsMissingNu'));
        end
    end
    if isscalar(R)
        if ~(-1 < R && R < 1)
            error(message('stats:copulaparam:BadScalarCorrelation'));
        end
    else
        err = 1;
        if all(diag(R) == 1)
            [~,err] = cholcov(R);
        end
        if err ~= 0
            error(message('stats:copulaparam:BadCorrelationMatrix'));
        end
    end
    switch type
    case 'kendall'
        param = sin(R.*pi./2);
    case 'spearman'
        param = 2.*sin(R.*pi./6);
        if isequal(family,'t')
            betaShape = internal.stats.tCopulaSpearmanUtil(nu);
            param = 2*betaincinv((param+1)/2,betaShape,betaShape) - 1;
        end
    end
    perfectCorr = (abs(R) == 1);
    param(perfectCorr) = R(perfectCorr);
    
case {'clayton' 'frank' 'gumbel'}
    r = varargin{1};
    if ~isscalar(r) || ~(-1 <= r && r <= 1)
        error(message('stats:copulaparam:BadCorrelation'));
    end
    switch family
    case 'clayton'
        if r < 0
            error(message('stats:copulaparam:NegativeCorrelation', family));
        end
        switch type
        case 'kendall'
            param = 2*r ./ (1-r);
        case 'spearman'
            % A quintic in terms of alpha/(2+alpha), forced through y=0 and
            % y=1 at x=0 (alpha=0) and x=1 (alpha=Inf). This is a smooth fit
            % to sample rank correlations computed from MC simulations.
            a = -0.1002; b = 0.1533; c = -0.5024; d = -0.05629;
            p = [a b c d -(a+b+c+d-1) 0];
            t = fzero(@(t) polyval(p,t)-r, [0 1]);
            param = 2*t ./ (1-t);
        end
    case 'frank'
        if r == 0
            param = 0;
        elseif abs(r) < 1
            % There's no closed form for alpha in terms of tau, so alpha has
            % to be determined numerically.
            switch type
            case 'kendall'
                param = fzero(@frankRootFunKendall,sign(r),[],r);
            case 'spearman'
                param = fzero(@frankRootFunSpearman,sign(r),[],r);
            end
        else
            param = sign(r).*Inf;
        end
    case 'gumbel'
        if r < 0
            error(message('stats:copulaparam:NegativeCorrelation', family));
        end
        switch type
        case 'kendall'
            param = 1 ./ (1-r);
        case 'spearman'
            % A quintic in terms of 1/alpha, forced through y=1 and y=0 at x=0
            % (alpha=1) and x=1 (alpha=Inf). This is a smooth fit to sample
            % rank correlations computed from MC simulations.
            a = -.2015; b = .4208; c = .2429; d = -1.453;
            p = [a b c d -(a+b+c+d+1) 1];
            t = fzero(@(t) polyval(p,t)-r, [0 1]);
            param = 1 ./ t;
        end
    end
end


function err = frankRootFunKendall(alpha,target)
if abs(alpha) < sqrt(realmin)
    r = 0;
else
    r = 1 + 4 .* (debye(alpha,1)-1) ./ alpha;
end
err = r - target;


function err = frankRootFunSpearman(alpha,target)
if abs(alpha) < sqrt(realmin)
    r = 0;
else
    r = 1 + 12 .* (debye(alpha,2)-debye(alpha,1)) ./ alpha;
end
err = r - target;
