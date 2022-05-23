classdef BurrDistribution < prob.ToolboxFittableParametricDistribution
%BurrDistribution Burr probability distribution.
%    An object of the BurrDistribution class represents a Burr
%    probability distribution with specific values for the ALPHA, C,
%    and K parameters. This distribution object can be created directly
%    using the MAKEDIST function or fit to data using the FITDIST function.
%
%    BurrDistribution methods:
%       cdf                   - Cumulative distribution function
%       icdf                  - Inverse cumulative distribution function
%       iqr                   - Interquartile range
%       mean                  - Mean
%       median                - Median
%       negloglik             - Negative log likelihood function
%       paramci               - Confidence intervals for parameters
%       pdf                   - Probability density function
%       proflik               - Profile likelihood function
%       random                - Random number generation
%       std                   - Standard deviation
%       truncate              - Truncation distribution to an interval
%       var                   - Variance
%
%    BurrDistribution properties:    
%       DistributionName      - Name of the distribution
%       alpha                 - Value of the alpha parameter (scale)
%       c                     - Value of the c parameter (first shape parameter)
%       k                     - Value of the k parameter (second shape parameter)
%       NumParameters         - Number of parameters
%       ParameterNames        - Names of parameters
%       ParameterDescription  - Descriptions of parameters
%       ParameterValues       - Vector of values of parameters
%       Truncation            - Two-element vector indicating truncation limits
%       IsTruncated           - Boolean flag indicating if distribution is truncated
%       ParameterCovariance   - Covariance matrix of estimated parameters
%       ParameterIsFixed      - Two-element boolean vector indicating fixed parameters
%       InputData             - Structure containing data used to fit the distribution
%
%    See also fitdist, makedist.

%    Copyright 2012-2018 The MathWorks, Inc.

    properties(Dependent=true)
%ALPHA ALPHA (scale) parameter
%    The ALPHA property represents the scale parameter of the
%    Burr distribution.
%
%    See also C, K.
        alpha

%C C (shape) parameter
%    The C property represents the first of two shape parameters of the
%    Burr distribution.
%
%    See also ALPHA, K.
        c

%K K (shape) parameter
%    The K property represents the second of two shape parameters of the
%    Burr distribution.
%
%    See also ALPHA, C.
        k
    end
    properties(GetAccess='public',Constant=true)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameBurr'));

%NumParameter Number of parameters.
%    NumParameters is the number of parameters in the distribution.
%
%    See also ParameterValues.
        NumParameters = 3;

%ParameterNames Parameter names.
%    ParameterNames is a cell array of strings containing the names of the
%    parameters of the probability distribution.
%
%    See also ParameterValues, ParameterDescription.
        ParameterNames = {'alpha' 'c' 'k'};

%ParameterDescription Parameter description.
%    ParameterNames is a cell array of strings containing short
%    descriptions of the parameters of the probability distribution.
%
%    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionScale')) ...
                                getString(message('stats:probdists:ParameterDescriptionShape1')) ...
                                getString(message('stats:probdists:ParameterDescriptionShape2'))};
    end
    properties(GetAccess='public',SetAccess='protected')
%ParameterValues Parameter values.
%    ParameterVales is a vector containing the values of the parameters of
%    the probability distribution.
%
%    See also ALPHA, C, K.
        ParameterValues
    end
    methods(Hidden)
        function pd = BurrDistribution(alpha,c,k)
            if nargin==0
                alpha = 1;
                c = 1;
                k = 1;
            end
            checkargs(alpha,c,k)
            
            pd.ParameterValues = [alpha c k];
            pd.ParameterIsFixed = [true true true];
            pd.ParameterCovariance = zeros(pd.NumParameters);
        end
    end
    methods
        function m = mean(this)
            requireScalar(this)
            if this.IsTruncated
                m = truncatedMoment(this,1);
                return
            end
            A = this.alpha;
            C = this.c;
            K = this.k;
            valid = A > 0 & C > 0 & K > 0;
            
            % ck used in determining need to calculate moments
            ck = C .* K;
            
            % Preallocate space for mean
            m = nan(size(A));
            hasmean = valid & ck > 1;
            % Infinite means
            m(valid & ck <= 1) = Inf;
            
            % Precalculate as many quantities as possible
            ci = 1./C;
            ac = A .* ci;
            k_ci = K - ci;
            
            % Greatly simplified moments for k == 1
            % using the identity that B(x, 1 - x) = Gamma(x)Gamma(1-x) = pi/sin(pi x)
            simple = (K == 1);
            simplemean = simple & hasmean;
            if any(simplemean(:))
                m(simplemean) = ac(simplemean).*pi./sin(pi.*ci(simplemean));
            end
            
            normal = ~simple;
            normalmean = normal & hasmean;
            if any(normalmean(:))
                m(normalmean) = ac(normalmean) .* beta(ci(normalmean), k_ci(normalmean));
            end
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            A = this.alpha;
            C = this.c;
            K = this.k;
            valid = A > 0 & C > 0 & K > 0;
            
            % ck used in determining need to calculate moments
            ck = C .* K;
            
            % Preallocate space for mean
            m = nan(size(A));
            hasmean = valid & ck > 1;
            % Infinite means
            m(valid & ck <= 1) = Inf;
            
            % Precalculate as many quantities as possible
            ci = 1./C;
            ac = A .* ci;
            k_ci = K - ci;
            
            % Greatly simplified moments for k == 1
            % using the identity that B(x, 1 - x) = Gamma(x)Gamma(1-x) = pi/sin(pi x)
            simple = (K == 1);
            simplemean = simple & hasmean;
            if any(simplemean(:))
                m(simplemean) = ac(simplemean).*pi./sin(pi.*ci(simplemean));
            end
            
            normal = ~simple;
            normalmean = normal & hasmean;
            if any(normalmean(:))
                m(normalmean) = ac(normalmean) .* beta(ci(normalmean), k_ci(normalmean));
            end
            
            v = nan(size(A));
            
            hasvar  = valid & ck > 2;
            % Infinite variances
            v(valid & ck <= 2) = Inf;
            
            simplevar = simple & hasvar;
            normalvar = normal & hasvar;
            
            
            if any(simplevar(:))
                v(simplevar) = 2 * A(simplevar) .*ac(simplevar).*pi./sin(2.*pi.*ci(simplevar)) - m(simplevar).^2;
            end
            if any(normalvar(:))
                v(normalvar) = 2 * A(normalvar) .* ac(normalvar) .* beta(ci(normalvar) .* 2, k_ci(normalvar) - ci(normalvar)) - m(normalvar).^2;
            end
        end
        function this = set.alpha(this,alpha)
            checkargs(alpha,this.c,this.k)
            this.ParameterValues(1) = alpha;
            this = invalidateFit(this);
        end
        function this = set.c(this,c)
            checkargs(this.alpha,c,this.k)
            this.ParameterValues(2) = c;
            this = invalidateFit(this);
        end
        function this = set.k(this,k)
            checkargs(this.alpha,this.c,k)
            this.ParameterValues(3) = k;
            this = invalidateFit(this);
        end
        function alpha = get.alpha(this)
            alpha = this.ParameterValues(1);
        end
        function c = get.c(this)
            c = this.ParameterValues(2);
        end
        function k = get.k(this)
            k = this.ParameterValues(3);
        end
    end
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the normal distribution to data.
%    Fitting requires the Statistics and Machine Learning Toolbox. You should call the FITDIST
%    function instead of calling this method directly.
%
%    See also FITDIST.
            [x,cens,freq,opt] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            p = burrfit(x,0.05,cens,freq,opt);
            [nll,cov] = burrlike(p,x,cens,freq);
            pd = prob.BurrDistribution.makeFitted(p,nll,cov,x,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = burrlike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = burrcdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = burrpdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = burrinv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = burrrnd(varargin{:});
        end
        function pd = makeFitted(p,nll,cov,x,cens,freq)
            pd = prob.BurrDistribution(p(1),p(2),p(3));
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.BurrDistribution');
            info.name = prob.BurrDistribution.DistributionName;
            info.code = 'burr';
            info.censoring = true;
            info.support = [0 Inf];
            info.optimopts = true;
            info.logci = [true true true];
            info.plim = [0 0 0;Inf Inf Inf];
        end
    end
end % classdef

function checkargs(alpha,c,k)
if ~(isscalar(alpha) && isnumeric(alpha) && isreal(alpha) && isfinite(alpha) && alpha>0)
    error(message('stats:probdists:PositiveParameter','ALPHA'))
end
if ~(isscalar(c) && isnumeric(c) && isreal(c) && isfinite(c) && c>0)
    error(message('stats:probdists:PositiveParameter','C'))
end
if ~(isscalar(k) && isnumeric(k) && isreal(k) && isfinite(k) && k>0)
    error(message('stats:probdists:PositiveParameter','K'))
end
end


function f = burrpdf(x, alpha, c, k)
%BURRPDF Burr probability density function (pdf).

[errorcode, x, alpha, c, k] = distchck(4, x, alpha, c, k);

if errorcode > 0
    error(message('stats:addburr:InputSizeMismatch'));
end

f = zeros(size(x));

alpha(alpha <= 0) = NaN;
c(c <= 0) = NaN;
k(k <= 0) = NaN;

% Specific values at 0
smallc = c < 1;
c1 = c == 1;
x0 = x == 0;

if (any(x0(:)))
    % For c < 1 it explodes
    f(smallc & x0) = Inf;

    % For c = 1 it converges
    c1x0 = c1 & x0;
    f(c1x0) = k(c1x0)./alpha(c1x0);

    % For c > 1, it is zero (no need to do anything)
end

normal = x > 0;
xa   = x(normal) ./alpha(normal);
xac  = xa.^c(normal);
xac1 = xac./xa;

% No formula simplification for k = 1
f(normal) = (k(normal).*c(normal)./alpha(normal)) .* xac1 ./ (xac + 1).^(k(normal) + 1);

% Find instances of overflow
overflow = isinf(xac);
if any(overflow(:))
    f(overflow) = 0;
end

% Pass input NaN values through to output.
f(isnan(x)) = NaN;

end

function F = burrcdf(x, alpha, c, k, tail)
%BURRCDF Burr cumulative distribution function (cdf).

[errorcode, x, alpha, c, k] = distchck(4, x, alpha, c, k);

if errorcode > 0
    error(message('stats:addburr:InputSizeMismatch'));
end

% Will force negative #s to evaluate to zero
x(x < 0) = 0; 
alpha(alpha <= 0) = NaN;
c(c <= 0) = NaN;
k(k <= 0) = NaN;

% Precalculate values needed in all formulas
xac  = (x./alpha).^c;

% If xac contains Inf, then use Pareto distribution
infxac = isinf(xac);

if nargin > 4 
    % Upper tail calculation
    if ~strcmpi(tail, 'upper')
        error(message('stats:cdf:UpperTailProblem'));
    else
        F(~infxac) = (1 + xac(~infxac)).^-k(~infxac);
        F(infxac) = (x(infxac)./alpha(infxac)).^-(c(infxac).*k(infxac));
    end
else
    F = zeros(size(x));

    % Much more simple, reliable formula, for k == 1
    % (reworked log-logistic distribution)
    simple = (k == 1);
    if any(simple(:))
        F(simple) = 1./(1 + 1./xac(simple));
    end

    % Normal evaluation
    normal = ~simple;
    if any(normal(:))
        xac1_nk = (1 + xac).^-k;
        F(normal&~infxac) = 1 - xac1_nk(normal&~infxac);
        F(normal&infxac) = 1 - (x(normal&infxac)./alpha(normal&infxac)).^-(c(normal&infxac).*k(normal&infxac));
        F_small = zeros(size(F));

        %Find small values of (x/alpha)^c, and use second degree Taylor expansion for any CDF smaller than 10e-11.
        %For CDF values within [10e-11, 10e-10], use a linear combination of the approxiamtion and the original formula.

        %Find small values of (x/alpha)^c, control the truncation error to be bounded by epsilon^3.
        epsilon = 10e-5;
        min_k1_ep = epsilon.*min(1./k,1);
        small = (xac <= min_k1_ep);
        if any(small(:))
            %For small (x/alpha)^c, use Taylor approxiamtion
            F_small(small) = k(small).*xac(small)-1/2*k(small).*(k(small)+1).*xac(small).^2;

            %Linear combination of the approximation and the original formula for CDF values within [10e-11, 10e-10]
            blended = normal & (xac1_nk>=(1-1e-10)) & (xac1_nk<=(1-1e-11)) & small;
            blending_weights = (xac1_nk(blended) - (1-1e-10))./(1e-10-1e-11);
            F(blended) = blending_weights .* F(blended) + (1 - blending_weights) .* F_small(blended);

            %For small CDF of values below 10e-11, use approximation only
            F(small & ~blended) = F_small(small & ~blended);
        end
    end
end
end

function x = burrinv(F, alpha, c, k)
%BURRINV Inverse of the Burr cumulative distribution function (cdf).

[errorcode, F, alpha, c, k] = distchck(4, F, alpha, c, k);

if errorcode > 0
    error(message('stats:addburr:InputSizeMismatch'));
end

F(F < 0 | F > 1) = NaN;
alpha(alpha <= 0) = NaN;
c(c <= 0) = NaN;
k(k <= 0) = NaN;

% Common factors
ci = c.^-1;

x = zeros(size(F));

% Greatly simplified formula for k = 1 (log-logistid dist.)
simple = (k == 1);
if any(simple(:))
    x(simple) = alpha(simple) .* (1./F(simple) - 1).^-ci(simple);
end


% General formula: first calculate x values for all Fs
normal = ~simple;
if any(normal(:))
    ki = 1./k;
    ci = 1./c;
    F1_nki = (1 - F).^(-ki);

    x(normal) = alpha(normal) .* (F1_nki(normal) - 1) .^ ci(normal);


    % Crossover range: blend the two formulas during this range
    % Use rescaled x-value (x/alpha)^c = F1_nki to avoid effects of alpha and c
    % i.e., use the value that will have 1 subtracted from it
    % to turn on the approximation
    blend_range = [1.000001324210863 1.000001523518778];
    small       = normal & F1_nki < blend_range(2); % All small values of x

    if any(small(:))
        blended     = normal & F1_nki > blend_range(1) & small; % Values over which the blend will happen
        % Linear fade from one to the other
        blending_weights = (F1_nki(blended) - blend_range(1))./(blend_range(2) - blend_range(1)); 

        % Calculate all the instances of the lower-tail approximation
        % that will be necessary
        x_small = zeros(size(x));
        x_small(small) =  alpha(small) .* (F(small).*ki(small).*(1 + F(small)/2.*(1 - ki(small))) .* F1_nki(small)) .^ (ci(small));

        % Use a weighted blend of the lower tail and the normal formula over the blend range
        x(blended) = blending_weights .* x(blended) + (1 - blending_weights) .* x_small(blended);

        % For anything smaller than the blend range, just use the lower tail approximation
        x(small & ~blended) = x_small(small & ~blended);
    end
end
end

function r = burrrnd(alpha, c, k, varargin)
%BURRRND Random arrays from the Burr distribution.

alpha(alpha <= 0) = NaN;
c(c <= 0) = NaN;
k(k <= 0) = NaN;

[err, sizeOut] = internal.stats.statsizechk(3, alpha, c, k, varargin{:});
if err > 0
    error(message('stats:addburr:InputSizeMismatch'));
end

% Use the inverse CDF, sampling from the uniform distribution
F = rand(sizeOut);

% (see burrinv above: don't care about accuracy for lower tail here,
% since the numbers are random anyway)
ki = 1./k;
ci = 1./c;
F1_nki = (1 - F).^(-ki);
r = alpha .* (F1_nki - 1) .^ ci;
end

function [nL,acov] = burrlike(parms,data,cens,freq)
%BURRLIKE Negative log-likelihood for the Burr distribution.
% Standardize x data
x = data(:);
if nargin < 3 || isempty(cens), cens = false(size(x)); else cens = logical(cens(:)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); else freq = freq(:); end

alpha = parms(1);
c = parms(2);
k = parms(3);

xa = x ./ alpha;
lxa = log(xa);

% Careful Logarithm: log(1 + x) = x for x < eps
l1_xac = log1p(xa .^ c);

% This has a tendency to overflow, but for large x,
% 1 + x $\approx$ x, and so log1p(xa^c) $\approx$ c*log(xa)
overflow = isinf(l1_xac);
l1_xac(overflow) = c.*log(xa(overflow));

% Normalize the censored flags, just in case they weren't
% logicals to begin with
uncens = ~cens;

% Frequently calculated sum
n_uncens = sum(freq(uncens));

% Uncensored data log-likelihoods
% Preallocate to prevent reallocation in case last elements are censored
L = zeros(size(x));
L(uncens) = (c - 1)*lxa(uncens) - (k + 1)*l1_xac(uncens);

% Censored data log-likelihoods
L(cens) = -k * l1_xac(cens);

% Add in precalcualted constants, sum and sign flip for the final value
nL = n_uncens * log(alpha/k/c) - sum(freq .* L);

if nargout > 1
    acov = burrcov(parms, data, cens, freq);
end
end

% ==== Burr fitting functions ====

function [phat,pci] = burrfit(x,alpha,cens,freq,opts)
%BURRFIT Parameter estimates and confidence intervals for Burr data.
%
% Due to the flexiblity of the Burr distribution, and its many limiting cases,
% this fit must follow a complicated procedure, outlined in the literature, and
% sketched below.
% 
% 1. The data must be fit to a Pareto (not generalized Pareto) distribution.
% 2. The data must be fit to a Weibull distribution 
% 3. The Weibull fit is used to compute a discriminator (delta), which
%    indicates if the Weibull is a better fit for the data than the Burr (if < 0).
% 4. If Delta > 0, calculate a Burr fit, with a small lower boundary for k. If k goes to that
%    lower boundary, throw out the fit. If it converges to a normal value of k,
%    then save that fit.
% 5. If Delta < 0, then do not calcualte a Burr fit, and save the Weibull fit instead
% 6. If a fit was saved, compare the likelihood of that fit to the likelihood of the
%    Pareto fit. Pick whichever one is greater.
% 7. If no fit was saved in 4&5, then just return the Pareto fit.
% 

%   References:
%      [1] Shao, Q. (2004) "Notes on maximum likelihood estimation for the
%          three-parameter Burr XII Distribution". Computational Statistics 
%          and Data Analysis. v. 45, pp. 675 - 687

% Standardize x data
x = x(:);
if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = false(size(x)); else cens = logical(cens(:)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); else freq = freq(:); end
if nargin < 5, opts = []; end

if any(x(:) <= 0)
    error(message('stats:addburr:BadData'));
end

% Check for identical data
if ~isscalar(x) && max(abs(diff(x(:))) ./ x(2:end)) <= sqrt(eps)

    warning(message('stats:addburr:NoDistinctX'))

    % Identical data
    phat = [x(1) Inf sum(~cens.*freq)/sum(freq)/log(2)];

    if nargout > 1
        pci = [phat; phat];
    end

    return
end

% Create flags to signify which distribution is the best fit
best_wbl = false;
best_par = false;

% Fit a pareto distribution to the data first: this is needed 
% in the k->0 limit, which always has a local maximum in likelihood
[phat_par, nL_par] = paretofit(x, cens, freq);

% Fit a Weibull distribution to the data: this is needed in the
% k->inf limit, which may be a better fit for the data than any
% traditional Burr distribution
[phat_wbl] = wblfit(x, alpha, cens, freq);
nL_wbl = wbllike(phat_wbl, x, cens, freq);
theta_wbl = phat_wbl(1);
beta_wbl = phat_wbl(2);

% Calculate the discriminator
x_theta = x ./theta_wbl;
x_theta_beta = x_theta.^beta_wbl;
delta_wbl = sum(freq.*(0.5 * x_theta_beta.^2 - (~cens).*x_theta_beta));

if delta_wbl > 0 % We want to calculate the Burr distribution
    % Signify to mle that we have an analytical gradient, should it choose to use it
    % and set high tolerance
    opts = statset(opts, 'GradObj', 'on');

    % Starting values:
    % Solve equations using the median and upper quartile in order
    % to have a good estimate of alpha and c, for better convergence
    qs = prctile(reps(x_theta(:), freq), [50 75]);

    % If the median and upper quartile are too close together, spread them
    % out. Otherwise the starting parameter estimate will be impossible
    dq = sqrt(eps(qs(1))) * qs(1);
    if (qs(2) - qs(1)) < dq
        if any(x_theta > qs(2))
            % Try to move the upper quartile to the next highest point
            qs(2) = min(x_theta(x_theta > qs(2)));
        elseif any(x_theta < qs(2))
            % All data are smushed together on the upper end. Move the median down instead
            qs(1) = max(x_theta(x_theta < qs(1)));
        end
        % If there are no data smaller than the median and no data larger than the upper quartile,
        % and both have the same value, then the data must all be identical.
        % This case has been eliminated previously
    end

    % For the rescaled data, a great estimate of a0 is ~1, even if the Burr distribution
    % is a better fit, because theta_wbl should be quite close to alpha.
    % Assuming that, there is only a single nonlinear equation to solve for c0:
    % (qs(2)/qs(1)).^c - qs(1).^c - 2)
    % This is a reduced form of a full equation (qs(1) is really scaled by alpha0, which vanishes
    % as we assume it's unity. The equation will only have a solution as long as qs(2)/qs(1) > qs(1)
    %
    % Thus, if the equation does not have an apparent solution, we scale a0 up until it does,
    % letting a0 = qs(2), and giving c0 a closed form. Otherwise, the equation holds, and we
    % can minimize directly.
    %
    if qs(1) >= qs(2)/qs(1)
        a0 = qs(1);
        c0 = log(3)/log(qs(2)/qs(1));
    else
        % It's possible to find a solution, so we must do the search.

        % We must bound our search to prevent invalid answers. When c = 0, the constraint
        % will always be negative, so that is a safe lower bound. The upper bound
        % can be found by finding a c that will generate a really large number, say,
        % sqrt(realmax). Since, in the case where we're doing the search, it is already
        % clear that qs(2)/qs(1) > qs(1), this number is the best shot we have
        % for picking something guaranteed to be positive
        cmax  = log(realmax)/(2*log(qs(2)/qs(1)));

        a0 = 1;
        c0 = fzero(@(c) (qs(2)/qs(1)).^c - qs(1).^c - 2, ...
                   [0 cmax], optimset(opts, 'Display', 'off'));
    end

    % Don't warn about reaching the iteration limit unless we really
    % need to: if we do worse than the Pareto, then there's no reason
    % to mention it, since we're about to throw a bigger error anyway
    % and that error is the reason that the optimization did not converge
    % within the limit
    wquery = warning('query', 'stats:mle:IterLimit');  % preserve current warning state
    warning('off', 'stats:mle:IterLimit');             % suppress wng
    [oldw_msg, oldw_id] = lastwarn();                  % preserve last warning
    lastwarn('');                                      % clear last warning
    cleaner = onCleanup(@() warning(wquery));          % will restore original warning state

    % Use the reduced log-likelihood to search for alpha and c, then
    % calculate a final k-hat from the results
    % In the search, use rescaled x-data: more stable, according
    % to Watkins, 1999.
    phat = mle(x_theta, 'nloglf', @burr_nloglf_obj, ...
        'censoring', cens, 'frequency', freq,  ...
        'start', [a0 c0], 'lowerbound', [0 0],   ...
        'alpha', alpha, 'options', opts);

    % Load up the last warning to see if MLE reached the iter limit.
    % If it did, lastwarn will report the right ID, even though
    % the warning was flagged "off".
    [~,lastw_id] = lastwarn();  
    % Restore lastwarn so that mle warn, if any, is temporarily ignored.
    % We will issue the warning later if warranted.
    lastwarn(oldw_msg, oldw_id);                     

    % K has a closed form, that we can now obtain from the other params
    kh = exp(logkhat(x_theta, phat(1), phat(2), cens, freq));
    phat = [phat kh];

    % Since all we did was rescale the data, the only thing affected by it
    % would be the scale parameter, alpha = phat(1). So now we scale it back
    phat(1) = phat(1) * theta_wbl;

    % Now that we have optimum paramter values, we also need to know
    % the minimum value of the neg-log-likelihood for the original
    % data with the true parameters. Since k-hat is always automatically
    % calculated from the other two parameters, we can use the shortcut
    % form of the log-likelihood.
    nL = burr_nloglf_obj(phat(1:2), x, cens, freq);

    % Now we need to see if our Burr distribution case is actually valid and optimal
    % We do this by seeing if the value of k it provides is "tiny", i.e., if
    % it hits some eps0
    % Instead of calculating out an exact theoretical value for eps0, as
    % outlined in (Shao, 2004), we can just practially choose an eps0 
    % below which it will not make sense to have a reasonable k. A good
    % practical choice is 1e-6;
    %
    % To see if it's optimal, we still have to compare its likelihood
    % to the Pareto likelihood. If either the khat value is too small
    % or the likelihood is too low (rather, the neg-log-likelihood is
    % too high), then we choose the pareto instead.
    eps0 = 1e-6;
    if kh <= eps0 || nL > nL_par
        best_par = true;
    elseif nargout > 1 
        % We're actually returning the burr distribution, 
        % so now calculate the cov's if asked for

        % Confidence interval estimates for parameters
        acov = burrcov(phat, x, cens, freq);

        % Use logarithm of parameters to calculate confidence intervals
        % this works better for two reasons:
        % 1. not generating negative estimates, which are nonsensical
        % 2. the log of the parameters seems to be more normally distributed
        pci = statparamci(phat, acov, alpha, true(1, 3));
    end

    % If we really do need to warn the user, then warn
    if ~best_par && strcmp(lastw_id, 'stats:mle:IterLimit') && strcmp(wquery.state, 'on')
        warning(message('stats:mle:IterLimit'));
    end

else
    % No need to look for the Burr distribution, since the Weibull fit
    % will always be better. Now we just need to compare the Weibull fit to 
    % the Pareto fit.
    if nL_par < nL_wbl
        best_par = true;
    else
        best_wbl = true;
    end
end

% If the best fit is not a Burr distribution, bomb out with the appropriate error
if best_par
    % Matlab generalized Pareto is param'd differently
    loc = sprintf('%g', phat_par(1));
    scale = sprintf('%g', phat_par(1)/phat_par(2)); 
    shape = sprintf('%g', 1/phat_par(2)); % Matlab generalized Pareto is param'd differently
    error(message('stats:addburr:ParetoBetter', shape, scale, loc));
elseif best_wbl
    scale = sprintf('%g', phat_wbl(1));
    shape = sprintf('%g', phat_wbl(2));
    error(message('stats:addburr:WeibullBetter', scale, shape));
end
end

function [repd] = reps(x, freq)
% Repeat values of x according to freq
np = sum(freq);
repd = zeros(1, np, class(x));
rx = 1;
for ix = 1:length(x)
    fx = freq(ix);
    repd(rx:rx + fx -1) = x(ix);
    rx = rx + fx;
end
end

function [phat_par, nL_par] = paretofit(x, cens, freq)
% Fit a true pareto distribution to the data. GP fit won't do: it may be too
% general, giving an unfarily high likelihood, when all you really want to do
% is fit a strict pareto to work as a limiting case of the Burr.

% The pareto is formulated with a CDF F(x) = 1 - (x0/x)^lambda

% Fortunately, the strict pareto has well-formed MLE estimates:
uncens = ~cens;
x0_hat = x(uncens);
x0_hat = min(x0_hat(:));

% If all data is censored, the likelihood is maximised at zero by x0 >= max(x)
% and lambda = undefined
if all(cens)
    phat_par = [ max(x) NaN ];
    nL_par = 0;
    return
end

lx = log(x);

% Precalculate commonly used values
uncens_freq = freq.*uncens;
n_uncens = sum(uncens_freq);
sum_logs = sum(freq .* (lx - log(x0_hat)).* (x > x0_hat));

lambda_hat = n_uncens ./ sum_logs ;

phat_par = [x0_hat lambda_hat];

% Negative log-likelihood with MLE params
nL_par = lambda_hat .* sum_logs + sum(uncens_freq .* lx) - log(lambda_hat) * n_uncens;
end

function acov = burrcov(phat, x, cens, freq)
%BURRCOV Construct analytical, closed-form version of log-likelihood hessian
%        and invert it to obtain the covariance matrix.

alpha = phat(1);
c = phat(2);
k = phat(3);

% Normalize the censored flags, just in case they weren't
% logicals to begin with
uncens = ~cens;

% Preallocate memory in only one vector. Use it for all the elements
d2Lvec = zeros(size(x));

% Use a separate vector to construct the base of the formula
base = zeros(size(x));

% Precalculate commonly used quantities
xa = x ./ alpha;
xac = xa .^ c;
lxa = log(xa);
xac1 = (1 + xac);
xac1_2 = xac1.^2;
xaq_2 = (1 + 1./xac).^2;

n_uncens = sum(freq(uncens));

% Find instances of overflow
overflow = isinf(xac);
safe = ~overflow;

% Diagonal elements

% alpha-alpha
base(safe) = -((1 + c)./xac(safe) + 1)./xaq_2(safe);
base(overflow) = -1;
d2Lvec(uncens) = base(uncens) .* (k+1) + 1;
d2Lvec(cens)   = base(cens)   .* k;
d2Laa = c./alpha.^2 .* sum(freq .* d2Lvec);

% c-c
base(safe) = xac(safe).*lxa(safe).^2./xac1_2(safe);
base(overflow) = 0;
d2Lvec = base .* k;
d2Lvec(uncens) = d2Lvec(uncens)  + base(uncens);
d2Lcc = -n_uncens./c.^2 - sum(freq .* d2Lvec);

% k-k
d2Lkk = -n_uncens./k.^2;

% Off-diagonal elements
base(safe) = xac(safe) .* (c.*lxa(safe) + xac1(safe)) ./ xac1_2(safe); 
base(overflow) = 1; 
d2Lvec(uncens) = (k + 1) .* base(uncens) - 1;
d2Lvec(cens)   =  k      .* base(cens);
d2Lac = sum(freq./alpha .* d2Lvec);

d2Lvec(safe) = xac(safe)./xac1(safe);
d2Lvec(overflow) = 1;
d2Lak =  (c./alpha) .* sum(freq.*d2Lvec);

d2Lvec = d2Lvec.*lxa;
d2Lck = -sum(freq.*d2Lvec);

% Form the empirical Fisher Information Matrix
Ihat = -[ d2Laa d2Lac d2Lak;
          d2Lac d2Lcc d2Lck;
          d2Lak d2Lck d2Lkk; ];

% Straight from mlecov
% Make sure the Hessian is pos def, refuse to compute the cov matrix if not.
[R,p] = chol(Ihat);
if p > 0
    warning(message('stats:addburr:NonPosDefHessian'));
    acov = NaN(length(phat));
else
    % The asymptotic cov matrix approximation is the negative inverse of the
    % Hessian.
    Rinv = inv(R);
    acov = Rinv*Rinv';
end
end

function [nL, dnL] = burr_nloglf_obj(parms, x, cens, freq)
%BURR_NLOGLF Objective function for Burr maximum likelihood.
%            with closed-form k derived from alpha, c

alpha = parms(1);
c = parms(2);

% K has a closed form
lkh = logkhat(x, alpha, c, cens, freq);

% Calculate the log likelihood using all 3 parameters.
% This approach uses unique simplifications that only apply to khat, and 
% therefore cannot be equated with the original burrlike
xa = x ./ alpha;
lxa = log(xa);

% Careful Logarithm: log(1 + x) = x for x < eps
l1_xac = log1p(xa .^ c);

% This has a tendency to overflow, but for large x,
% 1 + x $\approx$ x, and so log1p(xa^c) $\approx$ c*log(xa)
overflow = isinf(l1_xac);
l1_xac(overflow) = c.*log(xa(overflow));

% Normalize the censored flags, just in case they weren't
% logicals to begin with
uncens = ~cens;

% Frequently calculated sum
n_uncens = sum(freq(uncens));

% Uncensored data log-likelihoods
% Preallocate to prevent reallocation in case last elements are censored
L = zeros(size(x));
L(uncens) = (c - 1)*lxa(uncens) - l1_xac(uncens);

% Note that this is the only place k-hat factors in, as, unlike
% in the generic formula, the value of k-hat cancels with k*sum(l1_xac)
% to give n_uncens

% Add in precalcualted constants, sum and sign flip for the final value
nL = n_uncens * (log(alpha/c) - lkh + 1) - sum(freq .* L);

% Give the gradient in the only two variables that matter (alpha and c)
if nargout > 1
    kh = exp(lkh);
    % Normalize the censored flags, just in case they weren't
    % logicals to begin with
    uncens = ~cens;

    xa = x ./ alpha;
    xac = xa .^ c;

    lxa = log(xa);
    % Frequently calculated sum
    n_uncens = sum(freq(uncens));

    % Precalculate commonly used quotient
    xaq = 1 + 1./xac;

    % Preallocate to prevent reallocation in case last elements are censored
    dLda = zeros(size(x));
    dLdc = zeros(size(x));

    % Partial wrt alpha
    dLda(uncens) = ((kh + 1)./xaq(uncens) - 1);
    dLda(cens)   = kh./(xaq(cens));
    dLda = c/alpha * sum(freq .* dLda);

    % Partial wrt c
    dLdc(uncens) = lxa(uncens).*(1 - (kh + 1)./xaq(uncens));
    dLdc(cens) = -kh*lxa(cens)./xaq(cens);
    dLdc = n_uncens/c + sum(freq .* dLdc);

    % Partial wrt k: not calculated, but left here for reference
    % dLdk = n_uncens/k - sum(freq .* l1_xac);

    % Not fully vectorized to allow for random shapes in
    % dLda and dLdc (i.e., sums are taken for individual
    % components, not along some dimension at the very end
    dnL = -[dLda dLdc];
end
end

function lkh = logkhat(x, alpha, c, cens, freq)
% KHAT Closed form of the MLE estimate of the log of k (khat), as a
%      function of x, alpha, and c

xa = x ./ alpha;

% Careful Logarithm: log(1 + x) = x for x < eps
l1_xac = log1p(xa.^c);

% This has a tendency to overflow, but for large x,
% 1 + x $\approx$ x, and so log1p(xa^c) $\approx$ c*log(xa)
overflow = isinf(l1_xac);
l1_xac(overflow) = c.*log(xa(overflow));

k_denom = sum(freq .* l1_xac);

% In the case where the denominator is tiny, instead of exploding
% we can write an approximate form
if k_denom < eps
    % Since every term in the denominator is >=0, if the denominator is
    % vanishingly small, then it must mean that each term in it is also
    % vanishingly small.

    % This may be because of alpha, in which case, we have
    % a very clear approximation. But before we do that, we
    % have to make sure that using the unscaled data will not explode on us, and,
    % if it does, contain it cleverly (this could happen, for example, when
    % the values in x are large, c, is large, and a is even larger.

    % The trick we use to prevent explosion is just the saddlepoint
    % approximation
    lsxc = log(sum(freq.*(x.^c))); 
    if isinf(lsxc)
        [xmax, xmaxi] = max(x(:));
        lsxc = c*log(freq(xmaxi).*xmax);
    end

    lkh = log(sum(freq(~cens))) + c *log(alpha) - lsxc;
else
    lkh = log(sum(freq(~cens))) - log(k_denom);
end
end
