function s = addburr(s)
%ADDBURR Add the Burr distribution.

%   Copyright 1993-2018 The MathWorks, Inc.


j = length(s) + 1;
s(j).name = getString(message('stats:dfittool:NameBurr')); 
s(j).code = 'burr';
s(j).pnames = {'alpha' 'c' 'k'};
s(j).pdescription = {'scale' 'shape' 'shape'};
s(j).prequired = [false false false];
s(j).fitfunc = @burrfit;
s(j).likefunc = @burrlike;
s(j).cdffunc = @burrcdf;
s(j).pdffunc = @burrpdf;
s(j).invfunc = @burrinv;
s(j).statfunc = @burrstat;
s(j).randfunc = @burrrnd;
s(j).checkparam = @(p) all(p > 0);
s(j).cifunc = @(p,cv,a,x,c,f) statparamci(p,cv,a);
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = true;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = false;
s(j).uselogpp = false;
s(j).optimopts = true;
s(j).supportfunc = [];

% ==== Burr distribution functions ====


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

if nargin > 4 && strcmpi(tail, 'upper')
    % Upper tail calculation
    F = (1 + xac).^-k;
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
        F(normal) = 1 - xac1_nk(normal);
        F_small = zeros(size(F));

        % Crossover range: blend the two formulas during this range
        blend_comparison = xac1_nk; % Use the quantity subtracted from 1 as a guide for when to turn on
                                    % the approximation
        blend_range = [0.999993434530685  0.999982421806361]; % Note: left edge is bigger!
        small       = normal & blend_comparison > blend_range(2); % All small values of F
        if any(small(:))
            blended     = normal & blend_comparison < blend_range(1) & small; % Values over which the blend will happen
            % Linear fade from one to the other
            blending_weights = (blend_comparison(blended) - blend_range(2))./(blend_range(1) - blend_range(2));

            % Calculate all the instances of the lower-tail approximation
            % that will be necessary
            F_small(small) = xac(small) .* k(small) .* (xac(small).*(k(small)-1)./2 + 1) .* xac1_nk(small);

            % Use a weighted blend of the lower tail and the normal formula over the blend range
            F(blended) = blending_weights .* F(blended) + (1 - blending_weights) .* F_small(blended);

            % For anything smaller than the blend range, just use the lower tail approximation
            F(small & ~blended) = F_small(small & ~blended);
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


function [m,v] = burrstat(alpha, c, k)
%BURRSTAT Mean and variance for the Burr distribution.
% Moment r: alpha^r k Beta(r/c + 1, k - r/c) = alpha ^r r/c Beta(r/c, k - r/c)
% if ck > r. Else: Inf. Invalid params: NaN.

[errorcode, alpha, c, k] = distchck(3, alpha, c, k);

if errorcode > 0
    error(message('stats:addburr:InputSizeMismatch'));
end


% Any invalid values: don't bother, leave as NaNs
valid = alpha > 0 & c > 0 & k > 0;

% ck used in determining need to calculate moments
ck = c .* k;

% Preallocate space for mean
m = nan(size(alpha));
hasmean = valid & ck > 1;
% Infinite means
m(valid & ck <= 1) = Inf;

% Precalculate as many quantities as possible
ci = 1./c;
ac = alpha .* ci;
k_ci = k - ci;

% Greatly simplified moments for k == 1
% using the identity that B(x, 1 - x) = Gamma(x)Gamma(1-x) = pi/sin(pi x)
simple = (k == 1);
simplemean = simple & hasmean;
if any(simplemean(:))
    m(simplemean) = ac(simplemean).*pi./sin(pi.*ci(simplemean));
end

normal = ~simple;
normalmean = normal & hasmean;
if any(normalmean(:))
    m(normalmean) = ac(normalmean) .* beta(ci(normalmean), k_ci(normalmean));
end

% Only bother with the variance if it was asked for
if nargout > 1
    v = nan(size(alpha));

    hasvar  = valid & ck > 2;
    % Infinite variances
    v(valid & ck <= 2) = Inf;

    simplevar = simple & hasvar;
    normalvar = normal & hasvar;


    if any(simplevar(:))
        v(simplevar) = 2 * alpha(simplevar) .*ac(simplevar).*pi./sin(2.*pi.*ci(simplevar)) - m(simplevar).^2;
    end
    if any(normalvar(:))
        v(normalvar) = 2 * alpha(normalvar) .* ac(normalvar) .* beta(ci(normalvar) .* 2, k_ci(normalvar) - ci(normalvar)) - m(normalvar).^2;
    end
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
    phat = mlecustom(x_theta, 'nloglf', @burr_nloglf_obj, ...
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
    acov = Rinv*Rinv'; %#ok
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
