function s = addrice(s)
%ADDRICE Add the Rician distribution.

%   Copyright 1993-2018 The MathWorks, Inc.


j = length(s) + 1;
s(j).name = getString(message('stats:dfittool:NameRician'));
s(j).code = 'rician';
s(j).pnames = {'s' 'sigma'};
s(j).pdescription = {'noncentrality' 'scale'};
s(j).prequired = [false false];
s(j).fitfunc = @ricefit;
s(j).likefunc = @ricelike;
s(j).cdffunc = @ricecdf;
s(j).pdffunc = @ricepdf;
s(j).invfunc = @riceinv;
s(j).statfunc = @ricestat;
s(j).randfunc = @ricernd;
s(j).checkparam = @(p) p(1)>=0 & p(2)>0;
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

% ==== Rician distribution functions ====

% these distribution functions do not yet handle arrays of parameters

function y = ricepdf(x,s,sigma)
%RICEPDF Rician probability density function (pdf).
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

x(x<0) = 0;
sigsq = sigma.^2;
rsq = (x.^2 + s.^2)./(2.*sigsq);
z = x./sigsq;
expon = rsq - z.*s;
y = z .* exp(-expon) .* besseli(0,z.*s,1);
y(expon > (log(realmax(class(x)))-1)) = 0; % fix up 0*Inf


function p = ricecdf(x,s,sigma)
%RICECDF Rician cumulative distribution function (cdf).
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

x(x<0) = 0;
p = ncx2cdf((x./sigma).^2, 2, (s./sigma).^2);


function x = riceinv(p,s,sigma)
%RICEINV Inverse of the Rician cumulative distribution function (cdf).
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

x = sigma .* sqrt(ncx2inv(p, 2, (s./sigma).^2));


function r = ricernd(s,sigma,varargin)
%RICERND Random arrays from the Rician distribution.
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

[err, sizeOut] = internal.stats.statsizechk(2,s,sigma,varargin{:});
if err > 0
    error(message('stats:ricernd:InconsistentSizes'));
end

r = sigma .* sqrt(ncx2rnd(2, (s./sigma).^2, sizeOut));


function [m,v] = ricestat(s,sigma)
%RICESTAT Mean and variance for the Rician distribution.
s(s < 0) = NaN;
sigma(sigma <= 0) = NaN;

t = .5 .* (s./sigma).^2;
m = sigma.*sqrt(.5.*pi) .* ((1+t).*besseli(0,.5.*t,1) + t.*besseli(1,.5.*t,1));
v = 2.*sigma.^2 + s.^2 - m.^2;


function [nlogL,acov] = ricelike(params,data,cens,freq)
%RICELIKE Negative log-likelihood for the Rician distribution.
if nargin < 4 || isempty(freq), freq = ones(size(data)); end
if nargin < 3 || isempty(cens), cens = zeros(size(data)); end

nlogL = rice_nloglf(params, data, cens, freq);
if nargout > 1
    acov = mlecov(params, data, 'nloglf',@rice_nloglf, 'cens',cens, 'freq',freq);
end


% ==== Rician fitting functions ====

function [phat,pci] = ricefit(x,alpha,cens,freq,opts)
%RICEFIT Parameter estimates and confidence intervals for Rician data.

if nargin < 2 || isempty(alpha), alpha = .05; end
if nargin < 3 || isempty(cens), cens = zeros(size(x)); end
if nargin < 4 || isempty(freq), freq = ones(size(x)); end
if nargin < 5, opts = []; end

if any(x <= 0)
    error(message('stats:ricefit:BadData'));
end

% Moment estimators of the uncensored data as starting point
% E[x.^2] = s.^2 + 2.*sigma.^2
% E[x.^4] = s.^4 + 8.*s.^2.*sigma.^2 + 8.*sigma.^4
xsqunc = x(cens == 0).^2;
meanxsq = mean(xsqunc); meanx4th = mean(xsqunc.^2);
if meanxsq^2 < meanx4th && meanx4th < 2*meanxsq^2
    s4th = 2*meanxsq^2 - meanx4th;
    ssq = sqrt(s4th);
    sigsq = .5*(meanxsq - ssq);
    start = [sqrt(ssq) sqrt(sigsq)];
else
    start = cast([1 1],class(x));
end

% The default options include turning fminsearch's display off.  This
% function gives its own warning/error messages, and the caller can turn
% display on to get the text output from fminsearch if desired.
options = statset(statset('ricefit'), opts);
tolBnd = options.TolBnd;
options = optimset(options);

% Maximize the log-likelihood with respect to mu and sigma.
dfltOptions = struct('DerivativeCheck','off', 'HessMult',[], ...
    'HessPattern',ones(2,2), 'PrecondBandWidth',Inf, ...
    'TypicalX',ones(2,1,class(x)), 'MaxPCGIter',1, 'TolPCG',0.1);
funfcn = {'fungrad' 'ricefit' @rice_nloglf [] []};
[phat, ~, ~, err, output] = ...
         statsfminbx(funfcn, start, [-Inf; tolBnd], [Inf; Inf], ...
                     options, dfltOptions, 1, x, cens, freq);
if (err == 0)
    % fminsearch may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if output.funcCount >= options.MaxFunEvals
        warning(message('stats:ricefit:EvalLimit'));
    else
        warning(message('stats:ricefit:IterLimit'));
    end
elseif (err < 0)
    error(message('stats:ricefit:NoSolution'));
end

% Compute CIs using a normal approximation for phat.
if nargout > 1
    acov = mlecov(phat, x, 'nloglf',@rice_nloglf, 'cens',cens, 'freq',freq);
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';
    pci = norminv([probs probs], [phat; phat], [se; se]);
end


function [nll,ngrad] = rice_nloglf(parms, x, cens, freq)
%RICE_NLOGLF Objective function for Rician maximum likelihood.
s = parms(1);
sigma = parms(2);

theta = s/sigma;
z = x./sigma;
ztheta = z.*theta;
bess0 = besseli(0, ztheta, 1); % I0(z.*theta)*exp(-z.*theta)
rsq = (z.^2 + theta.^2)./2;
L = -rsq + log(bess0) + log(z./sigma) + ztheta;
ncen = sum(freq.*cens);
if ncen > 0
    cen = (cens == 1);
    zcen = z(cen);
    Q = marcumq(theta,zcen);
    L(cen) = log(Q);
end
nll = -sum(freq .* L);

if nargout > 1
    bess1 = besseli(1, ztheta, 1);
    dlogbess0 = bess1 ./ bess0;
    dL1 = (-theta + dlogbess0.*z) ./ sigma;
    dL2 = -2 * (1 - rsq + dlogbess0.*ztheta) ./ sigma;
    if ncen > 0
        t = exp(-rsq(cen) + ztheta(cen));
        dQdtheta = zcen.*bess1(cen).*t;
        dQdz = -zcen.*bess0(cen).*t;
        dtheta1 = 1./sigma;
        dtheta2 = -theta./sigma;
        % dz1 = 0;
        dz2 = -zcen./sigma;
        dL1(cen) = dQdtheta.*dtheta1 ./ Q;
        dL2(cen) = (dQdtheta.*dtheta2 + dQdz.*dz2) ./ Q;
    end
    ngrad = -[sum(freq .* dL1) sum(freq .* dL2)];
end


function Q = marcumq(a,b)
% Q = MARCUMQ(A,B) returns Marcum's "Q" function.

if isa(a,'single') || isa(b,'single')
   Q = NaN(size(b),'single');
else
   Q = NaN(size(b));
end
Q(a~=Inf & b==0) = 1;
Q(a~=Inf & b==Inf) = 0;
Q(a==Inf & b~=Inf) = 1;
z = (isnan(Q) & a==0 & b~=Inf);
if (any(z))
   Q(z) = exp((-b(z).^2)./2);
end

z = isnan(Q) & ~isnan(a) & ~isnan(b);
if (any(z(:)))
%    aa = (a(z).^2)./2;
   aa = (a.^2)./2;
   bb = (b(z).^2)./2;

   d = exp(-aa);
   h = d;
   f = bb.*exp(-bb);
   k = 1;
   delta = f .* h;
   sum = delta;
   j = (delta > sum.*eps(class(delta)));
   while any(j)
      d = aa.*d./k;
      h = h + d;
      f = bb.*f./(k+1);
      delta = f .* h;
      sum(j) = sum(j) + delta(j);
      j = (delta > sum.*eps(class(delta)));
      k = k + 1;
   end
   Q(z) = 1 - sum;
end
