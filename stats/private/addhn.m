function s = addhn(s)
%ADDhn Add the Half-normal distribution.

%   Copyright 2015-2018 The MathWorks, Inc.

j = length(s) + 1;
s(j).name = getString(message('stats:dfittool:NameHalfNormal')); 
s(j).code = 'half normal';
s(j).pnames = {'mu' 'sigma'};
s(j).pdescription = {'location' 'scale'};
s(j).prequired = [true false];
s(j).fitfunc = @hnfit;
s(j).likefunc = @hnlike;
s(j).cdffunc = @hncdf;
s(j).pdffunc = @hnpdf;
s(j).invfunc = @hninv;
s(j).statfunc = @hnstat;
s(j).randfunc = @hnrnd;
s(j).checkparam = @(p) p(2)>0;
s(j).cifunc = @stathnci;
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = false;
s(j).paramvec = true;
s(j).support = [-Inf Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = true;
s(j).uselogpp = false;
s(j).optimopts = false;
s(j).supportfunc = [];


% ==== Half normal distribution functions ====
% these distribution functions do not yet handle arrays of parameters

function [varargout] = hncdf(x,varargin)
%HNCDF Half-normal cumulative distribution function (cdf).

if nargin<1
   error(message('stats:addhn:TooFewInputsX'));
end
if nargin>1 && strcmpi(varargin{end},'upper')
    %Compute upper tail
    uflag=true;
    varargin(end)= [];
elseif nargin>1 && ischar(varargin{end})&& ~strcmpi(varargin{end},'upper')
    error(message('stats:cdf:UpperTailProblem'));
else
    uflag=false;
end
[varargout{1:max(1,nargout)}] = localhncdf(uflag,x,varargin{:});

function p = localhncdf(uflag,x,mu,sigma)
if nargin < 3
    mu = 0;
end
if nargin < 4
    sigma = 1;
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    z = (x-mu) ./ sigma;
catch ME
    error(message('stats:addhn:InputSizeMismatch'));
end

% The support is 0 <= (x-mu)/sigma, force zero below that.
z(z<0) = 0;

if uflag == true
    z = -z;
    p = 1+erf(z./sqrt(2));
else
    p = erf(z./sqrt(2));
end


function y = hnpdf(x,mu,sigma)
%HNPDF Half-normal probability density function (pdf).

if nargin<1
    error(message('stats:addhn:TooFewInputsX'));
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

try
    z = (x-mu) ./ sigma;
catch ME
    error(message('stats:addhn:InputSizeMismatch'));
end

y = sqrt(2/pi)./sigma.*exp(-0.5 * z.^2);

% The support is 0 <= (x-mu)/sigma, force zero below that.
y(z<0) = 0;


function x = hninv(p,mu,sigma)
%HNINV Inverse of the half-normal cumulative distribution function (cdf).

if nargin<1
    error(message('stats:addhn:TooFewInputsP'));
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end

% Return NaN for out of range parameters or probabilities.
sigma(sigma <= 0) = NaN;
p(p < 0 | 1 < p) = NaN;

x0 = erfinv(p);
try
    x = sqrt(2)*sigma.*x0 + mu;
catch
    error(message('stats:addhn:InputSizeMismatch'));
end


function r = hnrnd(mu,sigma,varargin)
%HNRND Random arrays from the half-normal distribution.

if nargin < 2
    error(message('stats:addhn:TooFewInputs'));
end

[err, sizeOut] = internal.stats.statsizechk(2, mu, sigma, varargin{:});
if err > 0
    error(message('stats:addhn:InputSizeMismatch'));
end

% Return NaN for elements corresponding to illegal parameter values.
sigma(sigma < 0) = NaN;

ty = internal.stats.dominantType(mu, sigma);
r = abs(randn(sizeOut,'like',ty)) .* sigma + mu;


function [m,v] = hnstat(mu,sigma)
%HNSTAT Mean and variance for the half-normal distribution.

if nargin < 2
    error(message('stats:addhn:TooFewInputs'));
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;
s2 = sigma .^ 2;

try
    m = mu + sqrt(2/pi)*sigma;
    v = (1-2/pi)*s2;
catch
    error(message('stats:addhn:InputSizeMismatch'));
end

% ==== Half normal fitting functions ====

function [phat,pci] = hnfit(x,mu,alpha)
%HNFIT Parameter estimates and confidence intervals for half-normal value data.
%   HNFIT does not estimate MU, and it must be assumed known.

% Illegal data return an error.
if any(x<mu)
    error(message('stats:addhn:BadDataLocation'));
end

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
end

n = size(x,1);

x = x-mu;
sigmahat = sqrt(sum(x.*x)./n);

if nargout < 2
    phat = [mu sigmahat];
else
    phat = [mu sigmahat];
    pci = stathnci(phat,alpha,x);
end


function parmci=stathnci(parmhat,alpha,x)
%STATSHNCI Half-normal parameter confidence interval.

% Number of observations
n = size(x,1);

mu = parmhat(1);
sigma = parmhat(2);

% Confidence interval
if n > 0
    chi2crit = chi2inv([alpha/2 1-alpha/2],n);
    parmci = [sigma*sqrt(n./chi2crit(2)); ...
            sigma*sqrt(n./chi2crit(1))];
    parmci = [[mu;mu] parmci];
else
    parmci = NaN(2,numel(parmhat),'like',x);
end


function [nlogL,acov] = hnlike(params,data)
%HNLIKE Negative log-likelihood for the half-normal distribution.

if nargin < 2
    error(message('stats:addhn:TooFewInputs'));
elseif numel(data) > length(data)
    error(message('stats:addhn:VectorRequired'));
end

mu = params(1);
sigma = params(2);

% Return NaN for out of range parameter or data.
sigma(sigma <= 0) = NaN;
data(data<mu) = NaN;
z = (data-mu) ./ sigma;

% Sum up the individual log-likelihood terms, and return the negative
% log-likelihood.
logL = -.5.*z.*z - log(sqrt(pi./2).*sigma);
nlogL = -sum(logL);

if nargout == 2
    nH = -sum(1 - 3.*z.*z);
    avar =  (sigma.^2) ./ nH;
    acov = [[0 0]; [0 avar]];
end