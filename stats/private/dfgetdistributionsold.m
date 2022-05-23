function [s,emsg] = dfgetdistributionsold(distname)
%DFGETDISTRIBUTIONSOLD Get old-style structure defining the distributions supported by dfittool


%   Copyright 2012-2015 The MathWorks, Inc.

% If a struct was passed in return
emsg = '';
if nargin>0 && isstruct(distname)
   return
end

% Get built-in distributions
s = getbuiltins;

% Get user-written distributions
try
    s = dfgetuserdists(s);
catch ME
    if nargout>=2
        emsg = ME.message;
        causes = ME.cause;
        for j=1:length(causes)
            emsg = sprintf('%s\n%s',emsg,causes{j}.message);
        end
    else
        rethrow(ME)
    end
end

% Sort by name
lowernames = lower(char(s.name));
[~, ind] = sortrows(lowernames);
s = s(ind);

if nargin>0 && ~isempty(distname)
   % Return only the distribution(s) requested, not all of them
   allnames = {s.code};
   distnum = strncmpi(distname, allnames, length(distname));
   s = s(distnum);
end


% ------------------------------------
function s = getbuiltins
%GETBUILTINS Get distributions functions provided off the shelf

ndists = 11;     % to be updated if distributions added or removed
s(ndists).name = '';

% Exponential distribution
j = 1;
s(j).name = getString(message('stats:dfittool:NameExponential'));      % distribution name
s(j).code = 'exponential';      % distribution code name
s(j).pnames = {'mu'};           % parameter names
s(j).pdescription = {'scale'};  % parameter descriptions
s(j).prequired = false;         % is a value required for this parameter?
s(j).fitfunc = @expfit;         % fitting function
s(j).likefunc = @explike;       % likelihood (and covariance) function
s(j).cdffunc = @expcdf;         % cdf function
s(j).pdffunc = @exppdf;         % pdf function
s(j).invfunc = @expinv;         % inverse cdf function
s(j).statfunc = @expstat;       % function to compute mean and var
s(j).randfunc = @exprnd;        % function to generate random numbers
s(j).checkparam = @(p) p>0;     % function to check for valid parameter
s(j).cifunc = @statexpci;       % function to compute conf intervals
s(j).loginvfunc = [];           % inverse cdf function on log scale, if any
s(j).logcdffunc = [];           % cdf function on log scale, if any
s(j).hasconfbounds = true;      % supports conf bnds for cdf and inverse
s(j).censoring = true;          % supports censoring
s(j).paramvec = true;           % returns fitted parameters as a vector
s(j).support = [0 Inf];         % range of x with positive density
s(j).closedbound = [true false];% is x at this boundary point acceptable
s(j).iscontinuous = true;       % is continuous, not discrete
s(j).islocscale = true;         % is location/scale family, no shape param
s(j).uselogpp = false;          % use log scale for probability plot
s(j).optimopts = false;         % fitting routine accepts an optim structure
s(j).supportfunc = [];          % compute support given specific parameters

% Extreme value
j = j+1;
s(j).name = getString(message('stats:dfittool:NameExtremeValue'));
s(j).code = 'extreme value';
s(j).pnames = {'mu' 'sigma'};
s(j).pdescription = {'location' 'scale'};
s(j).prequired = [false false];
s(j).fitfunc = @evfit;
s(j).likefunc = @evlike;
s(j).cdffunc = @evcdf;
s(j).pdffunc = @evpdf;
s(j).invfunc = @evinv;
s(j).statfunc = @evstat;
s(j).randfunc = @evrnd;
s(j).checkparam = @(p) p(2)>0;
s(j).cifunc = @(p,cv,a,x,c,f) statparamci(p,cv,a,[false,true]);
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = true;
s(j).censoring = true;
s(j).paramvec = true;
s(j).support = [-Inf Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = true;
s(j).uselogpp = false;
s(j).optimopts = true;
s(j).supportfunc = [];

% Gamma
j = j+1;
s(j).name = getString(message('stats:dfittool:NameGamma'));
s(j).code = 'gamma';
s(j).pnames = {'a' 'b'};
s(j).pdescription = {'shape' 'scale'};
s(j).prequired = [false false];
s(j).fitfunc = @gamfit;
s(j).likefunc = @gamlike;
s(j).cdffunc = @gamcdf;
s(j).pdffunc = @gampdf;
s(j).invfunc = @gaminv;
s(j).statfunc = @gamstat;
s(j).randfunc = @gamrnd;
s(j).checkparam = @(p) all(p>0);
s(j).cifunc = @(p,cv,a,x,c,f) statparamci(p,cv,a,[true,true]);
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = true;
s(j).censoring = true;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = false;
s(j).uselogpp = false;
s(j).optimopts = true;
s(j).supportfunc = [];

% Lognormal
j = j+1;
s(j).name = getString(message('stats:dfittool:NameLognormal'));
s(j).code = 'lognormal';
s(j).pnames = {'mu' 'sigma'};
s(j).pdescription = {'log location' 'log scale'};
s(j).prequired = [false false];
s(j).fitfunc = @lognfit;
s(j).likefunc = @lognlike;
s(j).cdffunc = @logncdf;
s(j).pdffunc = @lognpdf;
s(j).invfunc = @logninv;
s(j).statfunc = @lognstat;
s(j).randfunc = @lognrnd;
s(j).checkparam = @(p) p(2)>0;
s(j).cifunc = @(p,cv,a,x,c,f) reshape(statnormci(p,cv,a,x,c,f),[2 2]);
s(j).loginvfunc = @norminv;
s(j).logcdffunc = @normcdf;
s(j).hasconfbounds = true;
s(j).censoring = true;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = true;
s(j).uselogpp = true;
s(j).optimopts = true;
s(j).supportfunc = [];

% Normal
j = j+1;
s(j).name = getString(message('stats:dfittool:NameNormal'));
s(j).code = 'normal';
s(j).pnames = {'mu' 'sigma'};
s(j).pdescription = {'location' 'scale'};
s(j).prequired = [false false];
s(j).fitfunc = @normfit;
s(j).likefunc = @normlike;
s(j).cdffunc = @normcdf;
s(j).pdffunc = @normpdf;
s(j).invfunc = @norminv;
s(j).statfunc = @normstat;
s(j).randfunc = @normrnd;
s(j).checkparam = @(p) p(2)>0;
s(j).cifunc = @(p,cv,a,x,c,f) reshape(statnormci(p,cv,a,x,c,f),[2 2]);
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = true;
s(j).censoring = true;
s(j).paramvec = false;
s(j).support = [-Inf Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = true;
s(j).uselogpp = false;
s(j).optimopts = true;
s(j).supportfunc = [];

% Weibull
j = j+1;
s(j).name = getString(message('stats:dfittool:NameWeibull'));
s(j).code = 'weibull';
s(j).pnames = {'a' 'b'};
s(j).pdescription = {'scale' 'shape'};
s(j).prequired = [false false];
s(j).fitfunc = @wblfit;
s(j).likefunc = @wbllike;
s(j).cdffunc = @wblcdf;
s(j).pdffunc = @wblpdf;
s(j).invfunc = @wblinv;
s(j).statfunc = @wblstat;
s(j).randfunc = @wblrnd;
s(j).checkparam = @(p) all(p>0);
s(j).cifunc = @(p,cv,a,x,c,f) statparamci(p,cv,a,[true,true]);
s(j).loginvfunc = @evinv;
s(j).logcdffunc = @evcdf;
s(j).hasconfbounds = true;
s(j).censoring = true;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = true;
s(j).uselogpp = true;
s(j).optimopts = true;
s(j).supportfunc = [];

% Rayleigh
j = j+1;
s(j).name = getString(message('stats:dfittool:NameRayleigh'));
s(j).code = 'rayleigh';
s(j).pnames = {'b'};
s(j).pdescription = {'scale'};
s(j).prequired = false;
s(j).fitfunc = @raylfit;
s(j).likefunc = [];
s(j).cdffunc = @raylcdf;
s(j).pdffunc = @raylpdf;
s(j).invfunc = @raylinv;
s(j).statfunc = @raylstat;
s(j).randfunc = @raylrnd;
s(j).checkparam = @(p) p>0;
s(j).cifunc = @(p,cv,a,x,c,f) sqrt(2*numel(x)*p.^2 ./ chi2inv([1-a/2; a/2], 2*numel(x)));
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = false;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = true;
s(j).uselogpp = false;
s(j).optimopts = false;
s(j).supportfunc = [];

% Poisson
j = j+1;
s(j).name = getString(message('stats:dfittool:NamePoisson'));
s(j).code = 'poisson';
s(j).pnames = {'lambda'};
s(j).pdescription = {'mean'};
s(j).prequired = false;
s(j).fitfunc = @poissfit;
s(j).likefunc = [];
s(j).cdffunc = @poisscdf;
s(j).pdffunc = @poisspdf;
s(j).invfunc = @poissinv;
s(j).statfunc = @poisstat;
s(j).randfunc = @poissrnd;
s(j).checkparam = @(p) p>=0;
s(j).cifunc = @(p,cv,a,x,c,f) statpoisci(numel(x),p,a);
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = false;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [true false];
s(j).iscontinuous = false;
s(j).islocscale = false;
s(j).uselogpp = false;
s(j).optimopts = false;
s(j).supportfunc = [];

% Negative binomial
j = j+1;
s(j).name = getString(message('stats:dfittool:NameNegativeBinomial'));
s(j).code = 'negative binomial';
s(j).pnames = {'r' 'p'};
s(j).pdescription = {'successes' 'probability'};
s(j).prequired = [false false];
s(j).fitfunc = @nbinfit;
s(j).likefunc = @nbinlike;
s(j).cdffunc = @nbincdf;
s(j).pdffunc = @nbinpdf;
s(j).invfunc = @nbininv;
s(j).statfunc = @nbinstat;
s(j).randfunc = @nbinrnd;
s(j).checkparam = @(p) (p(1)>0 & 0<p(2) & p(2)<=1);
s(j).cifunc = @(p,cv,a,x,c,f) statparamci(p,cv,a);
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = false;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [true false];
s(j).iscontinuous = false;
s(j).islocscale = false;
s(j).uselogpp = false;
s(j).optimopts = true;
s(j).supportfunc = [];

% Beta
j = j+1;
s(j).name = getString(message('stats:dfittool:NameBeta'));
s(j).code = 'beta';
s(j).pnames = {'a' 'b'};
s(j).pdescription = {'a' 'b'};
s(j).prequired = [false false];
s(j).fitfunc = @betafit;
s(j).likefunc = @betalike;
s(j).cdffunc = @betacdf;
s(j).pdffunc = @betapdf;
s(j).invfunc = @betainv;
s(j).statfunc = @betastat;
s(j).randfunc = @betarnd;
s(j).cifunc = @(p,cv,a,x,c,f) statparamci(p,cv,a,[true,true]);
s(j).checkparam = @(p) all(p>0);
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = false;
s(j).paramvec = true;
s(j).support = [0 1];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = false;
s(j).uselogpp = false;
s(j).optimopts = false;
s(j).supportfunc = [];

% Binomial
j = j+1;
s(j).name = getString(message('stats:dfittool:NameBinomial'));
s(j).code = 'binomial';
s(j).pnames = {'N' 'p'};
s(j).pdescription = {'trials' 'probability'};
s(j).prequired = [true false];
s(j).fitfunc = @localbinofit;
s(j).likefunc = @localbinolike;
s(j).cdffunc = @binocdf;
s(j).pdffunc = @binopdf;
s(j).invfunc = @binoinv;
s(j).statfunc = @binostat;
s(j).randfunc = @binornd;
s(j).checkparam = @(p) (p(1)>=0 & p(1)==round(p(1)) & p(2)>=0 & p(2)<=1);
s(j).cifunc = @(p,cv,a,x,c,f) [[p(1);p(1)],statbinoci(p(1)*p(2)*numel(x),p(1)*numel(x),a)'];
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = false;
s(j).paramvec = true;
s(j).support = [0 Inf];
s(j).closedbound = [true false];
s(j).iscontinuous = false;
s(j).islocscale = false;
s(j).uselogpp = false;
s(j).optimopts = false;
s(j).supportfunc = @localbinosupport;

% Generalized extreme value
j = j+1;
s(j).name = getString(message('stats:dfittool:NameGeneralizedExtremeValue'));
s(j).code = 'generalized extreme value';
s(j).pnames = {'k' 'sigma' 'mu'};
s(j).pdescription = {'shape' 'scale' 'location'};
s(j).prequired = [false false false];
s(j).fitfunc = @gevfit;
s(j).likefunc = @gevlike;
s(j).cdffunc = @gevcdf;
s(j).pdffunc = @gevpdf;
s(j).invfunc = @gevinv;
s(j).statfunc = @gevstat;
s(j).randfunc = @gevrnd;
s(j).checkparam = @(p) p(2)>0;
s(j).cifunc = @(p,cv,a,x,c,f) statparamci(p,cv,a,[false,true,false]);
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = false;
s(j).paramvec = true;
s(j).support = [-Inf Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = false;
s(j).uselogpp = false;
s(j).optimopts = true;
s(j).supportfunc = @localgevsupport;

% Generalized Pareto
j = j+1;
s(j).name = getString(message('stats:dfittool:NameGeneralizedPareto'));
s(j).code = 'generalized pareto';
s(j).pnames = {'k' 'sigma' 'theta'};
s(j).pdescription = {'shape' 'scale' 'threshold'};
s(j).prequired = [false false true];
s(j).fitfunc = @localgpfit;
s(j).likefunc = @localgplike;
s(j).cdffunc = @gpcdf;
s(j).pdffunc = @gppdf;
s(j).invfunc = @gpinv;
s(j).statfunc = @gpstat;
s(j).randfunc = @gprnd;
s(j).checkparam = @(p) p(2)>0;
s(j).cifunc = @(p,cv,a,x,c,f) statparamci(p,cv,a,[false,true]);
s(j).loginvfunc = [];
s(j).logcdffunc = [];
s(j).hasconfbounds = false;
s(j).censoring = false;
s(j).paramvec = true;
s(j).support = [-Inf Inf];
s(j).closedbound = [false false];
s(j).iscontinuous = true;
s(j).islocscale = false;
s(j).uselogpp = false;
s(j).optimopts = true;
s(j).supportfunc = @localgpsupport;

s = addbisa(s);
s = addinvg(s);
s = addlogi(s);
s = addnaka(s);
s = addtls(s);
s = addrice(s);
s = addburr(s);
s = addhn(s);

% ------------ binomial function is a special case
function [phat,pci] = localbinofit(x,N,alpha)
%LOCALBINOFIT Version of binofit that operates on vectors

if nargin<3
    alpha = 0.05;
end
nx = length(x);
sumx = sum(x);
sumN = nx * N;
if nargout==2
   [phat,pci] = binofit(sumx,sumN,alpha);
else
   phat = binofit(sumx,sumN,alpha);
end

phat = [N phat];
if nargout==2
   pci = [NaN NaN; pci]';
end

function [nlogl,pcov] = localbinolike(params,x)
%LOCALBINOLIKE Binomial likelihood

N = params(1);
p = params(2);

nlogl = -sum(log(binopdf(x,N,p)));

pcov = [0 0;0 p*(1-p)/(N*length(x))];

function [range,closed] = localbinosupport(params)
range = [0, params(1)];
closed = [true true];

% ------------ generalized Pareto functions are a special case
function [phat,pci] = localgpfit(x,theta,alpha,varargin)
%LOCALGPFIT Version of gpfit that handles a fixed threshold param

if any(x<=theta)
    error(message('stats:gpfit:BadDataThreshold'));
end

if nargout < 2
    phat = [gpfit(x-theta,alpha,varargin{:}) theta];
else
    [phat,pci] = gpfit(x-theta,alpha);
    phat = [phat theta];
    pci = [pci [theta; theta]];
end

function [nlogL,acov] = localgplike(params,data)
%LOCALGPLIKE Version of gplike that handles a fixed threshold param

theta = params(3);
params = params(1:2);
if nargout < 2
    nlogL = gplike(params,data-theta);
else
    [nlogL,acov] = gplike(params,data-theta);
    acov = [acov [0; 0]; [0 0 0]];
end

function [range,closed] = localgpsupport(params)
k = params(1);
theta = params(3);
if k<0
    sigma = params(2);
    range = sort([theta, theta-sigma/k]);
    closed = [true true];
else
    range = [theta Inf];
    closed = [false false];
end

% ------------ generalized extreme value functions are a special case
function [range,closed] = localgevsupport(params)
k = params(1);
sigma = params(2);
mu = params(3);

if k==0
    range = [-Inf Inf];
    closed = [false false];
elseif k>0
    range = [mu-sigma/k, Inf];
    closed = [true false];
else
    range = [-Inf, mu-sigma/k];
    closed = [false true];
end
