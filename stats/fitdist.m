function [pd,gn,gl] = fitdist(x,distname,varargin)
%FITDIST Fit probability distribution to data.
%   PD = FITDIST(X,DISTNAME) fits the probability distribution DISTNAME to
%   the data in the column vector X, and returns an object PD representing
%   the fitted distribution.  PD is an object in a class derived from the
%   ProbabilityDistribution class.
%
%   DISTNAME can be 'kernel' to fit a nonparametric kernel-smoothing
%   distribution, or it can be any of the following parametric distribution
%   names:
%
%         'beta'                             Beta
%         'binomial'                         Binomial
%         'birnbaumsaunders'                 Birnbaum-Saunders
%         'burr'                             Burr Type XII
%         'exponential'                      Exponential
%         'extreme value', 'ev'              Extreme value
%         'gamma'                            Gamma
%         'generalized extreme value', 'gev' Generalized extreme value
%         'generalized pareto', 'gp'         Generalized Pareto
%         'half normal', 'hn'                Half Normal
%         'inversegaussian'                  Inverse Gaussian
%         'logistic'                         Logistic
%         'loglogistic'                      Log-logistic
%         'lognormal'                        Lognormal
%         'nakagami'                         Nakagami
%         'negative binomial', 'nbin'        Negative binomial
%         'normal'                           Normal
%         'poisson'                          Poisson
%         'rayleigh'                         Rayleigh
%         'rician'                           Rician
%         'stable'                           Stable
%         'tlocationscale'                   t location-scale
%         'weibull', 'wbl'                   Weibull
%
%   [PDCA,GN,GL] = FITDIST(X,DISTNAME,'BY',G) takes a grouping variable G,
%   fits the specified distribution to the X data from each group, and
%   returns a cell array PDCA of the fitted probability distribution
%   objects.  See "help groupingvariable" for more information.  G can also
%   be a cell array of multiple grouping variables.  GN is a cell array of
%   group labels.  GL is a cell array of grouping variable levels, with one
%   column for each grouping variable.
%
%   PD = FITDIST(..., 'NAME1',VALUE1,'NAME2',VALUE2,...) specifies optional
%   argument name/value pairs chosen from the following list. Argument
%   names are case insensitive and partial matches are allowed.
%
%      Name           Value
%      'censoring'    A boolean vector of the same size as X, containing
%                     ones when the corresponding elements of X are
%                     right-censored observations and zeros when the
%                     corresponding elements are exact observations.
%                     Default is all observations observed exactly.
%                     Censoring is not supported for all distributions.
%      'frequency'    A vector of the same size as X, containing
%                     non-negative integer frequencies for the
%                     corresponding elements in X.  Default is one
%                     observation per element of X.
%      'options'      A structure created by STATSET to specify control
%                     parameters for the iterative fitting algorithm
%
%   For the 'binomial' distribution only:
%      'ntrials'      A positive integer specifying the N parameter (number
%                     of trials).  Not allowed for other distributions.
%
%   For the 'generalized pareto' distribution only:
%      'theta'        The value of the THETA (threshold) parameter for
%                     the generalized Pareto distribution.  Default is 0.
%                     Not allowed for other distributions.
%
%   For the 'half normal' distribution only:
%      'mu'           The value of the MU (location) parameter for the 
%                     half normal distribution.  Default is 0.  Not allowed
%                     for other distributions.
%
%   For the 'kernel' distribution only:
%      'kernel'       The type of kernel smoother to use, chosen from among
%                     'normal' (default), 'box', 'triangle', and
%                     'epanechnikov'.
%      'support'      Either 'unbounded' (default) if the density can
%                     extend over the whole real line, or 'positive' to
%                     restrict it to positive values, or a two-element
%                     vector giving finite lower and upper limits for the
%                     support of the density.
%      'width'        The bandwidth of the kernel smoothing window.  The
%                     default is optimal for estimating normal densities,
%                     but you may want to choose a smaller value to reveal
%                     features such as multiple modes.
%
%   FITDIST treats NaNs as missing values, and removes them.
%
%   Examples:
%        % Fit MPG data using a kernel smooth density estimate
%        load carsmall
%        ksd = fitdist(MPG,'kernel')
%
%        % Fit separate Weibull distributions for each country of origin.
%        % Cell array is empty for countries with insufficient data.
%        wd = fitdist(MPG,'weibull', 'by',Origin)
%
%   See also MAKEDIST, GROUPINGVARIABLE, MLE.

%   Copyright 2008-2018 The MathWorks, Inc.

if nargin > 1
    distname = convertStringsToChars(distname);
end
[varargin{:}] = convertStringsToChars(varargin{:});

if nargin==0
    % Special case, list available distributions
    pd = prob.ProbabilityDistributionRegistry.list('fittable')';
    return
elseif nargin==1 && strcmpi(x,'-reset')
    % Special case, reset the list of distributions
    dfgetset('alldistributions','');       % clear all distributions
    prob.ProbabilityDistributionRegistry.refresh;
    return
end

% Error checking
narginchk(2,Inf);
internal.stats.checkNotTall(upper(mfilename),0,x,distname,varargin{:});

if ~isnumeric(x) || ~isvector(x) || size(x,2)~=1
    error(message('stats:fitdist:BadX'))
elseif any(x==Inf) || any(x==-Inf)
    error(message('stats:fitdist:InfiniteX'))
end
if ~ischar(distname) && ~isstruct(distname)
    error(message('stats:fitdist:BadDist'))
elseif ischar(distname) && strncmpi(distname,'kernel',max(1,length(distname)))
    distname = 'kernel';
end

% Process some args here; others are passed along
pnames = {'by' 'censoring' 'frequency'};
dflts  = {[]   []          []         };
[by,cens,freq,~,args] = internal.stats.parseArgs(pnames,dflts,varargin{:});

n = length(x);
if ~isempty(cens)
    if ~isvector(cens) || ~numel(cens)==n || (~islogical(cens) && ~all(isnan(cens) | cens==0 | cens==1))
        error(message('stats:ProbDistUnivParam:fit:BadCensoring'));
    end
end
if ~isempty(freq)
    if ~isvector(freq) || ~numel(freq)==n || ~all(isnan(freq) | (freq>=0 & freq==round(freq)))
        error(message('stats:ProbDistUnivParam:fit:BadFrequency'));
    end
end

if isstruct(distname)
    distname = distname.code;
end
try
    [dist,~,classname] = dfgetdistributions(distname);
    if isempty(dist)
        error(message('stats:ProbDistUnivParam:checkdistname:UnrecognizedName',distname));
    end
    if isempty(classname)
        % classname is empty if distname is unrecognized. 
        error(message('stats:ProbDistUnivParam:checkdistname:UnrecognizedName',distname));
    else
        fitter = eval(['@' classname '.fit']);
    end
catch ME
    if strcmp(ME.identifier,'stats:internal:getParamVal:BadValueListChoices')
        error(message('stats:ProbDistUnivParam:checkdistname:UnrecognizedName',distname));
    else
        rethrow(ME)
    end
end
if isfield(dist,'fittable')
    if ~dist.fittable
        error(message('stats:fitdist:NotFittable',dist.code))
    end
end
if ~isempty(cens) && any(cens)
    if ~dist.censoring
        error(message('stats:ProbDistUnivParam:fit:CensoringNotAllowed', dist.code));
    elseif (isempty(freq) && all(cens)) || (~isempty(freq) && all(cens(freq>0)))
        error(message('stats:fitdist:AllCensored'));
    end
end

% Fit distribution either singly or by group
if isempty(by)
    if nargout>1
       error(message('stats:fitdist:TooManyOutputs'));
    end
    pd = localfit(dist,fitter,x,cens,freq,args{:});
else
    [gidx,gn,gl] = internal.stats.mgrp2idx(by);
    ngroups = length(gn);
    if isempty(freq)
        freq = ones(size(x));
    end
    if isempty(cens)
        cens = false(size(x));
    end
    
    % Remove NaN and zero-frequency data
    freq(freq==0) = NaN;
    [badin,~,gidx,x,cens,freq] = internal.stats.removenan(gidx,x,cens,freq);
    if badin>0
        error(message('stats:fitdist:InputSizeMismatch'));
    end
   
    pd = cell(1,ngroups);
    for j=1:ngroups
        t = (gidx==j);
        xj = x(t);
        cj = cens(t);
        fj = freq(t);
        try
            pd{j} = localfit(dist,fitter,xj,cj,fj,args{:});
        catch myException
            switch(myException.identifier)
                % Some errors apply across all groups
                case {'stats:ProbDistUnivParam:checkdistname:UnrecognizedName' ...
                      'stats:ProbDistUnivParam:fit:NRequired' ...
                      'stats:ProbDistUnivParam:fit:BadThreshold' ...
                      'stats:internal:parseArgs:BadParamName'}
                    rethrow(myException);
 
                % For other or unanticipated errors, warn and continue
                otherwise
                    warning(message('stats:fitdist:FitError', gn{ j }, myException.message));
            end
        end
    end
end
end

function pd =localfit(spec,fitter,x,c,f,varargin)
prequired = spec.prequired;
nparams = sum(~prequired);
if isempty(f) && numel(x)<nparams
    error(message('stats:ProbDistUnivParam:fit:InsufficientData'));
elseif ~isempty(f) && sum(f>0)<nparams
    error(message('stats:ProbDistUnivParam:fit:InsufficientFreq'));
end

pd = feval(fitter,x,'cens',c,'freq',f,varargin{:});
end
