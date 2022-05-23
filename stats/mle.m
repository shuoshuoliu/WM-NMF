function [phat, pci] = mle(data,varargin)
%MLE Maximum likelihood estimation.
%   PHAT = MLE(DATA) returns maximum likelihood estimates (MLEs) for the
%   parameters of a normal distribution, computed using the sample data in
%   the vector DATA.
%
%   [PHAT, PCI] = MLE(DATA) returns MLEs and 95% confidence intervals for
%   the parameters.
%
%   [...] = MLE(DATA,'distribution',DIST) computes parameter estimates for
%   the distribution specified by DIST.  DIST is a string scalar or
%   character vector chosen from the following list of distributions:
%
%       'beta'                             Beta
%       'bernoulli'                        Bernoulli
%       'binomial'                         Binomial
%       'birnbaumsaunders'                 Birnbaum-Saunders
%       'burr'                             Burr Type XII
%       'discrete uniform' or 'unid'       Discrete uniform
%       'exponential'                      Exponential
%       'extreme value' or 'ev'            Extreme value
%       'gamma'                            Gamma
%       'generalized extreme value' 'gev'  Generalized extreme value
%       'generalized pareto' or 'gp'       Generalized Pareto
%       'geometric'                        Geometric
%       'half normal' or 'hn'              Half-normal
%       'inversegaussian'                  Inverse Gaussian
%       'logistic'                         Logistic
%       'loglogistic'                      Log-logistic
%       'lognormal'                        Lognormal
%       'nakagami'                         Nakagami
%       'negative binomial' or 'nbin'      Negative binomial
%       'normal'                           Normal
%       'poisson'                          Poisson
%       'rayleigh'                         Rayleigh
%       'rician'                           Rician
%       'stable'                           Stable
%       'tlocationscale'                   t location-scale
%       'uniform'                          Uniform
%       'weibull' or 'wbl'                 Weibull
%
%   [...] = MLE(DATA, ..., 'NAME1',VALUE1,'NAME2',VALUE2,...) specifies
%   optional argument name/value pairs chosen from the following list.
%   Argument names are case insensitive and partial matches are allowed.
%
%        Name           Value
%      'censoring'    A boolean vector of the same size as DATA,
%                     containing ones when the corresponding elements of
%                     DATA are right-censored observations and zeros when
%                     the corresponding elements are exact observations.
%                     Default is all observations observed exactly.
%                     Censoring is not supported for all distributions.
%      'frequency'    A vector of the same size as DATA, containing
%                     non-negative integer frequencies for the corresponding
%                     elements in DATA.  Default is one observation per
%                     element of DATA.
%      'alpha'        A value between 0 and 1 specifying a confidence level
%                     of 100*(1-alpha)% for PCI.  Default is alpha=0.05 for
%                     95% confidence.
%      'ntrials'      A scalar, or a vector of the same size as DATA,
%                     containing the total number of trials for the
%                     corresponding element of DATA.  Applies only to the
%                     binomial distribution.
%      'mu'           A scalar numeric value giving the location parameter
%                     of the half-normal distribution only.
%      'theta'        A scalar numeric value giving the threshold parameter
%                     of the generalized pareto distribuiotn only.
%      'options'      A structure created by a call to STATSET, containing
%                     numerical options for the fitting algorithm.  Not
%                     applicable to all distributions.
%
%   MLE can also fit a custom distribution that you define using
%   distribution functions, in one of three ways:
%
%   [...] = MLE(DATA,'pdf',PDF,'cdf',CDF,'start',START,...) returns MLEs
%   for the parameters of the distribution defined by the probability
%   density and cumulative distribution functions PDF and CDF.  PDF and CDF
%   are function handles created using @.  They accept as inputs a vector
%   of data and one or more individual distribution parameters, and return
%   vectors of probability density values and cumulative probability
%   values, respectively.  If the 'censoring' name/value pair is not
%   present, you may omit the 'cdf' name/value pair.  MLE computes the
%   estimates by numerically maximizing the distribution's log-likelihood,
%   and START is a vector containing initial values for the parameters.
%
%   [...] = MLE(DATA,'logpdf',LOGPDF,'logsf',LOGSF,'start',START,...)
%   returns MLEs for the parameters of the distribution defined by the log
%   probability density and log survival functions LOGPDF and LOGSF. LOGPDF
%   and LOGSF are function handles created using @.  They accept as inputs
%   a vector of data and one or more individual distribution parameters,
%   and return vectors of logged probability density values and logged
%   survival function values, respectively.  This form is sometimes more
%   robust to the choice of starting point than using PDF and CDF
%   functions.  If the 'censoring' name/value pair is not present, you may
%   omit the 'logsf' name/value pair.  START is a vector containing initial
%   values for the distribution's parameters.
%
%   [...] = MLE(DATA,'nloglf',NLOGLF,'start',START,...) returns MLEs for
%   the parameters of the distribution whose negative log-likelihood is
%   given by NLOGLF.  NLOGLF is a function handle specified using @, that
%   accepts the four input arguments
%      PARAMS - a vector of distribution parameter values
%      DATA   - a vector of data
%      CENS   - a boolean vector of censoring values
%      FREQ   - a vector of integer data frequencies
%   NLOGLF must accept all four arguments even if you do not supply the
%   'censoring' or 'frequency' name/value pairs (see above).  However,
%   NLOGLF can safely ignore its CENS and FREQ arguments in that case.
%   NLOGLF returns a scalar negative log-likelihood value and, optionally,
%   a negative log-likelihood gradient vector (see the 'GradObj' STATSET
%   parameter below).  START is a vector containing initial values
%   for the distribution's parameters.
%
%   PDF, CDF, LOGPDF, LOGSF, or NLOGLF can also be cell arrays whose first
%   element is a function handle as defined above, and whose remaining
%   elements are additional arguments to the function.  MLE places these
%   arguments at the end of the argument list in the function call.
%
%   The following optional argument name/value pairs are valid only when
%   'pdf' and 'cdf', 'logpdf' and 'logcdf', or 'nloglf' are given.
%
%      'lowerbound'   A vector the same size as START containing lower bounds
%                     for the distribution parameters.  Default is -Inf.
%      'upperbound'   A vector the same size as START containing upper bounds
%                     for the distribution parameters.  Default is Inf.
%      'optimfun'     A string, either 'fminsearch' or 'fmincon', naming the
%                     optimization function to be used in maximizing the
%                     likelihood.  Default is 'fminsearch'.  You may only
%                     specify 'fmincon' if the Optimization Toolbox is
%                     available.
%
%   When fitting a custom distribution, use the 'options' parameter to
%   control details of the maximum likelihood optimization.  See
%   STATSET('mlecustom') for parameter names and default values.  MLE
%   interprets the following STATSET parameters for custom distribution
%   fitting as follows:
%
%      'GradObj'      'on' or 'off', indicating whether or not FMINCON
%                     can expect the function provided with the 'nloglf'
%                     name/value pair to return the gradient vector of the
%                     negative log-likelihood as a second output.  Default
%                     is 'off'.  Ignored when using FMINSEARCH.
%      'DerivStep'    The relative difference used in finite difference
%                     derivative approximations when using FMINCON, and
%                     'GradObj' is 'off'.  May be a scalar, or the same
%                     size as START.  EPS^(1/3) by default.  Ignored when
%                     using FMINSEARCH.
%      'FunValCheck'  'on' or 'off', indicating whether or not MLE should
%                     check the values returned by the custom distribution
%                     functions for validity.  Default is 'on'.  A poor
%                     choice of starting point can sometimes cause these
%                     functions to return NaNs, infinite values, or out of
%                     range values if they are written without suitable
%                     error-checking.
%       'TolBnd'      An offset for upper and lower bounds when using
%                     FMINCON.  MLE treats upper and lower bounds as
%                     strict inequalities (i.e., open bounds).  With
%                     FMINCON, this is approximated by creating closed
%                     bounds inset from the specified upper and lower
%                     bounds by TolBnd.  Default is 1e-6.
%
%   See also FITDIST, MLECOV, STATSET.

%   [...] = MLE(..., 'optimOptions', OPTS) specifies additional control
%   parameters for FMINSERCH or FMINCON.  OPTS is a structure created by
%   OPTIMSET.  This is applicable only when fitting custom distributions.
%   Choices of OPTIMSET parameters depend on the optimization function.
%   See OPTIMSET('fminsearch') or OPTIMSET('fmincon').

%   When you supply distribution functions, MLE computes the parameter
%   estimates using an iterative maximization algorithm.  With some models
%   and data, a poor choice of starting point can cause MLE to converge to
%   a local optimum that is not the global maximizer, or to fail to
%   converge entirely. Even in cases for which the log-likelihood is
%   well-behaved near the global maximum, the choice of starting point is
%   often crucial to convergence of the algorithm.  In particular, if the
%   initial parameter values are far from the MLEs, underflow in the
%   distribution functions can lead to infinite log-likelihoods.

%   Copyright 1993-2017 The MathWorks, Inc.

% Check for the grandfathered syntax mle('distname',data,alpha,ntrials)
if nargin > 0
    data = convertStringsToChars(data);
end

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if ischar(data)
    if nargin < 2
        error(message('stats:mle:TooFewInputs'));
    elseif nargin <= 4
        dist = data;
        data = varargin{1};
        if nargin < 3
            alpha = [];
        else
            alpha = varargin{2};
        end
        if nargin < 4
            ntrials = [];
        else
            ntrials = varargin{3};
        end
        mu = [];
        theta = [];
        cens = [];
        freq = [];
        options = [];
    else
        error(message('stats:mle:InvalidSyntax'));
    end

else
    pnames = {'distribution' 'censoring' 'frequency' 'alpha' 'ntrials' ...
              'theta' 'mu' 'options' 'pdf' 'logpdf' 'nloglf' 'cdf' 'logsf'};
    dflts =  {[] [] [] [] [] [] [] [] [] [] []};
    % Ask for otherArgs here even though we don't use them, to allow for mlecustom's parameters
    [dist,cens,freq,alpha,ntrials,theta,mu,options,...
        pdf,logpdf,nloglf,cdfFcn,logsfFcn,~,~] = internal.stats.parseArgs(pnames, dflts, varargin{:});
end

if ~isvector(data)
    error(message('stats:mle:InvalidData'));
end

% Check sizes of censoring and frequency vectors.
if ~isempty(cens) && ~isequal(size(data), size(cens))
    error(message('stats:mle:CensSizeMismatch'));
elseif ~isempty(freq) && ~isequal(size(data), size(freq))
    error(message('stats:mle:FreqSizeMismatch'));
end

data = data(:);
cens = cens(:);
freq = freq(:);

if isempty(dist)
    % If a built-in distribution was not given, and a user-defined
    % distribution was given, then handle that.
    if ~(isempty(pdf) && isempty(logpdf) && isempty(nloglf))
        if nargout < 2
            phat = mlecustom(data,varargin{:});
        else
            [phat, pci] = mlecustom(data,varargin{:});
        end
        return

    % Otherwise default to a normal distribution.
    else
        if ~(isempty(cdfFcn) && isempty(logsfFcn))
            % In case someone provided a cdf expecting that to define a
            % distribution, warn about that
            warning(message('stats:mle:MissingPdf'));
        end
        dist = 'normal';
    end

elseif ischar(dist)
    % Remove any zero frequencies.
    if ~isempty(freq)
        zerowgts = find(freq == 0);
        if ~isempty(zerowgts)
            data(zerowgts) = [];
            if ~isempty(cens), cens(zerowgts) = []; end
            freq(zerowgts) = [];
        end
    end
    
    distNames = {'beta', 'bernoulli', 'binomial', 'exponential', 'extreme value', ...
                 'gamma', 'generalized extreme value', 'generalized pareto', ...
                 'geometric', 'half normal', 'lognormal', 'normal', 'negative binomial', ...
                 'poisson', 'rayleigh', 'discrete uniform', 'uniform', 'weibull'};

    i = find(strncmpi(dist, distNames, length(dist)));
    if numel(i) > 1
        error(message('stats:mle:AmbiguousDistName', dist));
    elseif numel(i) == 1
        dist = distNames{i};
    else % it may be an abbreviation that doesn't partially match the name
        dist = lower(dist);
    end
else
    error(message('stats:mle:InvalidDistName'));
end

if isempty(alpha), alpha = 0.05; end

if ~isempty(ntrials)
    if ~strncmpi(dist, 'binomial', length(dist))
        warning(message('stats:mle:NtrialsNotNeeded'));
    end
end

if ~isempty(mu)
    if ~any(ismember({'hn','half normal','halfnormal'},dist))
        warning(message('stats:mle:MuNotNeeded'));
    elseif any(data<mu)
        error(message('stats:mle:BadMu'));
    end
    data = data - mu;
elseif any(ismember({'hn','half normal','halfnormal'},dist)) && (any(data(:)<0))
    error(message('stats:mle:MuIsNeeded'));
end

if ~isempty(theta)
    if ~any(ismember({'gp','generalized pareto','generalizedpareto'},dist))
        warning(message('stats:mle:ThetaNotNeeded'));
    elseif any(data<theta)
        error(message('stats:mle:BadTheta'));
    end
    data = data - theta;
elseif any(ismember({'gp','generalized pareto','generalizedpareto'},dist)) && (any(data(:)<=0))
    error(message('stats:mle:ThetaIsNeeded'));
end

switch dist

case 'bernoulli'
    if ~isempty(cens)
        error(message('stats:mle:CensoringNotSupported', 'Bernoulli'));
    elseif any(data ~= 0 & data ~= 1)
        error(message('stats:mle:InvalidBernoulliData'));
    end
    if ~isempty(freq), data = expandInput(data,freq); end
    if nargout < 2
        phat = binofit(sum(data),numel(data));
    else
        [phat,pci] = binofit(sum(data),numel(data),alpha);
        pci = pci';
    end

case 'binomial'
    if isempty(ntrials)
        error(message('stats:mle:NtrialsNeeded'));
    elseif ~isempty(cens)
        error(message('stats:mle:CensoringNotSupported', 'binomial'));
    end
    ntrials = ntrials(:);
    if ~isempty(freq)
        data = expandInput(data,freq);
        if numel(ntrials) ~= 1
            ntrials = expandInput(ntrials,freq);
        end
    end
    n = numel(data);
    data = sum(data);
    if numel(ntrials) == 1
        ntrials = n .* ntrials;
    else
        ntrials = sum(ntrials);
    end
    if nargout < 2
        phat = binofit(data,ntrials);
    else
        [phat,pci] = binofit(data,ntrials,alpha);
        pci = pci';
    end

case 'geometric'
    if ~isempty(cens)
        error(message('stats:mle:CensoringNotSupported', 'geometric'));
    end
    if ~isempty(freq), data = expandInput(data,freq); end
    phat = 1 ./ (1 + mean(data));
    if nargout > 1
        n = numel(data);
        se = phat .* sqrt((1-phat) ./ n);
        pci = norminv([alpha/2; 1-alpha/2], [phat; phat], [se; se]);
    end

case {'discrete uniform', 'unid'}
    if ~isempty(cens)
        error(message('stats:mle:CensoringNotSupported', 'discrete uniform'));
    end
    if ~isempty(freq), data = expandInput(data,freq); end
    phat = max(data);
    if nargout > 1
        n = numel(data);
        pci = [phat; ceil(phat./alpha.^(1./n))];
    end

case 'uniform'
    if ~isempty(cens)
        error(message('stats:mle:CensoringNotSupported', 'uniform'));
    end
    if ~isempty(freq), data = expandInput(data,freq); end
    if nargout < 2
        [ahat,bhat] = unifit(data);
    else
        [ahat,bhat,aci,bci] = unifit(data,alpha);
        pci = [aci bci];
    end
    phat = [ahat bhat];

otherwise
    spec = dfgetdistributions(dist);
    if ~isscalar(spec)
         error(message('stats:mle:UnrecognizedDistName', dist));
    end
    if ~spec.censoring
        if ~isempty(cens) && ~all(~cens)
            error(message('stats:mle:CensoringNotSupported', dist));
        end
        if ~isempty(freq)
            data = expandInput(data,freq);
            freq = [];
        end
    end
    try
        pd = fitdist(data,dist,'cens',cens,'freq',freq,'options',options);
    catch myException
        switch(myException.identifier)
            case 'stats:ProbDistUnivParam:checkdistname:UnrecognizedName'
                error(message('stats:mle:UnrecognizedDistName', dist));
            otherwise
                rethrow(myException);
        end
    end
    
    % mle is more restrictive than fitdist because it recognizes only
    % parametric distributions, not others such as kernel distributions
    if isa(pd,'prob.ProbabilityDistribution') && ...
      ~isa(pd,'prob.ParametricDistribution')
        error(message('stats:mle:UnrecognizedDistName', dist));
    end

    phat = pd.Params(~pd.ParamIsFixed);
    if nargout>=2
        if isa(pd,'prob.FittableParametricDistribution')
            pci = paramci(pd,'alpha',alpha);
        else
            pci = paramci(pd,alpha);
        end
        pci = pci(:,~pd.ParamIsFixed);
    end
    % UMVUE to MLE conversions
    switch dist
        case {'lognormal', 'normal'}
            % If there was no censoring, fitdist estimated the UMVUE for sigsq.
            if isempty(cens)
                if ~isempty(freq)
                    n = sum(freq);
                else
                    n = numel(data);
                end
                phat(2) = phat(2).*sqrt((n-1)./n);
            end
    end
end


function expanded = expandInput(input,freq)
%EXPANDDATA Expand out an input vector using element frequencies.
if ~isequal(size(input),size(freq))
    error(message('stats:mle:FreqSizeMismatch'));
end
if ~all(freq==round(freq)) || any(freq<0)
    error(message('stats:mle:BadFrequency'));
end
i = cumsum(freq);
j = zeros(1, i(end));
j(i(1:end-1)+1) = 1;
j(1) = 1;
expanded = input(cumsum(j));
