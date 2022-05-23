function r = random(name,varargin)
%RANDOM Generate random arrays from a specified distribution.
%   R = RANDOM(NAME,A) returns an array of random numbers chosen from the
%   one-parameter probability distribution specified by NAME with parameter
%   values A.
%
%   R = RANDOM(NAME,A,B) or R = RANDOM(NAME,A,B,C) returns an array of random
%   numbers chosen from a two- or three-parameter probability distribution
%   with parameter values A, B (and C).
%
%   The size of R is the common size of the input arguments.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   R = RANDOM(NAME,A,M,N,...) or R = RANDOM(NAME,A,[M,N,...]) returns an
%   M-by-N-by-... array of random numbers for a one-parameter distribution.
%   Similarly, R = RANDOM(NAME,A,B,M,N,...) or R = RANDOM(NAME,A,B,[M,N,...]),
%   and R = RANDOM(NAME,A,B,C,M,N,...) or R = RANDOM(NAME,A,B,C,[M,N,...]),
%   return an M-by-N-by-... array of random numbers for a two- or
%   three-parameter distribution.
%
%
%   NAME can be:
%
%      'beta'  or 'Beta',
%      'bino'  or 'Binomial',
%      'burr'  or 'Burr',
%      'chi2'  or 'Chisquare',
%      'exp'   or 'Exponential',
%      'ev'    or 'Extreme Value',
%      'f'     or 'F',
%      'gam'   or 'Gamma',
%      'gev'   or 'Generalized Extreme Value',
%      'gp'    or 'Generalized Pareto',
%      'geo'   or 'Geometric',
%      'hn'    or 'Half Normal',
%      'hyge'  or 'Hypergeometric',
%      'logn'  or 'Lognormal',
%      'nbin'  or 'Negative Binomial',
%      'ncf'   or 'Noncentral F',
%      'nct'   or 'Noncentral t',
%      'ncx2'  or 'Noncentral Chi-square',
%      'norm'  or 'Normal',
%      'poiss' or 'Poisson',
%      'rayl'  or 'Rayleigh',
%      'stable'or 'Stable',
%      't'     or 'T',
%      'unif'  or 'Uniform',
%      'unid'  or 'Discrete Uniform',
%      'wbl'   or 'Weibull'.
%
%   Partial matches are allowed and case is ignored.
%
%   RANDOM is a generic function that accepts a distribution by name. It is
%   faster to use a more specialized function when possible, such as RANDN
%   or NORMRND for the normal distribution.
%
%   See also CDF, ICDF, MLE, PDF.

%   Copyright 1993-2018 The MathWorks, Inc.

if nargin > 0
    name = convertStringsToChars(name);
end

if ischar(name)
    name = internal.stats.getDistributionName(name);
else
    error(message('stats:random:IllegalDistribution'));
end

switch(nargin) % check for valid inputs
    case 1
        % nothing to check
    case 2
        if ~isnumeric(varargin{1})
            error(message('stats:random:NonNumeric'))
        end
    case 3
        if ~(isnumeric(varargin{1}) && isnumeric(varargin{2}))
            error(message('stats:random:NonNumeric'))
        end
    case 4
        if ~(isnumeric(varargin{1}) && isnumeric(varargin{2}) && isnumeric(varargin{3}))
            error(message('stats:random:NonNumeric'))
        end
    otherwise
        if nargin>1 && ~all(cellfun(@(x)isnumeric(x),varargin))
            % somewhat slower, so the first three cases avoid cellfun calls
            error(message('stats:random:NonNumeric'))
        end
end

% Determine, and call, the appropriate subroutine
switch name
    % The distributions following here could be processed generically
    % below (from the distribution registry), but they are included here
    % for faster performance. Also, the two most popular distributions are
    % listed first for faster performance.
    case 'normal'
        r = normrnd(varargin{:});
    case 'uniform'
        r = unifrnd(varargin{:});
    case 'beta'
        r = betarnd(varargin{:});
    case 'binomial'
        r = binornd(varargin{:});
    case 'extreme value'
        r = evrnd(varargin{:});
    case 'exponential'
        r = exprnd(varargin{:});
    case 'gamma'
        r = gamrnd(varargin{:});
    case 'generalized extreme value'
        r = gevrnd(varargin{:});
    case 'generalized pareto'
        r = gprnd(varargin{:});
    case 'lognormal'
        r = lognrnd(varargin{:});
    case 'negative binomial'
        r = nbinrnd(varargin{:});
    case 'poisson'
        r = poissrnd(varargin{:});
    case 'rayleigh'
        r = raylrnd(varargin{:});
    case 'weibull'
        r = wblrnd(varargin{:});

    % The distributions below are not in the registry and must be processed
    % here
    case 'chi-square'
        r = chi2rnd(varargin{:});
    case 'f'
        r = frnd(varargin{:});
    case 'geometric'
        r = geornd(varargin{:});
    case 'hypergeometric'
        r = hygernd(varargin{:});
    case 'noncentral f'
        r = ncfrnd(varargin{:});
    case 'noncentral t'
        r = nctrnd(varargin{:});
    case 'noncentral chi-square'
        r = ncx2rnd(varargin{:});
    case 't'
        r = trnd(varargin{:});
    case 'discrete uniform'
        r = unidrnd(varargin{:});

    % For other distribution names, try to find them in the registry
    otherwise
        spec = dfgetdistributions(name);
        if isempty(spec)
            try
                pd = makedist(name);
            catch ME
                % makedist will error if invalid distname: catch and rethrow
                % equivalent error.  If makedist fails for another reason,
                % pass along the error as-is.
                if strcmp(ME.identifier,'stats:ProbDistUnivParam:checkdistname:UnrecognizedName')
                    error(message('stats:random:BadDistribution', name));
                else
                    rethrow(ME)
                end
            end
            
            if ~isa(pd,'prob.ToolboxParametricDistribution')
                error(message('stats:random:BadDistribution', name));
            end
            
            % Split the input arguments into parameter values and size specifications
            nparam = pd.NumParameters;
            params = varargin(1:nparam);
            sizeargs = varargin(nparam+1:end);
            [err,~] = internal.stats.statsizechk(nparam,params{:},sizeargs{:});
            if err > 0
                error(message('stats:normrnd:InputSizeMismatch'));
            end
            r = pd.randfunc(params{:},sizeargs{:});
            return
        elseif length(spec)>1
            error(message('stats:random:AmbiguousDistribution', name));
        end
        
        if isempty(spec.randfunc)
            % Compute by inverting the cdf if necessary
            paramArgs = varargin(1:length(spec.pnames));
            sizeArgs = varargin(length(spec.pnames)+1:end);
            u = rand(sizeArgs{:});
            r = feval(spec.invfunc,u,paramArgs{:});
        else
            r = feval(spec.randfunc,varargin{:});
        end
end
