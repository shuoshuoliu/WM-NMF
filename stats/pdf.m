function y = pdf(name,x,varargin)
%PDF Density function for a specified distribution.
%   Y = PDF(NAME,X,A) returns an array of values of the probability density
%   function for the one-parameter probability distribution specified by NAME
%   with parameter values A, evaluated at the values in X.
%
%   Y = PDF(NAME,X,A,B) or Y = PDF(NAME,X,A,B,C) returns values of the
%   probability density function for a two- or three-parameter probability
%   distribution with parameter values A, B (and C).
%
%   The size of Y is the common size of the input arguments.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.  Each
%   element of Y contains the probability density evaluated at the
%   corresponding elements of the inputs.
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
%   PDF is a generic function that accepts a distribution by name. It is
%   faster to use a more specialized function when possible, such as
%   NORMPDF for the normal distribution.
%
%   See also CDF, ICDF, MLE, RANDOM.

%   Copyright 1993-2015 The MathWorks, Inc.

if nargin > 0
    name = convertStringsToChars(name);
end

if nargin<2,
    error(message('stats:pdf:TooFewInputs'));
end

if ischar(name)
    name = internal.stats.getDistributionName(name);
else
    error(message('stats:pdf:BadDistribution'));
end

if nargin<5
    c = 0;
else
    c = varargin{3};
end
if nargin<4
    b = 0;
else
    b = varargin{2};
end
if nargin<3
    a = 0;
else
    a = varargin{1};
end

switch(name)
    case 'beta'
        y = betapdf(x,a,b);
    case 'binomial'
        y = binopdf(x,a,b);
    case 'chi-square'
        y = chi2pdf(x,a);
    case 'exponential'
        y = exppdf(x,a);
    case 'extreme value'
        y = evpdf(x,a,b);
    case 'f'
        y = fpdf(x,a,b);
    case 'gamma'
        y = gampdf(x,a,b);
    case 'generalized extreme value'
        y = gevpdf(x,a,b,c);
    case 'Generalized Pareto'
        y = gppdf(x,a,b,c);
    case 'geometric'
        y = geopdf(x,a);
    case 'hypergeometric'
        y = hygepdf(x,a,b,c);
    case 'lognormal'
        y = lognpdf(x,a,b);
    case 'negative binomial'
        y = nbinpdf(x,a,b);
    case 'noncentral f'
        y = ncfpdf(x,a,b,c);
    case 'noncentral t'
        y = nctpdf(x,a,b);
    case 'noncentral chi-square'
        y = ncx2pdf(x,a,b);
    case 'normal'
        y = normpdf(x,a,b);
    case 'poisson'
        y = poisspdf(x,a);
    case 'rayleigh'
        y = raylpdf(x,a);
    case 't'
        y = tpdf(x,a);
    case 'discrete uniform'
        y = unidpdf(x,a);
    case 'uniform'
        y = unifpdf(x,a,b);
    case 'weibull'
        y = wblpdf(x,a,b);
        
    otherwise
        % Other valid distributions are known to makedist or the registry
        spec = dfgetdistributions(name);
        if length(spec)>1
            error(message('stats:pdf:AmbiguousDistName', name));
        elseif isscalar(spec)
            % Distribution is found in the toolbox, so it will have a pdffunc
            % that we can call directly.
            n = length(spec.prequired);
            y = feval(spec.pdffunc,x,varargin{1:min(n,end)});
        else
            % makedist will error if
            % (1) invalid distname
            % (2) pdf has been called with non-scalar parameters, the
            %     distribution class does not support this, and it
            %     issues an error of type 'stats:probdists:ScalarParameter'
            % (3) the constructor for the distribution object fails in
            %     some way other than (2)
            % In the future, we may process (2) to explain it in the
            % special context of this function, but for the present, as
            % in (1) and (3), we let the error ride.
            pd = makedist(name,varargin{:});
            y = pdf(pd,x);
        end
end
