function [H, pVal, ADStat, CV] = adtest(x, varargin)
% ADTEST Anderson-Darling goodness-of-fit hypothesis test.
%   H = ADTEST(X) performs an Anderson-Darling test to determine if a
%   random sample in the vector X could have come from a normal
%   distribution. Missing observations in X, indicated by NaN, are ignored.
%   The output H is the result of the hypothesis test:
%    H = 0 => Do not reject the null hypothesis at significance level 0.05.
%    H = 1 => Reject the null hypothesis at significance level 0.05.
%
%   [H,P] = ADTEST(X) also returns the P-value P.
%
%   [H,P,ADSTAT] = ADTEST(X) also returns the Anderson-Darling test
%   statistic ADSTAT.
%
%   [H,P,ADSTAT,CV] = ADTEST(X) returns the critical value of the test CV,
%   the cutoff value CV at the significance level alpha. The default for
%   alpha is 0.05.
%
%   [...] = ADTEST(X, 'PARAM1',val1, 'PARAM2',val2, ...) specifies one or
%   more of the following parameter name/value pairs to specify whether to
%   perform a simple or composite hypothesis test, to specify the
%   distribution being tested for, to control for the significance level,
%   and to specify whether to perform the test using Monte-Carlo
%   simulations:
%
%   'Distribution'  Distribution being tested for. Tests whether X could
%                   have come from the distribution specified by this
%                   parameter. Choices are:
%      - A ProbabilityDistribution object. In this case, all parameters of
%         the null distribution are specified, and X is tested against a
%         simple hypothesis.
%
%      - One of the strings: 'norm', 'exp', 'ev', 'logn', 'weibull'.
%        In this case, X is tested against a composite hypothesis for the
%        specified distribution family. (Default: 'norm').
%
%   'Alpha'         Significance level alpha for the test. Any numeric
%                   value between 0 and 1. The default is 0.05
%                   corresponding to the 5% significance level.
%
%   'MCTol'         Monte-Carlo tolerance value. In this case, an
%                   approximation for the P value is computed directly,
%                   using Monte-Carlo simulations.
%
%   'Asymptotic'    If true, use the limiting distribution of the
%                   Anderson-Darling statistic. An estimate computed by
%                   this method is likely to be more accurate than the
%                   low-statistics approximation for sample sizes greater
%                   than 120. This option is only valid if you pass a 
%                   ProbabilityDistribution for 'Distribution'. Default: false
%
%   The Anderson-Darling test statistic belongs to the family of Quadratic
%   Empirical Distribution Function statistics, which are based on the
%   weighted sum of the difference [Fn(x)-F(x)]^2 over the ordered sample
%   values  X1 < X2 < ... < Xn, where F is the hypothesized continuous
%   distribution and Fn is the empirical CDF based on the data sample with
%   n sample points.
 
% References:
% D'Agostino and Stephens, Goodness-Of-Fit Techniques, Marcel-Dekker, New
%   York, 1986. 
% Marsaglia, G; Marsaglia JCW; (2004) "Evaluating the Anderson Darling
%   distribution", Journal of Statistical Software, 9(2).
% M.A Stephens (1974) EDF Statistics for Goodness of Fit and Some
%   Comparisons, Journal of the American Statistical Association.
% A.N. Pettitt (1977) Testing the Normality of Several Independent Samples
%   Using the Anderson-Darling Statistic, Journal of the Royal Statistical
%   Society. Series C, 26(2).

% Parse arguments and check if parameter/value pairs are valid
if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

paramNames = {'Distribution', 'Alpha', 'MCTol', 'Asymptotic'};
dflts  =     {'normal' ,       0.05,     [],       false};

[valDistr, valAlpha, valMCTol, valAsymptotic] = ...
    internal.stats.parseArgs(paramNames, dflts, varargin{:});

% The 'Asymptotic' option is not valid for the composite distribution test
if valAsymptotic && ~(isa(valDistr, 'ProbDist') || isa(valDistr, 'prob.ProbabilityDistribution'))
    error(message('stats:adtest:BadAsymptoticWithProbDist'));
    
    % The 'Asymptotic' option is not valid if you want to run Monte-Carlo
    % simulations
elseif valAsymptotic && ~isempty(valMCTol)
    error(message('stats:adtest:BadAsymptoticWithMCTol'));
else
    isAsymptotic = valAsymptotic;
end

% Ensure the sample data is a real vector.
if ~isvector(x) || ~isreal(x)
    error(message('stats:adtest:BadXSample'));
end

% Remove missing values.
x = x(~isnan(x));

% Compute sample size n.
n = length(x);

% Ensure the significance level, ALPHA, is a scalar between 0 and 1.
if ~isscalar(valAlpha) || ~(valAlpha > 0 && valAlpha < 1)
    error(message('stats:adtest:BadAlpha'));
else
    alpha = valAlpha;
end

% Ensure the Monte-Carlo tolerance, MCTOL, is a numeric scalar greater than
% 0.
if ~isempty(valMCTol) && (~isscalar(valMCTol) || valMCTol <=0)
    % Invalid Monte-Carlo tolerance
    error(message('stats:adtest:BadMCTol'));
else
    mctol = valMCTol;
end

% Make sure the distribution specified is a ProbabilityDistribution object or a string
% specifying the distribution being tested for.
if internal.stats.isString(valDistr)
    DistrNames = {'normal', 'exponential', 'ev', 'lognormal', 'weibull', 'extreme value'};
    Distr = internal.stats.getParamVal(valDistr, DistrNames,'distribution name');
    % If the distribution specified is a valid string this is a composite
    % distribution test (i.e. parameters are unknown)
    isCompositeTest = true;
    
elseif isa(valDistr, 'ProbDist') || isa(valDistr, 'prob.ProbabilityDistribution')
    % If the distribution specified is a distribution fit object this is a
    % test against a specific distribution (parameters are completely
    % specified)
    isCompositeTest = false;
    Distr = valDistr;
else
    % Invalid distribution specified
    error(message('stats:adtest:BadDistribution'));
end

% Case 1: simple hypothesis test, which assumes that the parameters of the
% underlying distribution are known.
if ~isCompositeTest
    % Compute the hypothesized CDF at the values in x.
    z = cdf(Distr, x);
    
    % % Compute the Anderson-Darling statistic. If the data comes from a
    % hypothesized continuous distribution, it can be recast to a uniform
    % distribution. Thus the data Zi = F(Xi), 1<=i<=n, are assumed to come from a
    % uniform distruibution in [0,1].
    ADStat = ComputeADStat(z,n);
    
    if isempty(mctol)
        if isAsymptotic
            % Use the asymptotic distribution for the Anderson-Darling
            % test.
            if n<=120
                warning(message('stats:adtest:AsymptoticWithSmallN', n));
            end
            pVal = 1 - ADInf(ADStat);
        else
            if n == 1
                % For n=1 there is an exact analytical formula, valid for ad > log(4)-1,
                % which will always be the case
                pVal = 1 - sqrt(1-4*exp(-1-ADStat));
            else
                if n < 4
                    % You may get a P-value that will not be accurate to
                    % within 4 digits of precision over the entire
                    % distribution range using ADn. (This is mainly the
                    % case around the 33d percentile).
                    warning(message('stats:adtest:SmallSampleSize'));
                end
                % For n>=2 use formulas from Marsaglia
                % Compute the P-value based on this Anderson-Darling statistic
                % P-value = Pr(ADn >= d) = 1-Pr(ADn < d)
                pVal = 1 - ADn(n, ADStat);
            end
        end
        if nargout >  3
            % Compute critical value at specified percentage point
            
            if isAsymptotic
                [alphas, critVals] = findAsymptoticDistributionCriticalValues;
                % Make sure alpha is within the lookup table
                validateAlpha(alpha, alphas);
                i = find(alphas>alpha,1,'first');
                startVal = critVals(i-1);
                CV = fzero(@(ad)1-ADInf(ad)-alpha, startVal);
            else
                [alphas, CVs, sampleSizes] = findCriticalValues;
                % Make sure alpha is within the lookup table
                validateAlpha(alpha, alphas);
                [OneOverSampleSizes, LogAlphas] = meshgrid(1./sampleSizes, log(alphas));
                CV = interp2(OneOverSampleSizes, LogAlphas, CVs', 1./n, log(alpha));
            end
        end
    else
        % Compute the critical value and p-value on the fly using Monte-Carlo simulation.
        [CV, pVal] = adtestMC(ADStat, n, alpha, 'unif', mctol);
    end
    % "H = 0" implies that we "Do not reject the null hypothesis at the
    % significance level alpha," and "H = 1" implies that we "Reject the
    % null hypothesis at significance level alpha."
    H  =  (pVal < alpha);
else  % Composite hypothesis test
    % Case 2: Composite hypothesis test - one or more of the parameters
    % is estimated from the data set.
    
    % For the composite case, if the sample size is less than 4 we error
    % out
    if n<4
        error(message('stats:adtest:NotEnoughData'));
    end
    
    % If data come from a lognormal distribution, log(x) is normally
    % distributed
    if strcmp(Distr, 'lognormal')
        x = log(x);
        Distr = 'normal';
        
        % If data come from a Weibull distribution, log(x) has a type I extreme-value
        % distribution
    elseif strcmp(Distr, 'weibull')
        x = log(x);
        Distr = 'ev';
    end
    
    switch Distr
        case 'normal'
            if any(~isreal(x)) % This can occur when we do the test against the lognormal distribution
                
                % Data is not compatible with this distribution test. in
                % this case, set p-value to zero
                warning(message('stats:adtest:BadDataForDistribution'));
                ADStat = NaN;
            else
                z = normcdf(x, mean(x), std(x));
                ADStat = ComputeADStat(z,n);
            end
        case 'exponential'
            z = expcdf(x, mean(x)); % MLE estimate
            ADStat = ComputeADStat(z,n);
        case {'ev', 'extreme value'}
            if any(~isreal(x)) % This can occur when we do the test against the Weibull distribution
                warning(message('stats:adtest:BadDataForDistribution'));
                ADStat = NaN;
            else
                params = evfit(x);
                z = evcdf(x,params(1),params(2));
                ADStat = ComputeADStat(z,n);
            end
    end
    
    if isempty(mctol)
        switch Distr
            case 'normal', [alphas, CVs] = computeCriticalValues_norm(n);
            case 'exponential', [alphas, CVs] = computeCriticalValues_exp(n);
            case {'ev', 'extreme value'}, [alphas, CVs] = computeCriticalValues_ev(n);
        end
              
        % 1-D interpolation into the tabulated results. In the upper tail,
        % CV vs log(alpha) is approximately linear
        pp = pchip(log(alphas), CVs);
        CV = ppval(pp,log(alpha));
        
        % If alpha is not within the lookup table, throw a warning.
        % Hypothesis result is computed by comparing the p-value with
        % alpha, rather than CV with ADStat
        if alpha < alphas(1)
            CV = CVs(1);
            warning(message('stats:adtest:BadAlpha3'));
        elseif  alpha > alphas(end)
            CV = CVs(end);
            warning(message('stats:adtest:BadAlpha3'));
        end
        
        if ADStat > CVs(1)
            % P value is smaller than smallest tabulated value
            warning(message('stats:adtest:OutOfRangePLow',...
                sprintf('%g', alphas(1))));
            pVal = alphas(1);
        elseif ADStat < CVs(end)
            % P value is larger than largest tabulated value
            warning(message('stats:adtest:OutOfRangePHigh',...
                sprintf('%g', alphas(end))));
            pVal = alphas(end);
        elseif isnan(ADStat)
            % This happens, for example, in the test for exponentiality when there is negative data.
            pVal = 0;
        else
            % Find p-value by inverse interpolation
            i = find(ADStat>CVs,1,'first');
            logPVal = fzero(@(x)ppval(pp,x) - ADStat, log(alphas([i-1,i])));
            pVal = exp(logPVal);
        end
    else
        % Compute the critical value and p-value on the fly using Monte-Carlo simulation.
        [CV, pVal] = adtestMC(ADStat, n, alpha, Distr, mctol);
    end
    if isnan(ADStat)
        H = true;
    else
        if isempty(mctol)
            if (alpha < alphas(1) || alpha > alphas(end))
                H = pVal < alpha;
            else
                H = ADStat > CV;
            end
        else
            H = ADStat > CV;
        end
    end
end

%----------------Subfunctions--------------------------------------------
%------------------------------------------
function ADStat = ComputeADStat(z,n)
% Compute Anderson-Darling Statistic
% Sort the data and compute the statistic
z = reshape(z,n,1);
z = sort(z);
w = 2*(1:n) - 1;
ADStat = -w*(log(z)+ log(1-z(end:-1:1)))/n - n;

%------------------------------------------
function p = ADn(n,ad)
% Anderson-Darling distribution Pr(An<z) for arbitrary n.
x = adinf_short(ad);
p = x + errfix(n,x);

%------------------------------------------
function x = adinf_short(ad)
% Simplified method for computing ADInf(z), as specified in the Marsaglia
% paper. This is the quick and easy ADInf, rather than the full precision
% one. According to Marsaglia, precision beyond ~ 7 digits is likely to be
% wasted as errfix(n,x) only yields around 5-6 digits of accuracy
if any(ad<0)
    error(message('stats:adtest:BadAndersonDarlingStatistic'));
end

x = zeros(size(ad));
adl = ad(ad<2);
x(ad<2) = adl.^(-1/2).*exp(-1.2337./adl).*(2.00012 + (0.247105 - (0.0649821 - ...
    (0.0347962 - (0.0116720 - 0.00168691*adl).*adl).*adl).*adl).*adl);

adh = ad(ad>=2);
x(ad>=2) = exp(-exp(1.0776 - (2.30695 - (0.43424 - (0.082433 - (0.008056 - ...
    0.0003146.*adh).*adh).*adh).*adh).*adh));


%------------------------------------------
function e = errfix(n,x)
% Compute error function as described in Marsaglia's paper. For a given
% random variable Z with distribution ADn, ADn(Z) should be uniformly
% distributed in [0,1). errfix is the function that makes
% ADInf(Z)+errfix(n,ADInf(Z)) uniformly distributed in [0,1).
if any(x<0 | x>1)
    %Error function only defined between 0 and 1
    error(message('stats:adtest:BadValuesForErrorFunction'));
end

e = zeros(size(x));

% Left crossing depends on n
c = 0.01265 + 0.1757/n;

% Define function by intervals using 3 fixed functions g1, g2, g3.
xc1 = x(x<c)/c;
g1 = sqrt(xc1).*(1 - xc1).*(49*xc1 - 102);
e(x<c) = (0.0037/n^3 + 0.00078/n^2 + 0.00006/n)*g1;


xc2 = (x(x >= c & x < 0.8) - c)./(0.8 - c);
g2 = -0.00022633 + (6.54034 - (14.6538 - (14.458 - (8.259 -...
    1.91864.*xc2).*xc2).*xc2).*xc2).*xc2;
e(x >= c & x < 0.8) = (0.04213/n + 0.01365/n^2)*g2;


xc3 = x(x >= 0.8);
e(x >= 0.8) = 1/n*(-130.2137 + (745.2337 - (1705.091 - (1950.646 -...
    (1116.360 - 255.7844.*xc3).*xc3).*xc3).*xc3).*xc3).*xc3;


%------------------------------------------
function ad = ADInf(z)
% Evaluate the Anderson-Darling limit distribution, as presented in the
% Marsaglia paper
if z<0
    % Distribution is invalid for negative values
    error(message('stats:adtest:BadXForAsymptoticDistribution'));
end

% Return 0 below a certain threshold
if z < 0.02
    ad = 0;
    return;
% For an AD statistic z above the following threshold, the AD CDF is equal
% to 1 within floating-point precision, so set the CDF value to 1.
elseif z >=32.4
    ad = 1;
    return;
end

% Use N=500
N = 500;
n = 1:N;
K = 1/z*[1, ((4*n + 1).*cumprod((1/2 - n)./n))];
ADTerms = arrayfun(@(j)ADf(z,j),0:N);
ad = ADTerms*K';

%------------------------------------------
function f = ADf(z,j)
% This function is called by ADInf
% Term for series expansion for f(z,j). Use N=400.

N = 500;

% Compute t=tj=(4j+1)^2*pi^2/(8z)
t = (4*j + 1)^2*1.233700550136170/z;

% First 2 terms in recursive series
% c0=pi*exp(-t)/(sqrt(2t))
% c1=pi*sqrt(pi/2)*erfc(sqrt(t))
c0 = 2.221441469079183*exp(-t)/sqrt(t);
c1 = 3.937402486430604*erfc(sqrt(t));

r = z/8;
f = c0 + c1*r;
% Evaluate the recursion
for n=2:N
    c = 1/(n - 1)*((n - 3/2 - t)*c1 + t*c0);
    r = r*(z/8)*(1/n);
    fn = f + c*r;
    c0 = c1;
    c1 = c;
    if f==fn
        return;
    end
    f = fn;
end


%------------------------------------------
function [alphas, CVs] = computeCriticalValues_norm(n)

alphas = [0.0005 0.001 0.0015 0.002 0.005 0.01 0.025 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5...
    0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99];

% An improved version of the Petitt method for the composite normal case.
% The model used is A_n = A_inf (1+b_0/n+b_1/n^2), where different
% estimates b_0 and b_1 are used for each significance level. This allows
% us to model the entire range of the distribution.
CVs = ...
    [ 1.5649  1.4407  1.3699  1.3187  1.1556    1.0339    0.8733    0.7519    0.6308    0.5598    0.5092    0.4694    0.4366    0.4084...
    0.3835    0.3611    0.3405    0.3212    0.3029    0.2852    0.2679    0.2506    0.2330    0.2144...
    0.1935    0.1673    0.1296] +...
    [-0.9362 -0.9029  -0.8906  -0.8865  -0.8375   -0.7835   -0.6746   -0.5835   -0.4775   -0.4094   -0.3679   -0.3327   -0.3099   -0.2969...
    -0.2795   -0.2623   -0.2464   -0.2325   -0.2164   -0.1994   -0.1784   -0.1569   -0.1377   -0.1201...
    -0.0989   -0.0800   -0.0598]./n +...
    [-8.3249  -6.6022 -5.6461  -4.9685  -3.2208   -2.1647   -1.2460   -0.7803   -0.4627   -0.3672   -0.2833   -0.2349   -0.1442   -0.0229...
    0.0377    0.0817    0.1150    0.1583    0.1801    0.1887    0.1695    0.1513    0.1533    0.1724...
    0.2027    0.3158    0.6431]./n^2;

  
%------------------------------------------
function [alphas, CVs] = computeCriticalValues_exp(n)

alphas = [0.0005 0.001 0.0015 0.002  0.005 0.01 0.025 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5...
    0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99];

% An improved version of the Petitt method for the composite exponential
% case. The model used is A_n = A_inf (1+b_0/n), where different estimates
% b_0 are used for each significance level. This allows us to model the
% entire range of the distribution.
CVs = ...
    [3.2371    2.9303    2.7541    2.6307  2.2454    1.9621    1.5928    1.3223    1.0621    0.9153    0.8134    0.7355    0.6725    0.6194...
    0.5734    0.5326    0.4957    0.4617    0.4301    0.4001    0.3712    0.3428    0.3144    0.2849...
    0.2527    0.2131    0.1581]+...
    [1.6146    0.8716    0.4715    0.2066  -0.4682   -0.7691   -0.7388   -0.5758   -0.4036   -0.3142   -0.2564   -0.2152   -0.1845   -0.1607...
    -0.1409   -0.1239   -0.1084   -0.0942   -0.0807   -0.0674   -0.0537   -0.0401   -0.0261   -0.0116...
     0.0047    0.0275    0.0780]./n;


% ------------------------------------------
function [alphas, CVs] = computeCriticalValues_ev(n)

alphas = [0.0005 0.001 0.0015 0.002  0.005 0.01 0.025 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5...
    0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99];

% An improved version of the Petitt method for the composite extreme value
% case. The model used is A_n = A_inf (1+b_0/n^(1/2) ), where different
% estimates b_0 are used for each significance level. This allows us to
% model the entire range of the distribution.
CVs = ...
    [1.6473    1.5095    1.4301    1.3742  1.1974    1.0667    0.8961    0.7683    0.6416    0.5680    0.5156...
    0.4744    0.4405    0.4115    0.3858    0.3626    0.3415    0.3217...
    0.3029    0.2848    0.2672    0.2496    0.2315    0.2124    0.1909...
    0.1633    0.1223] +...
    [-0.7097   -0.5934   -0.5328   -0.4930  -0.3708  -0.2973   -0.2075   -0.1449   -0.0892   -0.0619   -0.0442...
    -0.0302  -0.0196   -0.0112   -0.0039    0.0024    0.0074    0.0122...
    0.0167   0.0207    0.0245    0.0282    0.0323    0.0371    0.0436...
    0.0549   0.0813]./n.^(1/2);


% ------------------------------------------
function [alphas, CVs, sampleSizes] = findCriticalValues
% Find rows of critical values at relevant significance levels.

alphas =  [0.0005  0.001  0.0015  0.002  0.005  0.01  0.025  0.05  0.10  0.15  0.20  0.25  0.30  0.35  0.40  0.45  0.5  0.55  0.60  0.65  0.70  0.75  0.80  0.85  0.90  0.95  0.99];

%0.0005   0.001     0.0015    0.002     0.005     0.01      0.025    0.05      0.10      0.15      0.20      0.25      0.30      0.35      0.40      0.45       0.50      0.55      0.60      0.65      0.70      0.75      0.80      0.85      0.90      0.95      0.99
CVs = [  7.2943    6.6014    6.1962    5.9088    4.9940    4.3033    3.3946   2.7142    2.0470    1.6682    1.4079    1.2130    1.0596    0.9353    0.8326     0.7465    0.6740    0.6126    0.5606    0.5170    0.4806    0.4508    0.4271    0.4091    NaN       NaN       NaN   ;...  %n=1
    7.6624    6.4955    5.9916    5.6682    4.7338    4.0740    3.2247   2.5920    1.9774    1.6368    1.4078    1.2329    1.0974    0.9873    0.8947     0.8150    0.7448    0.6820    0.6251    0.5727    0.5240    0.4779    0.4337    0.3903    0.3462    0.3030    0.2558;...  %n=2
    7.2278    6.3094    5.8569    5.5557    4.6578    4.0111    3.1763   2.5581    1.9620    1.6314    1.4079    1.2390    1.1065    0.9979    0.9060     0.8264    0.7560    0.6928    0.6350    0.5816    0.5314    0.4835    0.4371    0.3907    0.3424    0.2885    0.2255;...  %n=3
    7.0518    6.2208    5.7904    5.4993    4.6187    3.9788    3.1518   2.5414    1.9545    1.6288    1.4080    1.2416    1.1104    1.0025    0.9110     0.8315    0.7611    0.6977    0.6397    0.5859    0.5352    0.4868    0.4395    0.3920    0.3421    0.2845    0.2146;...  %n=4
    6.9550    6.1688    5.7507    5.4653    4.5949    3.9591    3.1370   2.5314    1.9501    1.6272    1.4080    1.2430    1.1126    1.0051    0.9138     0.8344    0.7640    0.7005    0.6424    0.5884    0.5375    0.4888    0.4411    0.3930    0.3424    0.2833    0.2097;...  %n=5
    6.8935    6.1345    5.7242    5.4426    4.5789    3.9459    3.1271   2.5248    1.9472    1.6262    1.4081    1.2439    1.1140    1.0067    0.9156     0.8362    0.7658    0.7023    0.6441    0.5901    0.5391    0.4901    0.4422    0.3938    0.3427    0.2828    0.2071;...  %n=6
    6.8509    6.1102    5.7053    5.4264    4.5674    3.9364    3.1201   2.5201    1.9451    1.6255    1.4081    1.2445    1.1149    1.0079    0.9168     0.8375    0.7671    0.7036    0.6454    0.5912    0.5401    0.4911    0.4430    0.3944    0.3430    0.2826    0.2056;...  %n=7
    6.8196    6.0920    5.6912    5.4142    4.5588    3.9293    3.1148   2.5166    1.9436    1.6249    1.4081    1.2450    1.1156    1.0087    0.9177     0.8384    0.7681    0.7045    0.6463    0.5921    0.5409    0.4918    0.4436    0.3949    0.3433    0.2825    0.2046;...  %n=8
    6.7486    6.0500    5.6582    5.3856    4.5384    3.9124    3.1024   2.5084    1.9400    1.6237    1.4081    1.2460    1.1171    1.0106    0.9197     0.8406    0.7702    0.7066    0.6483    0.5941    0.5428    0.4935    0.4451    0.3961    0.3441    0.2826    0.2029;...  %n=12
    6.7140    6.0292    5.6417    5.3713    4.5281    3.9040    3.0962   2.5044    1.9382    1.6230    1.4081    1.2465    1.1179    1.0115    0.9207     0.8416    0.7712    0.7077    0.6493    0.5950    0.5437    0.4944    0.4459    0.3968    0.3445    0.2827    0.2023;...  %n=16
    6.6801    6.0084    5.6252    5.3569    4.5178    3.8955    3.0900   2.5003    1.9365    1.6224    1.4081    1.2470    1.1186    1.0123    0.9217     0.8426    0.7723    0.7087    0.6503    0.5960    0.5446    0.4952    0.4466    0.3974    0.3450    0.2829    0.2019;...  %n=24
    6.6468    5.9877    5.6087    5.3425    4.5075    3.8869    3.0837   2.4963    1.9347    1.6218    1.4082    1.2474    1.1193    1.0132    0.9226     0.8436    0.7732    0.7097    0.6513    0.5969    0.5455    0.4960    0.4474    0.3980    0.3455    0.2832    0.2016;...  %n=48
    6.6634    5.9980    5.6169    5.3497    4.5127    3.8912    3.0868   2.4983    1.9356    1.6221    1.4081    1.2472    1.1190    1.0128    0.9222     0.8431    0.7728    0.7092    0.6508    0.5965    0.5451    0.4956    0.4470    0.3977    0.3453    0.2830    0.2017;...  %n=32
    6.6385    5.9825    5.6046    5.3389    4.5049    3.8848    3.0822   2.4953    1.9343    1.6217    1.4082    1.2475    1.1195    1.0134    0.9228     0.8438    0.7735    0.7099    0.6516    0.5972    0.5458    0.4962    0.4476    0.3982    0.3456    0.2833    0.2016;...  %n=64
    6.6318    5.9783    5.6012    5.3360    4.5028    3.8830    3.0809   2.4944    1.9339    1.6215    1.4082    1.2476    1.1197    1.0136    0.9230     0.8440    0.7737    0.7101    0.6517    0.5974    0.5459    0.4964    0.4477    0.3984    0.3457    0.2833    0.2015;...  %n=88
    6.6297    5.9770    5.6001    5.3350    4.5021    3.8825    3.0805   2.4942    1.9338    1.6215    1.4082    1.2476    1.1197    1.0136    0.9231     0.8441    0.7738    0.7102    0.6518    0.5974    0.5460    0.4965    0.4478    0.3984    0.3458    0.2834    0.2015;...  %n=100
    6.6262    5.9748    5.5984    5.3335    4.5010    3.8816    3.0798   2.4937    1.9336    1.6214    1.4082    1.2477    1.1198    1.0137    0.9232     0.8442    0.7739    0.7103    0.6519    0.5975    0.5461    0.4966    0.4479    0.3985    0.3458    0.2834    0.2015;...  %n=128
    6.6201    5.9709    5.5953    5.3308    4.4990    3.8800    3.0787   2.4930    1.9333    1.6213    1.4082    1.2478    1.1199    1.0139    0.9234     0.8443    0.7740    0.7105    0.6521    0.5977    0.5463    0.4967    0.4480    0.3986    0.3459    0.2834    0.2015;...  %n=256
    6.6127    5.9694    5.5955    5.3314    4.4982    3.87813   3.0775   2.4924    1.9330    1.6212    1.4082    1.2479    1.1201    1.0140    0.9235     0.8445    0.7742    0.7106    0.6523    0.5979    0.5464    0.4969    0.4481    0.3987    0.3460    0.2835    0.2015];    %n=Inf

sampleSizes = [1 2 3 4 5 6 7 8 12 16 24 32 48 64 88 100 128 256 Inf];

% ------------------------------------------
function [alphas, critVals] = findAsymptoticDistributionCriticalValues
% Return critical values for asymptotic distribution at relevant significance levels.

alphas =  [0.0005  0.001  0.0015  0.002  0.005  0.01  0.025  0.05  0.10  0.15  0.20  0.25  0.30  0.35  0.40  0.45  0.5  0.55  0.60  0.65  0.70  0.75  0.80  0.85  0.90  0.95  0.99];

           %0.0005            0.001             0.0015              0.002             0.005             0.01              0.025             0.05              0.10              0.15              0.20              0.25             0.30              0.35              0.40              0.45              0.50              0.55              0.60              0.65              0.70              0.75              0.80              0.85               0.90              0.95              0.99
critVals = [6.6127034546551   5.9694013422151   5.5954643397078     5.3313658857909   4.4981996466091   3.8781250216054   3.0774641787107   2.4923671600494   1.9329578327416   1.6212385363175   1.4081977005506   1.2478596347253  1.1200136586965   1.0140004020016   0.9235137094902   0.8445069178452   0.7742142410993   0.7106405935247   0.6522701010084   0.5978828157471   0.5464229310982   0.4968804113119   0.4481425895777   0.3987228486242    0.3460480234939   0.2835161264344   0.2014922164166]; %n=Inf

% ------------------------------------------
function validateAlpha(alpha, alphas)
% Make sure alpha is within the lookup table
if alpha< alphas(1) || alpha>alphas(end)
    error(message('stats:adtest:BadAlpha2', sprintf('%g', alphas(1)),...
        sprintf('%g', alphas(end))));
end

%------------------------------------------
function [crit, p] = adtestMC(ADStat, n, alpha, distr, mctol)
%ADTESTMC Simulated critical values and p-values for Anderson-Darling test.
%   [CRIT,P] = ADTESTMC(ADSTAT,N,ALPHA,DISTR,MCTOL) returns the critical
%   value CRIT and p-value P for the Anderson-Darling test of the null
%   hypothesis that data were drawn from a distribution in the family
%   DISTR, for a sample size N and confidence level 100*(1-ALPHA)%.  P is
%   the p-value for the observed value ADSTAT of the Anderson-Darling
%   statistic.  DISTR is 'norm', 'exp', 'ev' 'or 'unif'. ALPHA is a scalar
%   or vector.  ADTESTMC uses Monte-Carlo simulation to approximate CRIT
%   and P, and chooses the number of MC replications, MCREPS, large enough
%   to make the standard error for P, SQRT(P*(1-P)/MCREPS), less than
%   MCTOL.

vartol = mctol^2;
crit = 0;
p = 0;
mcRepsTot = 0;
mcRepsMin = 1000;

while true
    mcRepsOld = mcRepsTot;
    mcReps = ceil(mcRepsMin - mcRepsOld);
    ADstatMC = zeros(mcReps,1);
    
    switch distr
        % Simulate critical values for the normal
        case 'normal'
            mu0 = 0;
            sigma0 = 1;
            for rep = 1:length(ADstatMC)
                x = normrnd(mu0, sigma0,n,1);
                xCDF = sort(x);
                nullCDF = normcdf(xCDF, mean(x), std(x));
                w = 2*(1:n) - 1 ;
                ADstatMC(rep) = -w*(log(nullCDF)+ log(1-nullCDF(end:-1:1)))/n - n;
            end
        case 'exponential'
            beta0 = 1;
            for rep = 1:length(ADstatMC)
                x = exprnd(beta0,n,1);
                xCDF = sort(x);
                nullCDF = expcdf(xCDF, mean(x));
                w = 2*(1:n) - 1 ;
                ADstatMC(rep) = -w*(log(nullCDF)+ log(1-nullCDF(end:-1:1)))/n - n;
            end
        case {'ev', 'extreme value'}
            mu0 = 0;
            sigma0 = 1;
            for rep = 1:length(ADstatMC)
                x = evrnd(mu0,sigma0,n,1);
                pHat = evfit(x); %MLE Estimate
                xCDF = sort(x);
                nullCDF = evcdf(xCDF, pHat(1), pHat(2));
                w = 2*(1:n) - 1 ;
                ADstatMC(rep) = -w*(log(nullCDF)+ log(1-nullCDF(end:-1:1)))/n - n;
            end
        case 'unif'
            for rep = 1:length(ADstatMC)
                z = sort(rand(n,1));
                w = 2*(1:n) - 1 ;
                ADstatMC(rep) = -w*(log(z)+ log(1-z(end:-1:1)))/n - n;
            end
    end
    
    critMC = prctile(ADstatMC, 100*(1-alpha));
    pMC = sum(ADstatMC > ADStat)./mcReps;
    
    mcRepsTot = mcRepsOld + mcReps;
    crit = (mcRepsOld*crit + mcReps*critMC) / mcRepsTot;
    p = (mcRepsOld*p + mcReps*pMC) / mcRepsTot;
    
    % Compute a std err for p, with lower bound (1/N)*(1-1/N)/N when p==0.
    sepsq = max(p*(1-p)/mcRepsTot, 1/mcRepsTot^2);
    if sepsq < vartol
        break
    end
    % Based on the current estimate, find the number of trials needed to
    % make the MC std err less than the specified tolerance.
    mcRepsMin = 1.2 * (mcRepsTot*sepsq)/vartol;
end
