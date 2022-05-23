function [Sig, Mu, Mahal, Outliers, rsltStruct] = robustcov(X, varargin)
% ROBUSTCOV Compute a robust multivariate covariance and mean estimate.
%   SIG = ROBUSTCOV(X) computes the robust covariance matrix estimate of a
%   multivariate data set X, where X is a N-by-P matrix where each row is an
%   observation and each column a variable.  Supported algorithms are the
%   FASTMCD algorithm by Rousseuw and Leroy, the OGK algorithm by Maronna
%   and Zamar and the family of "concentration algorithms" by Olive and
%   Hawkins. For all supported methods, a reweighting for efficiency step
%   is performed. This step does not affect the robustness but improves
%   the efficiency of the estimator.
%
%   [SIG, MU] = ROBUSTCOV(X) also returns an estimate of the robust MCD mean.
%
%   [SIG, MU, MAH] = ROBUSTCOV(X) also returns the robust distances, computed
%   as the Mahalanobis distances of the observations using the robust
%   estimates of mean and covariance.
%
%   [SIG, MU, MAH, OUTLIERS] = ROBUSTCOV(X) also returns the indices of the
%   observations flagged as outliers in the data set.
%
%   [SIG, MU, MAH, OUTLIERS, S] = ROBUSTCOV(X) also returns a struct
%   containing information about the estimate.
%
%   [...] = ROBUSTCOV(X, 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to control the method being used
%   and the handling of special cases. Parameters are:
%
%    'Method'    - A string specifying which estimator to use. Choices are:
%           * 'fmcd'          - Use the FAST-MCD (Minimum Covariance Determinant)
%                               method by Rousseuw and Leroy.  This method
%                               looks for h observations out of N (where
%                               N/2 < h <= N) whose classical covariance
%                               matrix has the lowest possible determinant.
%                               The estimate is then the covariance matrix
%                               of the h points defined above, multiplied
%                               by a consistency factor to obtain
%                               consistency at the multivariate normal
%                               distribution, and by a correction factor to
%                               correct for bias at small samples. (This is
%                               the default.)
%
%          *  'ogk'           - Use the Orthogonalized Gnanadesikan-Kettenring (OGK)
%                               estimate by Maronna and Zamar. This
%                               estimate is a positive definite estimate of
%                               scatter starting from the Gnanadesikan and
%                               Kettering(GK) estimator, a pairwise robust
%                               scatter matrix that may be non-positive
%                               definite. The estimate uses a form of
%                               Principal Components, called an
%                               orthogonalization iteration, on the
%                               pairwise scatter matrix, replacing its
%                               eigenvalues, which could be negative, by
%                               robust variances. This procedure can be
%                               iterated for improved results, and
%                               convergence is usually obtained after 2 or
%                               3 iterations.
%
%        *  'olivehawkins'    - Use the "concentration algorithm" techniques
%                               proposed by Olive and Hawkins, a family of
%                               fast, consistent and highly
%                               outlier-resistant methods. The estimate is
%                               obtained by first generating trial
%                               estimates, or starts, and then using the
%                               concentration technique from each trial fit
%                               to obtain attractors. By default, two
%                               attractors are used. The first attractor is
%                               the DGK (Devlin-Gnanadesikan-Kettering)
%                               attractor, where the start is the classical
%                               estimator. The second attractor is the
%                               Median Ball (MB) attractor, where the start
%                               used is Mu=median(X) and Sigma=eye(P),
%                               i.e. the half set of data closest to
%                               median(X) in Euclidean distance. The MB
%                               attractor is used if the location estimator
%                               of the DGK attractor is outside of the
%                               median ball, and the attractor with the
%                               smallest determinant is used otherwise. The
%                               final mean estimate is the mean estimate of
%                               the chosen attractor, and the final
%                               covariance estimate is the covariance
%                               estimate of the chosen attractor,
%                               multiplied by a scaling factor to make the
%                               estimate consistent at the normal
%                               distribution.
%
%   For the 'fmcd' and 'olivehawkins' methods:
%     'OutlierFraction'      - A value between 0 (included) and 0.5.
%                              1-outlierfraction specifies the fraction of
%                              observations over which the covariance
%                              determinant is minimized. The default is
%                              0.5, for which the algorithm chooses a
%                              subsample of size h = ceil((N + P + 1)/2),
%                              where n is the number of observations and p
%                              is the number of dimensions. This is the
%                              value for which the maximum possible
%                              breakdown is achieved. This parameter
%                              controls the size of the subsets h over
%                              which the covariance determinant is
%                              minimized. The algorithm then chooses h ~
%                              (1-outlierfraction)*N observations per
%                              subset.
%    'NumTrials'             - For FMCD, this is the number of random subsamples of size
%                              P + 1 that are drawn from the dataset as
%                              starting points in the algorithm. The
%                              default is 500. 
%                              For 'olivehawkins', this is the number of
%                              trial fits, or attractors, to be used. By
%                              default, 'NumTrials' is equal to 2.
%
%   For the 'fmcd' method:
%    'BiasCorrection'        - A logical value specifying whether to use a
%                              small-sample correction factor in the
%                              covariance estimate. Default: true.
%
%   For the 'ogk' method:
%    'NumOGKIterations'      - An integer specifying the number of orthogonalization
%                              iterations. This is usually 1 or 2, and
%                              further steps are unlikely to improve the
%                              estimation. Default: 2.
%
%    'UnivariateEstimator'   - Function to use for computing the univariate robust
%                              location and scale estimates. Choices are:
%                                * 'tauscale' - Use the "tau-scale" estimate of Yohai and
%                                               Zamar.
%                                * 'qn'       - Use the Qn scale estimate 
%                                               of Croux and Rousseuw.
%
%   For the 'olivehawkins' method:
%    'ReweightingMethod'     - A string specifying how the reweighting for
%                              efficiency step is performed. Choices are:
%                                * 'rfch'       - RFCH uses two reweighting steps.
%                                                 This is a standard
%                                                 method of reweighting to
%                                                 improve efficiency. (This
%                                                 is the default.)
%                                *'rmvn'        - RMVN (reweighted multivariate normal) 
%                                                 uses two reweighting
%                                                 steps that can be useful
%                                                 for estimating the true
%                                                 covariance matrix under a
%                                                 variety of outlier
%                                                 configurations when the
%                                                 clean data are
%                                                 multivariate normal.
%
%    'NumConcentrationSteps'  - An integer k specifying the number of
%                               concentration steps to use for each attractor.
%                               Default: 10.
%
%    'StartMethod'            - A string, function handle, or cell array of strings and
%                               function handles specifying the start
%                               method to use for each attractor. The
%                               number of attractors used is equal to the
%                               length of the cell array. This option gives
%                               users more control over the algorithm and
%                               the ability to specify a custom number of
%                               attractors and starts. Choices are:
%                                 * 'classical'      - The start is the classical 
%                                                      estimator. 
%                                 * 'medianball'     - The start is the Median
%                                                       Ball estimator.
%                                 * 'elemental'      - The attractor is generated 
%                                                      by concentration where the
%                                                      start is a randomly
%                                                      selected elemental
%                                                      start: the classical
%                                                      estimator applied to
%                                                      a randomly selected
%                                                      "elemental set" of
%                                                      P + 1 cases.
%                                 * Function handle - Specify a function 
%                                                      that returns two output
%                                                      arguments used for
%                                                      computing the
%                                                      initial location and
%                                                      scatter.

%   Copyright 2015-2019 The MathWorks, Inc.

% References:
% P.J. Rousseuw and K. Van Driessen(1999), A fast algorithm for the minimum
%   covariance determinant estimator. Technometrics, 41. G.
% Pison, S. Van Aelst and G. Willems(2002), Small Sample Corrections for
%   LTS and MCD. Metrika 2002, 55.
% R. Maronna and R.H. Zamar(2002), Robust estimates of location and
%   dispersion for high dimensional datasets, Technometrics, 50, 2002.
% Olive, D.J. (2004), A Resistant Estimator of Multivariate Location and
%   Dispersion, Computational Statistics and Data Analysis, 46, 99-102.
% Olive, D.J., and Hawkins, D.M. (2010.) Robust multivariate location and
%   dispersion. Preprin

% ROBUSTCOV does not support integer, complex or sparse data.
if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

internal.stats.checkSupportedNumeric('X', X);

% ROBUSTCOV does not support ND arrays
if ~ismatrix(X)
    error(message('stats:robustcov:InputDim'));
end

% For univariate data, make X a column
if isvector(X)
    X = X(:);
end

% Remove missing values - Robust covariance algorithms are not valid when
% there are have missing values
NObs = size(X,1); % Number of observations before removing missing values.
rowswithnans = any(isnan(X), 2);
if any(rowswithnans)
    warning(message('stats:robustcov:MissingDataRemoved'));
    X = X(~rowswithnans, :);
end

% n points in p dimensional space
[n, p] = size(X);

if n < 2*p
    % The sample size should be at least twice the number of variables
    error(message('stats:robustcov:SmallSampleSize'));
end

% Parse arguments and check if parameter/value pairs are valid
paramNames = {'Method', 'OutlierFraction', 'NumTrials', 'BiasCorrection',...
    'NumOGKIterations', 'UnivariateEstimator', 'ReweightingMethod',...
    'NumConcentrationSteps', 'StartMethod'};

dflts  = {[], [], [], [], [], [], [], [], []};

[valMethod, valOutlierFraction, valNumTrials, valBiasCorrection, valNumOGKIterations,...
    valUnivariateEstimator, valReweightingMethod, valNumConcentrationSteps, valStartMethod] = ...
    internal.stats.parseArgs(paramNames, dflts, varargin{:});

useFMCDorOH = ~isempty(valOutlierFraction) || ~isempty(valNumTrials);
useFMCD = ~isempty(valBiasCorrection);
useOH = ~isempty(valReweightingMethod) || ~isempty(valNumConcentrationSteps) || ~isempty(valStartMethod);
useOGK = ~isempty(valNumOGKIterations) || ~isempty(valUnivariateEstimator);

% If 'Method' is empty, choose default method based on user-specified
% parameters.
if isempty(valMethod)
    if useFMCDorOH
        Method = 'fmcd';
        if useOH && ~useFMCD
            Method = 'olivehawkins';
        end
    elseif useFMCD
        Method = 'fmcd';
    elseif useOH
        Method = 'olivehawkins';
    elseif useOGK
        Method = 'ogk';
    else
        Method = 'fmcd';
    end
else
    % Validate string for Method input.
    MethodNames = {'fmcd', 'ogk', 'olivehawkins'};
    Method = internal.stats.getParamVal(valMethod, MethodNames,...
        '''Method''');
end

% Validate value for OutlierFraction parameter. This should be a value
% between 0 and 0.5.
if strcmp(Method, 'ogk') && ~isempty(valOutlierFraction)
    error(message('stats:robustcov:OutlierFractionUnsupported'));
else
    if isempty(valOutlierFraction)
        valOutlierFraction = 0.5;
        Alpha = 1 - valOutlierFraction;
    else
        if ~isreal(valOutlierFraction) || ~isscalar(valOutlierFraction) || ~(valOutlierFraction >= 0 && valOutlierFraction <= 0.5)
            error(message('stats:robustcov:BadOutlierFraction'));
        else
            Alpha = 1 - valOutlierFraction;
        end
    end
end

% For FMCD and Olive-Hawkins, compute h, the susbsample size over which
% concentration steps are performed
if ~strcmp(Method, 'ogk')
    n2 = floor((n+p+1)/2);
    h = floor(2*n2 - n + 2*(n-n2)*Alpha);
end

if strcmpi(Method, 'olivehawkins')
    if ~iscell(valStartMethod)
        % Case where start is specified as a single string or function
        % handle, or empty.
        valStartMethod = {valStartMethod};
    end
    emptyStartMethod = all(cellfun(@isempty, valStartMethod));
else
    if ~isempty(valStartMethod)
        error(message('stats:robustcov:StartMethodUnsupported'));
    end
end

% Deafult 'NumTrials' option.
if any(strcmpi(Method, {'fmcd', 'olivehawkins'}))
    if isempty(valNumTrials)
        % Set default values for NumTrials
        switch Method
            case 'fmcd'
                NumTrials = 500;
            case 'olivehawkins'
                if emptyStartMethod
                    % Default start method for OH is {'medianball',
                    % 'classical'}
                    NumTrials = 2;
                else
                    NumTrials = numel(valStartMethod);
                end
        end
    else
        % Validate value for NumTrials parameter.
        if ~isscalar(valNumTrials) || ~isreal(valNumTrials) || ~(valNumTrials > 0)...
                || round(valNumTrials) ~= valNumTrials
            error(message('stats:robustcov:BadNumTrials'));
        else
            NumTrials = valNumTrials;
        end
    end
else % OGK
    if ~isempty(valNumTrials)
        % Outlier Fraction parameter is not supported for OGK.
        error(message('stats:robustcov:NumTrialsUnsupported'));
    end
end

% Validate BiasCorrection parameter
if strcmpi(Method, 'fmcd')
    if isempty(valBiasCorrection)
        valBiasCorrection = true;
    end
    % Validate boolean value for BiasCorrection parameter.
    useBiasCorrection = internal.stats.parseOnOff(valBiasCorrection, '''BiasCorrection''');
else
    if ~isempty(valBiasCorrection)
        error(message('stats:robustcov:BiasCorrectionUnsupported'));
    end
end

% Validate values for 'NumOGKIterations' and 'UnivariateEstimator' parameters.
if strcmpi(Method, 'ogk')
    if isempty(valNumOGKIterations)
        NumOGKIterations = 2;
    elseif ~isscalar(valNumOGKIterations) || ~isreal(valNumOGKIterations)...
            || ~(valNumOGKIterations > 0) || round(valNumOGKIterations) ~= valNumOGKIterations
        error(message('stats:robustcov:BadNumOGKIterations'));
    else
        NumOGKIterations = valNumOGKIterations;
    end
    if ~isempty(valUnivariateEstimator)
          UnivariateEstimatorNames = {'tauscale', 'qn'};
          valUnivariateEstimator = internal.stats.getParamVal(valUnivariateEstimator,...
              UnivariateEstimatorNames, '''UnivariateEstimator''');
    else
        valUnivariateEstimator = 'tauscale';
    end
    robustUnivariateFcn = str2func(valUnivariateEstimator);
else
    if ~isempty(valNumOGKIterations)
        error(message('stats:robustcov:NumOGKIterationsUnsupported'));
    end
    if ~isempty(valUnivariateEstimator)
        error(message('stats:robustcov:UnivariateEstimatorUnsupported'));
    end
end

% Validate values for 'StartMethod', 'ReweightingMethod' and
% 'NumConcentrationSteps' parameters.
StartMethodNames = {'classical', 'medianball', 'elemental'};

if strcmp(Method, 'olivehawkins')
    if isempty(valNumConcentrationSteps)
        k = 10;
    else
        % Validate value for NumConcentrationSteps parameter.
        if ~isscalar(valNumConcentrationSteps) || ~isreal(valNumConcentrationSteps)...
                || ~(valNumConcentrationSteps > 0) || round(valNumConcentrationSteps) ~= valNumConcentrationSteps
            error(message('stats:robustcov:BadNumConcentrationSteps'));
        else
            k = valNumConcentrationSteps;
        end
    end
    
    if isempty(valReweightingMethod)
        reweightMethod = 'rfch';
    else
        % Validate value for ReweightingMethod parameter.
        ReweightNames = {'rfch', 'rmvn'};
        reweightMethod = internal.stats.getParamVal(valReweightingMethod, ReweightNames,...
            '''ReweightingMethod''');
    end
    
    if emptyStartMethod
        start = {'classical', 'medianball'};
        numStarts = 2;
    else
        numStarts = numel(valStartMethod);
        if numStarts > 1
            % If medianball or classical are repeated more than once, remove them.
            for s = {'classical', 'medianball'}
                sInds = find(strcmpi(s, valStartMethod));
                if numel(sInds) > 1
                    valStartMethod(sInds(2:end)) = [];
                end
            end
        end
        start = valStartMethod;
        
        invalidStarts = find(~cellfun(@(x) any(x),...
            cellfun(@(x)any(strcmpi(StartMethodNames, x) | isa(x, 'function_handle')),...
            start, 'uniformoutput',false)), 1);
        
        if ~isempty(invalidStarts)
            error(message('stats:robustcov:InvalidStartMethod'));
        end
    end
    NumAttractors = numel(start);
    
    % Check if one of the chosen starts is 'elemental'
    hasElementalStart = any(strcmpi(start, 'elemental'));
    
    if NumAttractors==1
        % If elemental attractor is used and NTrials > 1, elemental attractor
        % should be repeated.
        if NumTrials > 1
            if ischar(start{:}) && hasElementalStart
                start = repmat(start, 1, NumTrials);
            else
                if NumTrials ~= numStarts
                    error(message('stats:robustcov:NumTrialsNumStartsNotEqual'));
                end
                NumTrials = NumAttractors;
            end
        end
    else
        if NumTrials ~= numStarts
            if hasElementalStart 
               % One (or more) of the starts is elemental, the others are
               % not. Expand elemental start to equal the number of trials
               % NumTrails.
                numElementalStarts = NumTrials-numStarts;
                if numElementalStarts < 0
                    error(message('stats:robustcov:NumTrialsNumStartsNotEqual'));
                end
                start = [start, repmat({'elemental'}, 1, numElementalStarts)];
                NumAttractors = numel(start);
            else
                error(message('stats:robustcov:NumTrialsNumStartsNotEqual'));
            end
        end
        NumTrials = NumAttractors;
    end
else
    % These parameters should error for the other methods
    if ~isempty(valNumConcentrationSteps)
        error(message('stats:robustcov:NumConcentrationStepsUnsupported'));
    end
    if ~isempty(valReweightingMethod)
        error(message('stats:robustcov:ReweightingMethodUnsupported'));
    end
end

switch Method
    case 'olivehawkins'
        % Steps 1 and 2 - For each iteration, compute the initial mean and
        % covariance using the specified start method. Then perform k
        % concentration steps starting from the initial estimate. The final
        % mean and covariance estimates are the attractor for this start.
        
        % For special case h==n we do not perform concentration steps as
        % the estimate is itself.
        isFullSubset = h==n;
        hasSingularSubset = false;
        
        % Compute Median Ball Attractor if it is specified in start. The
        % Median Ball attractor will be used in the criterion for choosing
        % whether to include or exclude every other attractor, so it is
        % computed up front. 
        % If a singular subset is found during a consentration step, its
        % covariance determinant is zero and this will be the final
        % estimate (as no other attractor can possibly have smaller
        % covariance determinant).
        hasMBAttractor = any(strcmp('medianball', start));
        if hasMBAttractor
            % Compute the coordinatewise median of the data, as well as the
            % median value of the median ball distances. These values will
            % be needed for picking the final estimate based on the Median
            % Ball criterion.
            med_X = median(X);
            medianBallDists = sum(bsxfun(@minus, X, med_X).^2, 2); % pdist2(X, med_X, 'mahal', I_p);
            medianMedianBallDists = median(medianBallDists);
            
            [t, V] = computeStartEstimate(X, 'medianball');
            if ~isFullSubset
                for c = 1:k % Do concentration steps
                    if ~hasSingularSubset
                        [t, V, ~, hasSingularSubset] = computeCStep(X, t, V, h);
                    end
                end
            end
            tBest = t;
            VBest = V;
            detMin = det(V);
        else
            detMin = Inf;
        end
        
        % Compute attractors for the other starts.
        for N=1:NumTrials
            if ~hasSingularSubset % No singular subset was found for the Median Ball attractor or using another attractor.
                currentStart = start{N};
                if ~strcmp(currentStart, 'medianball')
                    [t, V] = computeStartEstimate(X, currentStart);
                    if ~isFullSubset
                        for c = 1:k  % Do concentration steps for k iterations.
                            if ~hasSingularSubset
                                [t, V, ~, hasSingularSubset] = computeCStep(X, t, V, h);
                            end
                        end
                    end
                    detV = det(V);
                    if hasMBAttractor
                        % If the Median Ball is one of the attractors, then
                        % check if the Euclidean distance from the current
                        % location estimate to the coordinatewise median is
                        % smaller than the median of the median ball
                        % distances. If not, discard it. 
                        % If a singular subset was found and its mean is
                        % within the median ball distances, this will be
                        % the final estimate.
                        if (t - med_X)*(t - med_X)' <= medianMedianBallDists
                            if detV < detMin
                                tBest = t;
                                VBest = V;
                                detMin = detV;
                            end
                        end
                    else
                        % If the Median Ball attractor is not one of the
                        % attractors, choose the attractor with the smallest
                        % determinant.
                        if detV < detMin
                            tBest = t;
                            VBest = V;
                            detMin = detV;
                        end
                    end
                end
            end
        end
        
        if hasSingularSubset % A singular subset was found. Do not apply reweighting.
            warning(message('stats:robustcov:SingularHSubsetFoundOH'));
            Mu = tBest;
            Sig = VBest;
        else
            
            t_fch = tBest;
            
            % Apply reweighting to make the covariance estimate consistent at
            % the normal model.
            V_fch = median(localSquaredMahal(X, tBest, VBest))/chi2inv(0.5, p)*VBest;
            
            % Now apply reweight for efficiency step.
            switch reweightMethod
                case 'rfch'
                    % The RFCH estimator uses two standard "reweight for
                    % efficiency" steps.
                    
                    % Step 1. Apply classical estimator to n1 points s.t.
                    % Di(T,V)^2<=chi2(p, 0.975)
                    n1_inds = localSquaredMahal(X, t_fch, V_fch) <= chi2inv(0.975, p);
                    [Sig1, Mu1] = localMeanCov(X(n1_inds, :));
                    
                    % Multiply covariance estimate by a consistency factor.
                    Sig1_bar = median(localSquaredMahal(X, Mu1, Sig1))/chi2inv(0.5, p)*Sig1;
                    
                    % Step 2. Apply classical estimator to n2 points s.t.
                    % Di(T,V)^2<=chi2(p, 0.975), where (T,V) is the estimate in
                    % step 1.
                    n2_inds = localSquaredMahal(X, Mu1, Sig1_bar) <= chi2inv(0.975, p);
                    [Sig2_bar, t_RFCH] = localMeanCov(X(n2_inds, :));
                    
                    % Multiply covariance estimate by a consistency factor.
                    V_RFCH = median(localSquaredMahal(X, t_RFCH, Sig2_bar))/chi2inv(0.5, p)*Sig2_bar;
                    
                    % Return final estimates.
                    Mu = t_RFCH;
                    Sig = V_RFCH;
                    
                case 'rmvn'
                    % The RMVN estimator uses two reweighting steps. If the
                    % bulk of the data is normally distributed, then RMVN gives
                    % useful covariance estimates under a variety of outlier
                    % configurations.
                    
                    % Step 1. Apply classical estimator to n1 points s.t.
                    % Di(T,V)^2<=chi2(p, 0.975)
                    n1_inds = localSquaredMahal(X, t_fch, V_fch) <= chi2inv(0.975, p);
                    [Sig1, Mu1] = localMeanCov(X(n1_inds, :));
                    n1 = sum(n1_inds);
                    
                    % Find q1 = min{0.5(0.975)n/n1, 0.995}. Apply reweighting
                    % step for consistency at normal distributions.
                    q1 = min(0.5*(0.975)*n/n1, 0.995);
                    Sig1_bar = median(localSquaredMahal(X, Mu1, Sig1))/chi2inv(q1, p)*Sig1;
                    
                    % Step 2. Let (t_RMVN, Sig2) be the classical estimator
                    % applied to the n2 cases with Di(Mu1, Sig1_bar)^2<=chi2(p, 0.975)
                    n2_inds = localSquaredMahal(X, Mu1, Sig1_bar) <= chi2inv(0.975, p);
                    [Sig2, t_RMVN] = localMeanCov(X(n2_inds, :));
                    n2 = sum(n2_inds);
                    
                    % Let q2 = min{0.5(0.975)n/n2, 0.995}. Apply new
                    % reweighting step to obtain final estimate.
                    q2 = min(0.5*(0.975)*n/n2, 0.995);
                    V_RMVN = median(localSquaredMahal(X, t_RMVN, Sig2))/chi2inv(q2, p)*Sig2;
                    
                    % Return final estimates.
                    Mu = t_RMVN;
                    Sig = V_RMVN;
            end
        end
        if nargout > 4
            rsltStruct.OutlierFraction = valOutlierFraction;
            rsltStruct.reweightingMethod = reweightMethod;
            rsltStruct.NumTrials = NumTrials;
        end
        
    case 'ogk'
        % OGK method is invalid for univariate data
        if p==1
            error(message('stats:robustcov:UnivariateDataOGK'));
        end
        
        % Initialize
        Z = X;
        A = eye(p);
        
        for N=1:NumOGKIterations
            % Step 1. Make the estimate scale invariant using the robust
            % univariate scale estimate on all predictors.
            Sig = robustUnivariateFcn(Z);
            D = diag(Sig);
            Y = bsxfun(@rdivide, Z, Sig);%Z*diag(1./Sig');
            
            % Extend definition to include zero scales. if sig(Xj) = 0 in
            % the above step, set Yj=0.
            Y(:, Sig==0) = 0;
            
            % Step 2. Compute the "Correlation matrix" , using a robust
            % estimate of the covariance v(Yj, Yk) of 2 random variables Yj
            % and Yk. Here the GK estimate is used. Estimate the covariance
            % of random variables X and Y as:
            % v(X,Y)=1/4(sig(X+Y)^2-sig(X-Y)^2).
            U = zeros(p);
            for i=1:p
                Yi = Y(:,i);
                for j=i+1:p
                    Yj = Y(:,j);
                    Sig_ipj = robustUnivariateFcn(Yi+Yj);
                    Sig_imj = robustUnivariateFcn(Yi-Yj);
                    U(i,j) = 1/4*(Sig_ipj^2 - Sig_imj^2);
                end
            end
            % Compute the correlation matrix such that U(i,j) = U(j,i) and
            % U(i,i)=1, for all i=1...p
            U = U+U'+eye(p);
            
            % Step 3. Compute the eigenvalue decoomposition of U. E is the
            % matrix of eigenvectors such that U=E*L*E', where L is the
            % diagonal matrix of eigenvalues
            [E, ~] = eig(U);
            
            % Step 4. Project Y onto the space of eigenvectors E.
            Z = Y*E;
            A = A*D*E;
        end
        
        % Do a robust location and scale estimate on the (approximately
        % uncorrelated) Zj's
        [V_Z, t_Z] = robustUnivariateFcn(Z);
        
        % Apply a reweighted version of OGK for improved estimate.
        
        % 1. Compute Mahalanobis distances d(Xi)=mahal(Xi,t,V) for each
        % column of X, where t and V are the basic OGK estimates of the
        % mean and covariance.
        d_mahal = sum(bsxfun(@rdivide, bsxfun(@minus, Z, t_Z), V_Z).^2 , 2);
        
        % If estimate includes zero scales, distances are zero.
        d_mahal(:, any(V_Z == 0)) = 0;

        % 2. Compute hard rejection threshold used in the hard rejection
        % weight function W(d) = I(d<=d0) used for discarding observations
        % whose Mahalanobis distances are above this value.
        
        % Specify B, where Chi2(B,p) the B-quantile of the chi-square
        % distribution with p degrees of freedom (simulations show best
        % results for B = 0.9).
        B = 0.9;
        d0 = chi2inv(B,p)*median(d_mahal)/chi2inv(0.5,p);
        
        % 3. Weighted mean and covariance
        hardLimIdx = d_mahal<=d0;
        t_w = sum(X(hardLimIdx, :))./ sum(hardLimIdx);
        Xc = bsxfun(@minus, X(hardLimIdx,:), t_w);
        V_w = (Xc'*Xc)/sum(hardLimIdx);
        
        % Robust location and scale estimates are the reweighting step
        % estimates.
        Mu = t_w;
        Sig = V_w;
        
        if nargout > 4
            rsltStruct.NumOGKIterations = NumOGKIterations;
            rsltStruct.UnivariateScale = func2str(robustUnivariateFcn);
        end
        
    case 'fmcd'
        % Compute coordinatewise median and median absolute deviation.
        MEDd = median(X);
        MADd = mad(X, 1);
        
        tol = 10*eps(max(MADd));
        
        if any(MADd < tol) % zero h-order statistic for one or more variables
            % In this case, the data is certain to have an h-subset for
            % which the covariance is singular. Do not carry out the
            % algorithm and return covariance of singular h-subset.
            warning(message('stats:robustcov:ZeroHOrderStatistic'));
            firstSingInds = find(MADd < tol, 1, 'first');
            ObsWithZeroMADInds = abs(X(:,firstSingInds) - MEDd(firstSingInds)) < tol;
           
            [Sig, Mu] = localMeanCov(X(ObsWithZeroMADInds, :));
            
            Mahal = zeros(NObs,1);  % Observations with missing values are set to NaN.
            if ~ZeroCovDet(Sig) % Check for singularity of the covariance.
                Mahal(~rowswithnans, :) = sqrt(localSquaredMahal(X, Mu, Sig));
                Mahal(rowswithnans, :) = NaN;
            end
            Outliers = false(NObs, 1);
            Outliers(~rowswithnans, :) = Mahal > sqrt(chi2inv(0.975, p));
            
            if nargout > 4
                % Return results in a struct as the fifth output argument.
                rsltStruct.Mu = Mu;
                rsltStruct.Sigma = Sig;
                rsltStruct.Method = Method;
                rsltStruct.Distances = Mahal;
                rsltStruct.Outliers = Outliers;
                rsltStruct.BiasCorrection = useBiasCorrection;
                rsltStruct.NumTrials = NumTrials;
                rsltStruct.OutlierFraction = valOutlierFraction;
            end
            return;
        end
                
        % Standardize the data using univariate robust estimates for improved
        % estimate (i.e. use the median and median absolute deviation). We
        % can do this as the MCD estimate is affine equivariant.
        Xs = bsxfun(@rdivide, bsxfun(@minus, X, MEDd), MADd);
        
        % Special case 1: univariate data.
        if p==1 && h~=n
            % For the univariate case, compute the MCD estimate by the
            % exact algorithm of Rousseuw and Leroy (1987, pp171-172).
            
            % First, sort the data.
            Xs = Xs(:);
            x = sort(Xs);
            
            hM = zeros(n-h+1,1); % Initialize h*means vector
            SQ = zeros(n-h+1,1); % Initialize sum of squares (or h*variances) vetors
            
            % Initialze values for first subset.
            hM1 = sum(x(1:h));
            hM(1) = hM1;
            SQ(1) = sum(x(1:h).^2) - hM1^2/h;
            
            % Update sum of squares values SQ recursively
            for j=2:n-h+1
                hM(j) = hM(j-1) - x(j-1) + x(j+h-1);
                SQ(j) = SQ(j-1) - x(j-1)^2 + x(j+h-1)^2 + hM(j-1)^2/h - hM(j)^2/h;
            end
            
            % Find subset(s) for which SQ (i.e subset variances) is
            % minimized.
            [~, hminIdx] = min(SQ, [], 1);
            
            t_MCD = hM(hminIdx)/h;
            V_full = SQ(hminIdx)/(h-1);
            
            % Apply consistency factor
            c0 = computeConsistencyFactor(h, n, p);
            V_MCD = c0*V_full;
            
            % If we use small-sample corrections, multiply by a correction
            % factor to make the estimate unbiased for small samples.
            if useBiasCorrection
                correction_factor = applyCorrectionFactor(n, p, Alpha);
                V_MCD = correction_factor*V_MCD;
            end
            
            % Apply reweighting for efficiency step
            % Compute weighted mean and covariance
            hardLimIdx = (Xs - t_MCD).^2 / V_MCD <= chi2inv(0.975,p);
            
            Xw = Xs(hardLimIdx);
            t_w = sum(Xw) / numel(Xw);
            V_w = sum((Xw - t_w).^2)/(sum(hardLimIdx)-1);
          
            % Introduce a new factor c1 to obtain consistency at the normal
            % model (i.e. c1*Sigma is a consistent estimate of the reweighted
            % robust estimate).
            g = sum(hardLimIdx);
            c1 = computeConsistencyFactor(g, n, p);
            V_MCDw = c1*V_w;
            
            % If we use small-sample corrections, multiply again by a
            % correction factor to make the reweighted estimate unbiased for
            % small samples.
            if useBiasCorrection
                correction_factor = applyCorrectionFactor_rew(n, p, Alpha);
                V_MCDw = correction_factor*V_MCDw;
            end
            
            % Robust location and scale estimates are the reweighting step
            % estimates.
            % Rescale the dataset.
            Mu = t_w*MADd + MEDd;
            Sig = V_MCDw*MADd*MADd;
            
            sqMahal = (X - Mu).^2 / Sig;
            Mahal(~rowswithnans, :) = sqrt(sqMahal);
            Mahal(rowswithnans, :) = NaN;
            Outliers = false(NObs, 1);
            Outliers(~rowswithnans, :) = sqMahal > chi2inv(0.975, p);
            
            if nargout > 4
                % Return results in a struct as the fifth output argument.
                rsltStruct.Mu = Mu;
                rsltStruct.Sigma = Sig;
                rsltStruct.Method = Method;
                rsltStruct.Distances = Mahal;
                rsltStruct.Outliers = Outliers;
                rsltStruct.BiasCorrection = useBiasCorrection;
                rsltStruct.NumTrials = NumTrials;
                rsltStruct.OutlierFraction = valOutlierFraction;
            end
            return;
        end
        
        % Special case 2: outlier fraction is 0. In this case h=n, and the
        % reweighted classical estimates are returned.
        if h==n
            t_MCD = sum(Xs,1)/n;
            V_MCD = cov(Xs);
            
            % Compute reweighted mean and covariance
            [t_w, V_MCDw] = WHardlimMeanCov(Xs, t_MCD, V_MCD);
            
        else
            % In the following lines, the FAST-MCD algorithm covers cases
            % h < n, p > 1.
            
            % The following section implements partitioning and nested
            % extensions, as described in section 3.3 of the Rousseuw-Van
            % Driessen paper. When N > Nmerged, generate a nested system of 5
            % subsets of size 300. Together they form a merged set, which is in
            % turn a proper subset of the full dataset of size n.
            Nmerged = 1500;
            if n <= 600
                % When n <= 600, use the full dataset and do not apply
                % partitioning.
                apply_pooling = false;
                n_sub = n;
                n_groups = 1;
                Part = {1:n};
            else
                apply_pooling = true;
                n_sub = 300;
                if n >= Nmerged
                    % n>=Nmerged
                    % Draw 1500 observations, one by one, without replacement. The first
                    % 300 observations are put in the first subset, and so on, for a total
                    % of 5 subsets.
                    n_groups = 5;
                    mergedSetInds = randsample(n, Nmerged);
                    Part = cell(n_groups,1);
                    
                    % h-subset size to apply in final ieration.
                    h_merged = ceil(h/n*Nmerged);
                    for k=1:n_groups
                        Part{k} = mergedSetInds((k-1)*n_sub + 1:k*n_sub);
                    end
                else %600 < n < Nmerged
                    % When 600 < n < 1500, partition the data into at most 4
                    % subsets of size 300 observations so that each observation
                    % belongs to a subset and so that each subset has about the
                    % same size. For example, 601 is split as 300+301 and 900
                    % as 450+450. For n=901, use 300+300+301; for n=1499, use
                    % 375+375+375+374.
                    n_groups = floor((n-1)/n_sub);
                    siz = floor(n/n_groups);
                    R = rem(n, n_groups);
                    mergedSetInds = randsample(n, n);
                    Part = cell(n_groups,1);
                    
                    % h-subset size to apply in final ieration.
                    h_merged = h;
                    
                    for k=1:n_groups
                        if k<=R
                            Part{k} = mergedSetInds((k-1)*(siz+1) + 1:k*(siz+1));
                        else
                            Part{k} = mergedSetInds((k-1)*siz + 1:k*siz);
                        end
                    end
                end
            end
            
            % hasSingularSubset is a check for singularity of an h-subset
            % of the original dataset. If a singular h-subset is found,
            % then the solution returns a singular covariance matrix (with
            % of course, zero determinant). This singular covariance matrix
            % is returned as the final FMCD solution.
            hasSingularSubset = false;            
            
            % For each subset in the dataset (only one subset if n<=600),
            % repeat NumTrials/n_groups times: construct an initial subset of
            % size h_sub=n_sub*(h/n), carry out 2 C-steps using n_sub and
            % h_sub, and keep the 10 best results.
            Niter = ceil(NumTrials/n_groups);
            h_sub = ceil(h/n*n_sub);
            valsToKeep = cell(n_groups, 1);
                        
            for Ng = 1:n_groups
                DataSubset = Xs(Part{Ng}, :);
                [valsToKeep{Ng}, hasSingularSubset] = runNiterIterationsOnInitialEstimates(DataSubset,...
                    Niter, h_sub);
            end
            
            if ~hasSingularSubset
                % Merge all returned results into one cell array (i.e. pool the
                % results).
                mergedValsToKeep = reshape([valsToKeep{:}], min(Niter, 10)*n_groups, 1);
                if apply_pooling
                    % Do second iterative phase if data partitioning is
                    % used: in the merged set, for each of the 50 returned
                    % solutions, carry out 2 C-steps (using h_merged and
                    % Nmerged), and keep the 10 best solutions.
                    [finalMergedValsToKeep, hasSingularSubset] = runNiterIterationsOnCurrentEstimates(Xs(mergedSetInds,:),...
                        numel(mergedValsToKeep), h_merged, mergedValsToKeep);
                else
                    finalMergedValsToKeep = mergedValsToKeep;
                end
            else
                finalMergedValsToKeep = valsToKeep;
            end
            
            if ~hasSingularSubset
                % Final step: Apply C-steps to convergence on the full dataset,
                % using the best results from the previous phase.
                [finalSolution, hasSingularSubset] = runNiterIterationsOnFinalEstimates(Xs, h, finalMergedValsToKeep);
            else
                warning(message('stats:robustcov:SingularHSubsetFound'));
                finalSolution = finalMergedValsToKeep{:};
            end
            
            % Return the full raw MCD estimates.
            t_full = finalSolution{1};
            V_full = finalSolution{2};
            
            % Multiply the covariance estimate by a factor c0 to obtain
            % consistency at the normal model (i.e. c0*Sigma is a consistent
            % estimate of the covariance matrix).
            t_MCD = t_full;
            c0 = computeConsistencyFactor(h, n, p);
            V_MCD = c0*V_full;
            
            % If we use small-sample corrections, multiply by a correction
            % factor to make the estimate unbiased for small samples.
            if useBiasCorrection
                correction_factor = applyCorrectionFactor(n, p, Alpha);
                V_MCD = correction_factor*V_MCD;
            end
            
            if ~hasSingularSubset
                % Only apply reweigting if no singular subset solution has been found.
                
                % Now apply a reweighted version of MCD for improved estimate. This
                % improves the efficiency without affecting the robustness of the
                % estimate. Also return the number of observations retained as
                % inliers.
                [t_w, V_MCDw, g] = WHardlimMeanCov(Xs, t_MCD, V_MCD);
                
                % Introduce a new factor c1 to obtain consistency at the normal
                % model (i.e. c1*Sigma is a consistent estimate of the reweighted
                % robust estimate).
                c1 = computeConsistencyFactor(g, n, p);
                V_MCDw = c1*V_MCDw;
            else
                t_w = t_MCD;
                V_MCDw = V_MCD;
            end
            % If we use small-sample corrections, multiply again by a
            % correction factor to make the reweighted estimate unbiased for
            % small samples.
            if useBiasCorrection
                correction_factor = applyCorrectionFactor_rew(n, p, Alpha);
                V_MCDw = correction_factor*V_MCDw;
            end
        end
        
        % Rescale the dataset.
        Mu = t_w.*MADd + MEDd;
        Sig = V_MCDw.*MADd(ones(p,1), :).*MADd(ones(p,1), :)';
        
        if nargout > 4
            rsltStruct.BiasCorrection = useBiasCorrection;
            rsltStruct.NumTrials = NumTrials;
            rsltStruct.OutlierFraction = valOutlierFraction;
        end
end

% Return final robust distances and indices of observations returned as
% outliers.
if nargout > 2
    Mahal = zeros(NObs,1);  % Observations with missing values are set to NaN.
    if ~ZeroCovDet(Sig) % Check for singularity of the covariance.
        Mahal(~rowswithnans, :) = sqrt(localSquaredMahal(X, Mu, Sig));
        Mahal(rowswithnans, :) = NaN;
    end    
    if nargout > 3
        Outliers = false(NObs, 1);
        Outliers(~rowswithnans, :) = Mahal(~rowswithnans, :) > sqrt(chi2inv(0.975, p));
    end
    if nargout > 4
        % Return results in a struct as the fifth output argument.
        rsltStruct.Mu = Mu;
        rsltStruct.Sigma = Sig;
        rsltStruct.Method = Method;
        rsltStruct.Distances = Mahal;
        rsltStruct.Outliers = Outliers;
    end
end

%----------------Subfunctions--------------------------------------------
%------------------------------------------
function [t, V] = computeStartEstimate(X, currentStart)

if ischar(currentStart)
    switch currentStart
        case 'medianball'
            % The Median Ball (MB) estimator uses the classical estimator
            % computed from the cases with Di(med(X), eye(p)) <=
            % med(Di(med(X), eye(p)))
            medianBallDists = sum(bsxfun(@minus, X, median(X)).^2, 2);
            medianBallObsInds = medianBallDists <= median(medianBallDists); % pdist2(X, med_X, 'mahal', I_p);
            [V, t] = localMeanCov(X(medianBallObsInds, :));
            
        case 'classical'
            % The classical estimator (or DGK estimator) uses the classical
            % estimator as the start.
            [V, t] = localMeanCov(X);
            
        case 'elemental'
            % The elemental estimator draws (p+1) samples at random from
            % the dataset. The estimate is evaluated over this subset.
            [t, V] = elementalSubsetEstimate(X);
            
        otherwise
            % User-defined start
            startFcn = str2func(currentStart);
            [t, V] = startFcn(X);
    end
else %Start is a user-defined function handle.
    startFcn = currentStart;
    [t, V] = startFcn(X);
end

%------------------------------------------
function [Tinit, Sinit, h] = elementalSubsetEstimate(X)
% Draw a random initial (p+1)-subset J of X. Compute T0=mean(J), S0=cov(J).
% If S0 is singular (det(S0)=0), add another random observation until it is
% non-singular. Return initial mean and covariance estimates.

[n, p] = size(X);
h = randsample(n, p+1);

% Initial subset
J = X(h, :);

% Test for positive definiteness of covariance matrix using CHOLCOV
C = cov(J);
hasZeroCovDet = ZeroCovDet(C);

numNewPoints = 0;
while hasZeroCovDet && numNewPoints < n-p  % Add observations until det(cov(J))>0
    hnew = randsample(n, 1);
    h = unique([h; hnew]);
    J = X(h, :);
    C = cov(J);
    hasZeroCovDet = ZeroCovDet(C);
    numNewPoints = numNewPoints + 1;
end

Tinit = sum(J,1) / size(J, 1);
Sinit = C;

%------------------------------------------
function [hasZeroCovDet, C] = ZeroCovDet(S)
% Test for positive definiteness of covariance matrix using CHOL. Return
% Cholesky factor when it is needed for computing Squared Mahal distances.
[C, p] = chol(S);
hasZeroCovDet = p > 0;

%------------------------------------------
function [T, S, sortedInds, hasSingularSubset] = computeCStep(X, T, S, h)
% Perform a C-step, or "concentration" step. We "concentrate" on the h
% observations with smallest relative distances given the pair (Told,
% Sold). Return Tnew and Snew, the mean and covariance matrix based on
% these h smallest observations.

[hasSingularSubset, C] = ZeroCovDet(S);

if hasSingularSubset
    % If the current estimate is not positive definite, then a solution
    % with minimum covariance determiannt is found since det(cov(Xnew=0)).
    % We thus return the current estimate (T, S) as the new estimate.
    sortedInds = 1:h;
else
    % 1. Compute the distances dold(i)=dist(xi), i=1,..,n.
    Dold = sum((bsxfun(@minus, X, T) / C).^2, 2); %localSquaredMahal(X, T, S);
    
    % 2. Sort these distances, and return sorted indices.
    [~, sorted] = sort(Dold);
    
    % 3. Retain a new set of the h smallest observations.
    sortedInds = sorted(1:h);
    
    Xnew = X(sortedInds, :);
    
    % 4. Return new mean and covariance.
    [S, T] = localMeanCov(Xnew);
end

%------------------------------------------
function [S, T] = localMeanCov(X)
% Local function for computing classical mean and covariance. 
n = size(X, 1);
T = sum(X, 1)/n;
Xc = bsxfun(@minus, X, T);
S = (Xc'*Xc)/(n-1);

%------------------------------------------
function d_mahal = localSquaredMahal(X, T, S)
% Local computation of Mahalanobis distances
C = chol(S);
d_mahal = sum((bsxfun(@minus, X, T) / C).^2, 2);

%------------------------------------------
function c = computeConsistencyFactor(h, n, p)
% Subfunction for computing consistency factors, both for the raw and
% reweighted estimates.
if h==n
    c = 1;
else
    quan = h/n;
    c = quan/chi2cdf(chi2inv(quan, p), p+2);
end

%------------------------------------------
function [tw, Cw, nI] = WHardlimMeanCov(X, t, C)
% Subfunction for computing reweighted mean and covariance estimates based
% on the hard limit weight thresold: w = I(d_mahal(x,T,C) <= sqrt(chi2(p, 0.975)))
% Also return the number of inliers as this will be needed.

% 1. Compute Squared Mahalanobis distances d(Xi)=mahal(Xi, t, V)^ 2 for each
% row of X, where t and V are the raw MCD estimates of the mean and
% covariance.
sqMalahDists = localSquaredMahal(X, t, C);

% 2. Compute hard rejection threshold used in the hard rejection
% weight function W(d) = I(d<=d0) used for discarding observations
% whose Mahalanobis distances are above this value.
rej = 0.975; % Rejection threshold.
p = size(X,2);
chi2pv = chi2inv(rej,p); %0.975

% 3. Compute weighted mean and covariance
hardLimIdx = sqMalahDists<=chi2pv;
nI = sum(hardLimIdx);
tw =  sum(X(hardLimIdx, :))./ nI;
Xc = bsxfun(@minus, X(hardLimIdx,:), tw);
Cw = (Xc'*Xc)/(nI - 1);

%------------------------------------------
function [sig, mu] = tauscale(X) %#ok<DEFNU>
% Compute robust univariate location and scale estimates using the "tau
% scale" of Yohai and Zamar (1998). This is a truncated standard deviation
% of x, where x is a univariate sample to be estimated. Also return a
% weighted mean mu. Input is a data matrix X. Compute tau-scale estimates
% along the columns of X.
[n, p] = size(X);

% If X is a vector make it a column vector in all cases.
if n==1
    X = X';
    n = p;
end

c1 = 4.5;
W_c1 = @(x)((1 - (x/c1).^2).^2).*(abs(x) <= c1);
MAD = mad(X, 1);
MED = median(X, 1);
W = W_c1(bsxfun(@rdivide, bsxfun(@minus, X, MED), MAD));

% For univariate estimates with zero MAD, set weights to 1.
W(:, MAD == 0) = 1;

% Location estimate
mu = sum(X.*W)./sum(W);

% Univariate dispersion estimate
c2 = 3;
rho_c = @(x)min(x.^2, c2^2);
sig = sqrt(sum(rho_c(bsxfun(@rdivide, bsxfun(@minus, X, mu), MAD))).*(MAD.^2/n));

%------------------------------------------
function [sig, mu] = qn(X) %#ok<DEFNU>
% Compute robust univariate location and scale estimates using the "Qn"
% estimate of Rousseuw and Croux (1993). The Qn estimator is defined as:
%   Qn = d*{|xi-xj|; i<j}(k)
% i.e. the k-th quantile of all pairwise distances between two data points,
% where k=nchoosek(h,2) and h=floor(n/2)+1 is roughly half the number of
% observations. This estimator is considered as an efficient alternative to
% MAD, with 50% breakdown and high efficiency at the normal distribution.
% If the robust location is needed, the median is returned.
% Input is a data matrix X. Compute Qn estimates along the columns of X.
% Reference:
% P.J. Rousseeuw, C. Croux (1993) Alternatives to the Median Absolute
% Deviation, JASA.
[n, p] = size(X);

% If X is a vector make it a column vector in all cases.
if n==1
    X = X';
    n = p;
    p = 1;
end

h = floor(n/2)+1;
Qn = zeros(1,p);

if n > 1
    for j=1:p
        Xj = X(:,j);
        D = abs(bsxfun(@minus, Xj, Xj'));
        sortedDistVals = sort(D(triu(ones(n),1) > 0));
        % Return Qn estimate, multiplied by a consistency factor to obtsin
        % consistency at the normal distribution.
        Qn(j) = 2.21914*sortedDistVals(h*(h-1)/2);
    end
end

% Also multiply by c_n, a small sample correction factor.
if n <= 12
    c_n = [0, 0.399356, 0.99365, 0.51321, 0.84401, 0.6122,...
        0.85877, 0.66993, 0.87344, 0.72014, 0.88906, 0.75743];
    Qn = c_n(n)*Qn;
    
else
    if mod(n,2) == 1 % odd n
        Qn = n/(n+1.4)*Qn;
    else % even n
        Qn = n/(n+3.8)*Qn;
    end
end
sig = Qn;
mu = median(X);

%------------------------------------------
function [valsToKeep, hasSingularSubset] = runNiterIterationsOnInitialEstimates(DataSubset,...
    Niter, h)
% Subfunction for computing a small number of initial C-steps (here 2) on a
% subsample of the data. Return the 10 best results in a cell array, where
% each element contains the mean, covariance and objective (determinant)
% of that sample.
numLowestDet = min(Niter, 10); % This handles the edge case where nIter < 10.
valsInNiter = cell(Niter, 1);

hasSingularSubset = false;
for N=1:Niter
    
    % Start from an initial (p+1) random subset
    [T0, S0] = elementalSubsetEstimate(DataSubset);
    
    % Carry out a small number of C-steps.
    numCsteps = 2;
    for c=1:numCsteps
        if ~hasSingularSubset
            [T0, S0, ~, hasSingularSubset] = computeCStep(DataSubset, T0, S0, h);
        end
    end
    
    % Save current estimates for later comparison.
    detS0 = det(S0);
    valsInNiter{N} = {T0, S0, detS0};
    if hasSingularSubset
        % A singular subset was found. Exit program and return this
        % estimate as the final solution.
        valsToKeep = valsInNiter{N};
        return;
    end
end

% Keep the estimates with lowest objective function values.
[~, sortedDetInds] = sort(cellfun(@(x)x{3}, valsInNiter));
valsToKeep = valsInNiter(sortedDetInds(1:numLowestDet));

%------------------------------------------
function [valsToKeep, hasSingularSubset] = runNiterIterationsOnCurrentEstimates(DataSubset,...
    Niter, h, currentSolutions)
% For a subset of the data set (a disjoint subset of size n_s in the case
% of nested partitioning with n >= 600, the full data set otherwise for
% small n), from an initial estimate of mean and covariance, perform C=2
% concentration steps based on n_s and h_s = n_s*(h/n). Repeat Niter times
% and retain the best = 10 estimates with lowest covariance determinant.
numLowestDet = min(Niter,10); % This handles the edge case where nIter < 10.

hasSingularSubset = false;
for N=1:Niter
    T0 = currentSolutions{N}{1};
    S0 = currentSolutions{N}{2};
    
    % Carry out C-steps
    numCsteps = 2;
    for c=1:numCsteps
        if ~hasSingularSubset
            [T0, S0, ~, hasSingularSubset] = computeCStep(DataSubset, T0, S0, h);
        end
    end
    detS0 = det(S0);
    
    % Replace the new value for the covariance determinant in the current
    % solutions
    currentSolutions{N}{3} = detS0;
    if hasSingularSubset
        % A singular subset was found. Exit program and return this
        % estimate as the final solution.
        valsToKeep = currentSolutions{N};
        return;
    end
end

[~, lowestDetsInds] = sort(cellfun(@(x)x{3}, currentSolutions));
valsToKeep = currentSolutions(lowestDetsInds(1:numLowestDet));

%------------------------------------------
function [valsToKeep, hasSingularSubset] = runNiterIterationsOnFinalEstimates(DataSubset, h,...
    currentSolutions)
% For the final results with lowest values of covariance determinant, carry
% out C-steps until convergence on the final data set.

numCstepsf = 100;
bestObj = Inf;
Niter = length(currentSolutions);

hasSingularSubset = false;

for N=1:Niter
    T0 = currentSolutions{N}{1};
    S0 = currentSolutions{N}{2};
    detS0 = currentSolutions{N}{3};
    
    [T0, S0] = computeCStep(DataSubset, T0, S0, h);
    detS = det(S0);
    c = 1;
    while detS ~= detS0 && c <= numCstepsf
        detS0 = detS;
        if ~hasSingularSubset
            [T0, S0, ~, hasSingularSubset] = computeCStep(DataSubset, T0, S0, h);
        end
        detS = det(S0);
        c = c+1;
    end
    
    if detS < bestObj
        bestObj = detS;
        T = T0;
        S = S0;
    end
    if hasSingularSubset
        % A singular subset was found. Exit program and return this
        % estimate as the final solution.
        valsToKeep = currentSolutions{N};
        return;
    end
end

valsToKeep = {T, S, bestObj};


% %------------------------------------------
% Following are subfunctions for calculating small-sample correction
% factors. They are calculated using a function that was fit to
% approximate the correction factors obtained from simulations. For fixed p
% and a given value of alpha, the model used as a function of n is:
%  f_alpha_p(n) = 1 + Gamma / n^Beta
% For each p and alpha are corresponding parameters Gamma = Gamma(alpha, p)
% and Beta = Beta(alpha, p).
% For p = 1,2, the correction factor is c_alpha_p(n) = 1 / f_alpha_p(n).
% For p > 2, we first solve the system of equations:
%
%  1 + N_alpha_3 / P^K_alpha_3 = 1 + Gamma_alpha_p / (2*p^2)^Beta_alpha_p
%  1 + N_alpha_5 / P^K_alpha_5 = 1 + Gamma_alpha_p / (3*p^2)^Beta_alpha_p
%
% to obtain estimates for Gamma_alpha_p, Beta_alpha_p (after rewriting the
% above equation into a linear system by taking logarithms). The
% corresponding correction factor is c_alpha_p(n) = 1 / f_alpha_p(n), where
% f_alpha_p(n) = 1 + Gamma_alpha_p / n^Beta_alpha_p
% Empirically, it was found that for a given n, p, the estimate is
% approximately linear as a function of alpha when 0.5 <= alpha <= 0.875
% and for 0.875 <= alpha <= 1 (the correction factor is 1 for alpha = 1).
%
% Similarly, the reweighted MCD estimates are not unbiased at small samples
% even when the consistency factor is included. A similar procedure is
% therefore carried out on the reweighted estimates.

%------------------------------------------
function corr_factor_alpha = applyCorrectionFactor(n, p, Alpha)
% From the simulations, it was found emprically that for fixed n and p the
% mean function is approximately a linear function of alpha. The values
% beween 0.5 and 0.875 are determined by linear interpolation. There is no
% correction factor for alpha=1 and values between 0.875 and 1 are also
% determined using linear interpolation.

corr_factor500 = applyCorrectionFactor_alpha500(n,p);
corr_factor875 = applyCorrectionFactor_alpha875(n,p);
if (Alpha >= 0.5) && (Alpha <=0.875)
    OneOverCorr_factor_alpha = interp1([0.5, 0.875], [corr_factor500, corr_factor875], Alpha);
else
    OneOverCorr_factor_alpha = interp1([0.875, 1], [corr_factor875, 1], Alpha);
end

corr_factor_alpha = 1 / OneOverCorr_factor_alpha;

%------------------------------------------
function corr_factor_alpha = applyCorrectionFactor_rew(n, p, Alpha)
% Similar procedure to 'applyCorrectionFactor' subfunction for the
% reweighted estimates.

corr_factor500 = applyCorrectionFactor_alpha500_rew(n,p);
corr_factor875 = applyCorrectionFactor_alpha875_rew(n,p);
if (Alpha >= 0.5) && (Alpha <=0.875)
    OneOverCorr_factor_alpha = interp1([0.5, 0.875], [corr_factor500, corr_factor875], Alpha);
else
    OneOverCorr_factor_alpha = interp1([0.875, 1], [corr_factor875, 1], Alpha);
end

corr_factor_alpha = 1 / OneOverCorr_factor_alpha;

%------------------------------------------
function correction_factor_alpha500_p_n = applyCorrectionFactor_alpha500(n, p)
% Small sample correction factors for alpha = 0.5.
% Note: TO DO: Run more simulations for more precise results.

if p==1
    Coeffs_alpha_p_500 = [0.262024211897096; 0.604756680630497];
    fp_500_n = 1 - exp(Coeffs_alpha_p_500(1))/ (n^Coeffs_alpha_p_500(2));
elseif p==2
    Coeffs_alpha_p_500 = [0.723382493037778; 0.698065251669507];
    fp_500_n = 1 - exp(Coeffs_alpha_p_500(1))/ (n^Coeffs_alpha_p_500(2));
else
    B_alpha500_p = [0.356026735711780 - 1.262633369321510*log(p);...
        0.059599305599129 - 1.289079914403870*log(p)];
    A_alpha500_p = [1  -log(2*p^2); 1  -log(3*p^2)];
    
    Coeffs_alpha_p_500_q23 = A_alpha500_p \ B_alpha500_p;
    Gam_alpha500_p = -exp(Coeffs_alpha_p_500_q23(1));
    Bet_alpha500_p = Coeffs_alpha_p_500_q23(2);
    
    fp_500_n = 1 + Gam_alpha500_p/(n^Bet_alpha500_p);
    
end
correction_factor_alpha500_p_n = fp_500_n;


%------------------------------------------
function correction_factor_alpha875_p_n = applyCorrectionFactor_alpha875(n, p)
% Small sample correction factors for alpha = 0.875.
% Note: TO DO: Run more simulations for more precise results.

if p==1
    Coeffs_alpha_p_500 = [-0.351584646688712; 1.01646567502486];
    fp_875_n = 1 - exp(Coeffs_alpha_p_500(1))/ (n^Coeffs_alpha_p_500(2));
elseif p==2
    Coeffs_alpha_p_875 = [0.446537815635445; 1.06690782995919];
    fp_875_n = 1 - exp(Coeffs_alpha_p_875(1))/ (n^Coeffs_alpha_p_875(2));
else
    B_alpha875_p = [-0.787063511268168 - 1.111925412787940*log(p);...
        -1.223355411552061 - 1.096493291498110*log(p)];
    
    A_alpha875_p = [1  -log(2*p^2); 1  -log(3*p^2)];
    
    Coeffs_alpha_p_875_q23 = A_alpha875_p \ B_alpha875_p;
    Gam_alpha875_p = -exp(Coeffs_alpha_p_875_q23(1));
    Bet_alpha875_p = Coeffs_alpha_p_875_q23(2);
    
    fp_875_n = 1 + Gam_alpha875_p/(n^Bet_alpha875_p);
    
end
correction_factor_alpha875_p_n = fp_875_n;

%------------------------------------------
function correction_factor_alpha500_p_n_rew = applyCorrectionFactor_alpha500_rew(n, p)
% Small sample correction factors for alpha = 0.5 - reweighted estimates

% Note: TO DO: Run more simulations for more precise results.

if p==1
    Coeffs_alpha_p_500 = [1.11098143415027; 1.5182890270453];
    fp_500_n = 1 - exp(Coeffs_alpha_p_500(1))/ (n^Coeffs_alpha_p_500(2));
elseif p==2
    Coeffs_alpha_p_500 = [1.963349082029735; 1.472294837772169];
    fp_500_n = 1 - exp(Coeffs_alpha_p_500(1))/ (n^Coeffs_alpha_p_500(2));
else
    B_alpha500_p = [0.028029212868299 - 1.67659883081926*log(p);...
        -1.316758095133297 - 1.35968562893582*log(p)];
    
    A_alpha500_p = [1  -log(2*p^2); 1  -log(3*p^2)];
    
    Coeffs_alpha_p_500_q23 = A_alpha500_p \ B_alpha500_p;
    Gam_alpha500_p = -exp(Coeffs_alpha_p_500_q23(1));
    Bet_alpha500_p = Coeffs_alpha_p_500_q23(2);
    
    fp_500_n = 1 + Gam_alpha500_p/(n^Bet_alpha500_p);
end

correction_factor_alpha500_p_n_rew = fp_500_n;


%------------------------------------------
function correction_factor_alpha875_p_n_rew = applyCorrectionFactor_alpha875_rew(n, p)
% Small sample correction factors for alpha = 0.875 - reweighted estimates

if p==1
    Coeffs_alpha_p_875 = [-0.66046776772861; 0.88939595831888];
    fp_875_n = 1 - exp(Coeffs_alpha_p_875(1))/ (n^Coeffs_alpha_p_875(2));
elseif p==2
    Coeffs_alpha_p_875 = [0.79473550581058; 1.10081930350091];
    fp_875_n = 1 - exp(Coeffs_alpha_p_875(1))/ (n^Coeffs_alpha_p_875(2));
else
    B_alpha875_p = [-0.607919580335715 - 1.259944832222920*log(p);...
        -1.067721154423836 - 1.251590042571330*log(p)];
    
    A_alpha875_p = [1  -log(2*p^2); 1  -log(3*p^2)];
    
    Coeffs_alpha_p_875_q23 = A_alpha875_p \ B_alpha875_p;
    Gam_alpha875_p = -exp(Coeffs_alpha_p_875_q23(1));
    Bet_alpha875_p = Coeffs_alpha_p_875_q23(2);
    
    fp_875_n = 1 + Gam_alpha875_p/(n^Bet_alpha875_p);
end

correction_factor_alpha875_p_n_rew = fp_875_n;
