function [B,stats] = lasso(X,Y,varargin)
%LASSO Perform lasso or elastic net regularization for linear regression.
%   [B,STATS] = lasso(X,Y,...) Performs L1-constrained linear least  
%   squares fits (lasso) or L1- and L2-constrained fits (elastic net)
%   relating the predictors in X to the responses in Y. The default is a
%   lasso fit, or constraint on the L1-norm of the coefficients B.
%
%   Positional parameters:
%
%     X                A numeric matrix (dimension, say, NxP)
%     Y                A numeric vector of length N
%   
%   Optional input parameters:  
%
%     'Weights'        Observation weights.  Must be a vector of non-negative
%                      values, of the same length as columns of X.  At least
%                      two values must be positive. (default ones(N,1) or 
%                      equivalently (1/N)*ones(N,1)).
%     'Alpha'          Elastic net mixing value, or the relative balance
%                      between L2 and L1 penalty (default 1, range (0,1]).
%                      Alpha=1 ==> lasso, otherwise elastic net.
%                      Alpha near zero ==> nearly ridge regression.
%     'NumLambda'      The number of lambda values to use, if the parameter
%                      'Lambda' is not supplied (default 100).  Ignored
%                      if 'Lambda' is supplied.  LASSO may return fewer
%                      fits than specified by 'NumLambda' if the residual
%                      error of the fits drops below a threshold percentage 
%                      of the variance of Y.
%     'LambdaRatio'    Ratio between the minimum value and maximum value of
%                      lambda to generate, if the  parameter "Lambda" is not 
%                      supplied.  Legal range is [0,1). Default is 0.0001.
%                      If 'LambdaRatio' is zero, LASSO will generate its
%                      default sequence of lambda values but replace the
%                      smallest value in this sequence with the value zero.
%                      'LambdaRatio' is ignored if 'Lambda' is supplied.
%     'Lambda'         Lambda values. Will be returned in return argument
%                      STATS in ascending order. The default is to have LASSO
%                      generate a sequence of lambda values, based on 'NumLambda'
%                      and 'LambdaRatio'. LASSO will generate a sequence, based
%                      on the values in X and Y, such that the largest LAMBDA                 
%                      value is just sufficient to produce all zero coefficients B.
%                      You may supply a vector of real, non-negative values of 
%                      lambda for LASSO to use, in place of its default sequence.
%                      If you supply a value for 'Lambda', 'NumLambda' and 
%                      'LambdaRatio' are ignored.
%     'DFmax'          Maximum number of non-zero coefficients in the model.
%                      Can be useful with large numbers of predictors.
%                      Results only for lambda values that satisfy this
%                      degree of sparseness will be returned. Default is
%                      to not limit the number of non-zero coefficients.
%     'Intercept'      Flag for fitting the model with an intercept term.
%                      The default value is true, which indicates to include
%                      the intercept term in the model. If the 'Intercept'
%                      value is false, the returned intercept value is 0.
%     'Standardize'    Flag for centering and scaling X prior to fitting the 
%                      model sequence. This setting determines whether the 
%                      regularization is applied to the coefficients on the 
%                      standardized scale or the original scale. The results 
%                      are always presented on the original data scale. The 
%                      default value is true, which indicates to center and 
%                      scale X. When the 'Intercept' value is false, the 
%                      'Standardize' value is set to false.
%                      Note: X and Y are always centered when the 'Intercept' 
%                      value is true.
%     'RelTol'         Convergence threshold for coordinate descent algorithm.
%                      The coordinate descent iterations will terminate
%                      when the relative change in the size of the
%                      estimated coefficients B drops below this threshold.
%                      Default: 1e-4. Legal range is (0,1).
%     'CV'             If present, indicates the method used to compute MSE.
%                      When 'CV' is a positive integer K, LASSO uses K-fold
%                      cross-validation.  Set 'CV' to a cross-validation 
%                      partition, created using CVPARTITION, to use other
%                      forms of cross-validation. You cannot use a
%                      'Leaveout' partition with LASSO.                
%                      When 'CV' is 'resubstitution', LASSO uses X and Y 
%                      both to fit the model and to estimate the mean 
%                      squared errors, without cross-validation.  
%                      The default is 'resubstitution'.
%     'MCReps'         A positive integer indicating the number of Monte-Carlo
%                      repetitions for cross-validation.  The default value is 1.
%                      If 'CV' is 'resubstitution' or a cvpartition of type
%                      'resubstitution', 'MCReps' must be 1.  If 'CV' is a
%                      cvpartition of type 'holdout', then 'MCReps' must be
%                      greater than one.
%     'MaxIter'        Maximum number of iterations allowed.  Default is 1e5.
%     'PredictorNames' A string/cell array of names for the predictor variables,
%                      in the order in which they appear in X. 
%                      Default: {}
%     'Options'        A structure that contains options specifying whether to
%                      conduct cross-validation evaluations in parallel, and
%                      options specifying how to use random numbers when computing
%                      cross validation partitions. This argument can be created
%                      by a call to STATSET. CROSSVAL uses the following fields:
%                        'UseParallel'
%                        'UseSubstreams'
%                        'Streams'
%                      For information on these fields see PARALLELSTATS.
%                      NOTE: If supplied, 'Streams' must be of length one.
%   
%   Return values:
%     B                The fitted coefficients for each model. 
%                      B will have dimension PxL, where 
%                      P = size(X,2) is the number of predictors, and
%                      L = length(lambda).
%     STATS            STATS is a struct that contains information about the
%                      sequence of model fits corresponding to the columns
%                      of B. STATS contains the following fields:
%
%       'Intercept'    The intercept term for each model. Dimension 1xL.
%       'Lambda'       The sequence of lambda penalties used, in ascending order. 
%                      Dimension 1xL.
%       'Alpha'        The elastic net mixing value that was used.
%       'DF'           The number of nonzero coefficients in B for each
%                      value of lambda. Dimension 1xL.
%       'MSE'          The mean squared error of the fitted model for each
%                      value of lambda. If cross-validation was performed,
%                      the values for 'MSE' represent Mean Prediction
%                      Squared Error for each value of lambda, as calculated 
%                      by cross-validation. Otherwise, 'MSE' is the mean
%                      sum of squared residuals obtained from the model
%                      with B and STATS.Intercept.
%
%     If cross-validation was performed, STATS also includes the following
%     fields:
%
%       'SE'           The standard error of MSE for each lambda, as
%                      calculated during cross-validation. Dimension 1xL.
%       'LambdaMinMSE' The lambda value with minimum MSE. Scalar.
%       'Lambda1SE'    The largest lambda such that MSE is within 
%                      one standard error of the minimum. Scalar.
%       'IndexMinMSE'  The index of Lambda with value LambdaMinMSE.
%       'Index1SE'     The index of Lambda with value Lambda1SE.
%
%     Examples:
%
%        % (1) Run the lasso on data obtained from the 1985 Auto Imports Database 
%        % of the UCI repository.  
%        % http://archive.ics.uci.edu/ml/machine-learning-databases/autos/imports-85.names
%        load imports-85;
%        Description
%
%        % Extract Price as the response variable and extract non-categorical
%        % variables related to auto construction and performance
%        %
%        X = X(~any(isnan(X(:,1:16)),2),:);
%        Y = X(:,16);
%        Y = log(Y);
%        X = X(:,3:15);
%        predictorNames = {'wheel-base' 'length' 'width' 'height' ...
%            'curb-weight' 'engine-size' 'bore' 'stroke' 'compression-ratio' ...
%            'horsepower' 'peak-rpm' 'city-mpg' 'highway-mpg'};
%
%        % Compute the default sequence of lasso fits.
%        [B,S] = lasso(X,Y,'CV',10,'PredictorNames',predictorNames);
%
%        % Display a trace plot of the lasso fits.
%        axTrace = lassoPlot(B,S);
%        % Display the sequence of cross-validated predictive MSEs.
%        axCV = lassoPlot(B,S,'PlotType','CV');
%        % Look at the kind of fit information returned by lasso.
%        S
%
%        % What variables are in the model corresponding to minimum 
%        % cross-validated MSE, and in the sparsest model within one 
%        % standard error of that minimum.
%        minMSEModel = S.PredictorNames(B(:,S.IndexMinMSE)~=0)
%        sparseModel = S.PredictorNames(B(:,S.Index1SE)~=0)
%
%        % Fit the sparse model and examine residuals.
%        fitSparse = S.Intercept(S.Index1SE) + X*B(:,S.Index1SE);
%        corr(fitSparse,Y-fitSparse)
%        figure
%        plot(fitSparse,Y-fitSparse,'o')
%
%        % Consider a slightly richer model. A model with 6 variables may be a 
%        % reasonable alternative.  Find the index for a corresponding fit.
%        df6index = min(find(S.DF==6));
%        fitDF6 = S.Intercept(df6index) + X*B(:,df6index);
%        corr(fitDF6,Y-fitDF6)
%        plot(fitDF6,Y-fitDF6,'o')         
%         
%        % (2) Run lasso on some random data with 250 predictors
%        %
%        n = 1000; p = 250;
%        X = randn(n,p);
%        beta = randn(p,1); beta0 = randn;
%        Y = beta0 + X*beta + randn(n,1);
%        lambda = 0:.01:.5;
%        [B,S] = lasso(X,Y,'Lambda',lambda);
%        lassoPlot(B,S);
%
%        % compare against OLS
%        %
%        figure
%        bls = [ones(size(X,1),1) X] \ Y;
%        plot(bls,[S.Intercept; B],'.');
%
%        % Run the same lasso fit but restricting the number of
%        % non-zero coefficients in the fitted model.
%        %
%        [B2,S2] = lasso(X,Y,'Lambda',lambda,'DFmax',12);
%
%   See also lassoPlot, ridge, parallelstats.

%   References: 
%   [1] Tibshirani, R. (1996) Regression shrinkage and selection
%       via the lasso. Journal of the Royal Statistical Society,
%       Series B, Vol 58, No. 1, pp. 267-288.
%   [2] Zou, H. and T. Hastie. (2005) Regularization and variable
%       selection via the elastic net. Journal of the Royal Statistical
%       Society, Series B, Vol. 67, No. 2, pp. 301-320.
%   [3] Friedman, J., R. Tibshirani, and T. Hastie. (2010) Regularization
%       paths for generalized linear models via coordinate descent.
%       Journal of Statistical Software, Vol 33, No. 1,
%       http://www.jstatsoft.org/v33/i01.
%   [4] Hastie, T., R. Tibshirani, and J. Friedman. (2008) The Elements
%       of Statistical Learning, 2nd edition, Springer, New York.

%   Copyright 2011-2020 The MathWorks, Inc.

% -------------------------------------- 
% Sanity check the positional parameters
% --------------------------------------

% X a real 2D matrix
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if ~ismatrix(X) || ~isreal(X) 
    error(message('stats:lasso:XnotaReal2DMatrix'));
end

if size(X,1) < 2
    error(message('stats:lasso:TooFewObservations'));
end

% Y a vector, same length as the columns of X
if ~isvector(Y) || ~isreal(Y) || size(X,1) ~= length(Y)
    error(message('stats:lasso:YnotaConformingVector'));
end

% If Y is a row vector, convert to a column vector
if isrow(Y)
    Y = Y';
end

% This screen (okrows) selects all the predictions and response we can use, 
% but there may be further reductions for zero observation weights.
% Defer further checking of X,Y until weights are pre-screened.
okrows = all(isfinite(X),2) & all(isfinite(Y),2);

% --------------------------------------------------------------------
% Parse and process the optional parameters
% --------------------------------------------------------------------

LRdefault = 0.0001;

pnames = { 'weights' 'alpha' 'numlambda' 'lambdaratio' 'lambda' ...
    'dfmax' 'standardize' 'reltol' 'cv' 'mcreps' 'maxiter'...
    'predictornames' 'options' 'intercept'};
dflts  = { []        1       100       LRdefault     []      ...
     []      true          1e-4    'resubstitution'  1  1e5 ...
     {}               []       true};
[weights, alpha, nLambda, lambdaRatio, lambda, ...
    dfmax, standardize, reltol, cvp, mcreps, maxIter,predictorNames, ParOptions, withConstant] ...
     = internal.stats.parseArgs(pnames, dflts, varargin{:});
validateattributes(withConstant,{'logical'},{'scalar'},'lasso','Intercept');

% === 'Alpha' parameter ===

% Require 0 < alpha <= 1.
% "0" would correspond to ridge, "1" is lasso.
if ~isscalar(alpha) || ~isreal(alpha) || alpha <= 0 || alpha > 1
    error(message('stats:lasso:InvalidAlpha'))
end

% === 'Weights' parameter ===

observationWeights = ~isempty(weights);
if observationWeights
    % This screen works on weights prior to stripping NaNs and Infs.
    if ~isvector(weights) || ~isreal(weights) || size(X,1) ~= length(weights) || ...
            ~all(isfinite(weights)) || any(weights<0) || sum(weights>0) < 2
        error(message('stats:lasso:InvalidObservationWeights'));
    end
    
    okrows = okrows & weights(:)>0; % remove nan here as well
    weights = weights(okrows);
    
end

% Add type checks
internal.stats.checkSupportedNumeric('X',X,false,false,false,false);
internal.stats.checkSupportedNumeric('Y',Y,false,false,false,false);
internal.stats.checkSupportedNumeric('Weights',weights,false,false,false,false);
internal.stats.checkSupportedNumeric('Alpha',alpha,false,false,false,false);
internal.stats.checkSupportedNumeric('Lambda',lambda,false,false,false,false);
internal.stats.checkSupportedNumeric('LambdaRatio',lambdaRatio,false,false,false,false);

% Cast everything to type of X
if isa(X,'single')
    Y = single(Y);
    weights = single(weights);
    alpha = single(alpha);
    lambda = single(lambda);
    lambdaRatio = single(lambdaRatio);
else
    Y = double(Y);
    weights = double(weights);
    alpha = double(alpha);
    lambda = double(lambda);
    lambdaRatio = double(lambdaRatio);
end

% Remove observations with NaNs and Infs in the predictor or response
% or with zero observation weight.
X = X(okrows,:);
Y = Y(okrows);

% We need at least two observations after stripping NaNs and Infs.
if size(X,1) < 2
    error(message('stats:lasso:TooFewObservationsAfterNaNs'));
end

% If X has any constant columns, we want to exclude them from the
% coordinate descent calculations.  The corresponding coefficients
% will be returned as zero.
constantPredictors = (range(X)==0);
ever_active = ~constantPredictors;

[~,P] = size(X);

% === 'Standardize' option ===

% Require a logical value.
if ~isscalar(standardize) || (~islogical(standardize) && standardize~=0 && standardize~=1)
    error(message('stats:lasso:InvalidStandardize'))
end

if ~withConstant && standardize
    warning(message('stats:lasso:StandardizeIgnored'));
end

weights = standardizeW(weights);
[X0,Y0,muX,sigmaX,muY] = standardizeXY(X, Y, weights, withConstant, standardize, constantPredictors);

% === 'Lambda' sequence ===

% lambdaMax is the penalty term (lambda) beyond which coefficients
% are guaranteed to be all zero.  If the command line does not provide
% a lambda sequence, we use lambdaMax in constructing the default 
% lambda sequence.  We always skip computation with lambda > lambdaMax
% because we know a priori that the computed coefficients will be zero.
%
% nullMSE is the mse of the fit using just a constant term.
% It is used to terminate the (ever-less penalized) fits when it becomes
% clear that we are overfitting.

[lambdaMax,nullMSE]=computeLambdaMax(X0,Y0,weights,alpha);

% Used with nullMSE (calculated below) to terminate
% (ever-less penalized) fits when overfitting is detected.
userSuppliedLambda = true;

if isempty(lambda)
    
    % Used with nullMSE (calculated below) to terminate 
    % (ever-less penalized) fits when overfitting is detected.
    userSuppliedLambda = false;
    
    % Sanity-check of 'NumLambda', should be positive integer.
    if ~isscalar(nLambda) || ~isreal(nLambda) || ~isfinite(nLambda) || nLambda < 1
        error(message('stats:lasso:InvalidNumLambda'));
    else
        nLambda = floor(nLambda);
    end
    
    % Sanity-checking of LambdaRatio, should be in [0,1).
    if ~isscalar(lambdaRatio) || ~isreal(lambdaRatio) || lambdaRatio <0 || lambdaRatio >= 1
        error(message('stats:lasso:InvalidLambdaRatio'));
    end
    
    if nLambda==1
        lambda = lambdaMax;
    else
        % Fill in a number "nLambda" of smaller values, on a log scale.
        if lambdaRatio==0
                lambdaRatio = LRdefault;
                addZeroLambda = true;
        else
            addZeroLambda = false;
        end
        lambdaMin = lambdaMax * lambdaRatio;
        loghi = log(lambdaMax);
        loglo = log(lambdaMin);
        lambda = exp(linspace(loghi,loglo,nLambda));
        if addZeroLambda
            lambda(end) = 0;
        else
            lambda(end) = lambdaMin;
        end
    end
    
else

    % Sanity check on user-supplied lambda sequence.  Should be non-neg real.
    if ~isreal(lambda) || any(lambda < 0)
        error(message('stats:lasso:NegativeLambda'));
    end

    lambda = sort(lambda(:),1,'descend');
    
end

% === 'RelTol' parameter ===
%
if ~isscalar(reltol) || ~isreal(reltol) || reltol <= 0 || reltol >= 1 || isnan(reltol)
    error(message('stats:lasso:InvalidRelTol'));
end

% === 'DFmax' parameter ===
%
% DFmax is #non-zero coefficients 
% DFmax should map to an integer in [1,P] but we truncate if .gt. P
%
if isempty(dfmax)
    dfmax = P;
else
    if ~isscalar(dfmax) || ~isreal(dfmax)
        error(message('stats:lasso:DFmaxBadType'));
    end
    try
        dfmax = uint32(dfmax);
    catch ME
        mm = message('stats:lasso:DFmaxBadType');
        throwAsCaller(MException(mm.Identifier,'%s',getString(mm)));
    end
    if dfmax < 1
        error(message('stats:lasso:DFmaxNotAnIndex'));
    else
        dfmax = min(dfmax,P);
    end
end

% === 'Mcreps' parameter ===
%
if ~isscalar(mcreps) || ~isreal(mcreps) || ~isfinite(mcreps) || mcreps < 1
    error(message('stats:lasso:MCRepsBadType'));
end
mcreps = fix(mcreps);

% === 'MaxIter' parameter ===
validateattributes(maxIter, {'numeric'},...
    {'scalar','positive','finite','integer'},...
    mfilename,'''MaxIter'' parameter');


% === 'CV' parameter ===
%
if isnumeric(cvp) && isscalar(cvp) && (cvp==round(cvp)) && (0<cvp)
    % cvp is a kfold value. Create a cvpartition to pass to crossval. 
    if (cvp > size(X,1))
        error(message('stats:lasso:InvalidCVforX'));
    end
    cvp = cvpartition(size(X,1),'Kfold',cvp);
elseif isa(cvp,'cvpartition')
    if strcmpi(cvp.Type,'resubstitution')
        cvp = 'resubstitution';
    elseif strcmpi(cvp.Type,'leaveout')
        error(message('stats:lasso:InvalidCVtype'));
    elseif strcmpi(cvp.Type,'holdout') && mcreps<=1
        error(message('stats:lasso:InvalidMCReps'));
    end
elseif strncmpi(cvp,'resubstitution',length(cvp))
    % This may have been set as the default, or may have been
    % provided at the command line.  In case it's the latter, we
    % expand abbreviations.
    cvp = 'resubstitution';
else
    error(message('stats:lasso:InvalidCVtype'));
end
if strcmp(cvp,'resubstitution') && mcreps ~= 1
    error(message('stats:lasso:InvalidMCReps'));
end

if isa(cvp,'cvpartition')
    if (cvp.N ~= size(X,1)) || (min(cvp.TrainSize) < 2)
        % We need partitions that match the total number of observations
        % (after stripping NaNs and Infs and zero observation weights), and
        % we need training sets with at least 2 usable observations.
        error(message('stats:lasso:TooFewObservationsForCrossval'));
    end
end

% === 'PredictorNames' parameter ===
%
% If PredictorNames is not supplied or is supplied as empty, we just 
% leave it that way. Otherwise, confirm that it is a cell array of strings.
%
if ~isempty(predictorNames) 
    if ~iscellstr(predictorNames) || length(predictorNames(:)) ~= size(X,2)
        error(message('stats:lasso:InvalidPredictorNames'));
    else
        predictorNames = predictorNames(:)';
    end
end

% === 'Options' parameter ===
% The 'Options' parameter is passed to crossval for handling.
% crossval will do sanity checking.

% --------------------
% Lasso model fits
% --------------------

% The struct 'stats' will comprise the second return argument.
% Put place holders for ever-present fields to secure the order
% we want in the struct.
stats = struct();
stats.Intercept      = [];
stats.Lambda         = [];
stats.Alpha          = alpha;
stats.DF             = [];
stats.MSE            = [];
stats.PredictorNames = predictorNames;

[B,Intercept,lambda,mse] = ...
    lassoFit(X,Y,X0,Y0,muX,sigmaX,muY,...
    weights,lambda,alpha,dfmax,standardize,reltol,lambdaMax,...
    ever_active,userSuppliedLambda,nullMSE,maxIter,withConstant);

% Store the number of non-zero coefficients for each lambda.
df = sum(B~=0,1);

% ---------------------------------------------------------
% If requested, use cross-validation to calculate 
% Prediction Mean Squared Error for each lambda.
% ---------------------------------------------------------

if ~isequal(cvp,'resubstitution')   
    % Replace dfmax with P, the number of predictors supplied at the
    % command line. dfmax might cause one fold to return empty values, 
    % because no lambda satisfies the dfmax criteria, while other folds 
    % return a numeric value. The lambda sequence has already been 
    % pruned by dfmax, if appropriate, in the call to lassoFit above.
    cvfun = @(Xtrain,Ytrain,Xtest,Ytest) lassoFitAndPredict( ...
        Xtrain,Ytrain,Xtest,Ytest,...
        lambda,alpha,P,standardize,reltol,ever_active,true,maxIter,withConstant,constantPredictors);
    if isempty(weights)
        weights = nan(size(X,1),1);
    end
    
    cvMSE = crossval(cvfun,[weights(:) X],Y, ...
        'Partition',cvp,'Mcreps',mcreps,'Options',ParOptions);
    mse = mean(cvMSE);
    se  = std(cvMSE) / sqrt(size(cvMSE,1));
    minMSE = min(mse);
    minIx = find(mse==minMSE,1);
    lambdaMin = lambda(minIx);
    minplus1 = mse(minIx) + se(minIx);
    seIx = find((mse(1:minIx) <= minplus1),1,'first');
    if isempty(seIx)
        lambdaSE = [];
    else
        lambdaSE = lambda(seIx);
    end
    
    % Deposit cross-validation results in struct for return value.
    stats.SE           = se;
    stats.LambdaMinMSE = lambdaMin;
    stats.Lambda1SE    = lambdaSE;
    stats.IndexMinMSE  = minIx;
    stats.Index1SE     = seIx;
end

% ------------------------------------------
% Order results by ascending lambda
% ------------------------------------------

nLambda = length(lambda);
reverseIndices = nLambda:-1:1;
lambda = lambda(reverseIndices);
lambda = reshape(lambda,1,nLambda);
B = B(:,reverseIndices);
Intercept = Intercept(reverseIndices);
df = df(reverseIndices);
mse = mse(reverseIndices);
if ~isequal(cvp,'resubstitution')
    stats.SE          = stats.SE(reverseIndices);
    stats.IndexMinMSE = nLambda - stats.IndexMinMSE + 1;
    stats.Index1SE    = nLambda - stats.Index1SE + 1;
end

stats.Intercept = Intercept;
stats.Lambda = lambda;
stats.DF = df;
stats.MSE = mse;
   
end % lasso

% ------------------------------------------
% SUBFUNCTIONS 
% ------------------------------------------

% ===================================================
%                 lassoFitAndPredict() 
% ===================================================

function mse = lassoFitAndPredict(Xtrain,Ytrain,Xtest,Ytest, ...
    lambda,alpha,dfmax,standardize,reltol,ever_active,userSuppliedLambda,maxIter,withConstant,constantPredictors)


trainWeights = Xtrain(:,1)';
if any(isnan(trainWeights))
    trainWeights = [];
end
trainWeights = standardizeW(trainWeights);
Xtrain = Xtrain(:,2:end);

trainWeights = standardizeW(trainWeights);
[Xtrain0,Ytrain0,muXtrain,sigmaXtrain,muYtrain] = standardizeXY(Xtrain, Ytrain, trainWeights, withConstant, standardize, constantPredictors);
[lambdaMaxtrain, nullMSEtrain] = computeLambdaMax(Xtrain0, Ytrain0, trainWeights, alpha);

[B,Intercept] = lassoFit(Xtrain,Ytrain,Xtrain0,Ytrain0,muXtrain,sigmaXtrain,muYtrain, ...
    trainWeights,lambda,alpha,dfmax,standardize,reltol,lambdaMaxtrain,ever_active,userSuppliedLambda,nullMSEtrain,maxIter,withConstant);

testWeights = Xtest(:,1)';
if any(isnan(testWeights))
    testWeights = ones(1,size(Xtest,1));
end
testWeights = standardizeW(testWeights);
Xtest = Xtest(:,2:end);
yfit = Intercept + Xtest*B;
mse = testWeights*((Ytest-yfit).^2) / sum(testWeights);
end

% ===================================================
%                 lassoFit() 
% ===================================================
function [B,Intercept,lambda,mspe] = ...
    lassoFit(X,Y,X0,Y0,muX,sigmaX,muY,weights,lambda,alpha,dfmax,standardize,reltol,lambdaMax,ever_active,userSuppliedLambda,nullMSE,maxIter,withConstant)
%
% ------------------------------------------------------
% Perform model fit for each lambda and the given alpha
% ------------------------------------------------------

[N,P] = size(X);
nLambda = length(lambda);

% If X has any constant columns, we want to exclude them from the
% coordinate descent calculations.  The corresponding coefficients
% will be returned as zero.
constantPredictors = (range(X)==0);
ever_active = ever_active & ~constantPredictors;

% === standardization and weights ===

observationWeights = ~isempty(weights);

% If using observation weights, make a weighted copy of the 
% predictor matrix, to save time in the weighted partial regressions.
if observationWeights
    wX0 = X0 .* weights';
    totalweight = cast(1,'like',X0);
else
    wX0 = X0;
    totalweight = cast(N,'like',X0);
end

% b will be the current coefficient estimate, iteratively updated.
% Because we retain b from one value of lambda to the next,
% we get a de facto warm start.
b = zeros(P,1,'like',X0);

r = Y0; % initial value of r = Y-X*0 = Y

% Record r and b used for threshold screen
b_ts_old = [];
r_ts = [];

% Preallocate the returned matrix of coefficients, B, and the intercepts.
B = zeros(P,nLambda,'like',X0);

active = false(1,P);

for i = 1:nLambda
    
    lam = lambda(i);
    if lam >= lambdaMax
        continue;
    end
    threshold = lam * alpha;
    
    % Denominator in coordinate descent update
    if withConstant
        if standardize
            if observationWeights
                % withConstant, standardize, observationWeights
                shrinkFactor = weights*(X0.^2) + lam*(1 - alpha);
            else % ~observationWeights
                % withConstant, standardize, ~observationWeights
                shrinkFactor = repmat(1 + lam*(1 - alpha),1,P);
            end
        else % ~standardize
            if observationWeights
                % withConstant, ~standardize, observationWeights
                shrinkFactor = weights*(X0.^2) + lam*(1 - alpha);
            else % ~observationWeights
                % withConstant, ~standardize, ~observationWeights
                shrinkFactor = (1/N) * sum(X0.^2) + lam*(1 - alpha);
            end
        end
    else % ~withConstant
        % standardize is ignored
        shrinkFactor = (1/N) * sum(X0.^2) + lam*(1 - alpha);
    end

    % Iterative coordinate descent until converged
    for numIter = 1:maxIter
        
        bold = b;
        
        [b,active,r] = cdescentCycle(X0,wX0,Y0, ...
            b,active,totalweight,shrinkFactor,threshold,r);
        
        if norm( (b-bold) ./ (1.0 + abs(bold)), Inf ) < reltol
            % Cycling over the active set converged.
            % Do one full pass through the predictors.
            % If there is no predictor added to the active set, we're done.
            % Otherwise, resume the coordinate descent iterations.
            bold = b;
            [potentially_active,b_ts_old,r_ts]  = ...
                thresholdScreen(X0,wX0,Y0,b,ever_active,threshold,b_ts_old,r_ts);
            if any(potentially_active)
                new_active = active | potentially_active;
                [b,new_active,r] = cdescentCycle(X0,wX0,Y0, ...
                    b,new_active,totalweight,shrinkFactor,threshold,r);
            else
                new_active = active;
            end

            if isequal(new_active, active)
                break
            else
                active = new_active;
            end
            
            if norm( (b-bold) ./ (1.0 + abs(bold)), Inf ) < reltol
                break
            end
        end

        if numIter == maxIter
            warning(message('stats:lasso:MaxIterReached',num2str(lam)));
        end
    end
    
    B(:,i) = b;
    
    % Halt if maximum model size ('DFmax') has been met or exceeded.
    if sum(active) > dfmax
        % truncate B and lambda output arguments
        lambda = lambda(1:(i-1));
        B = B(:,1:(i-1));
        break
    end
    
    % Halt if we have exceeded a threshold on the percent of
    % residual variance left unexplained.
    if ~userSuppliedLambda
        % Calculate mse of the current fit
        bsig = b ./ sigmaX';
        fit = (muY-muX*bsig) + X*bsig;
        residuals = Y - fit;
        if ~observationWeights
            mspe = mean(residuals.^2);
        else
            % This line relies on the weights having been normalized.
            mspe = weights * (residuals.^2);
        end
        if mspe < 1.0e-3 * nullMSE
            lambda = lambda(1:i);
            B = B(:,1:i);
            break
        end
    end
    
end % of lambda sequence

% ------------------------------------------
% Unwind the centering and scaling (if any)
% ------------------------------------------

B = B ./ sigmaX';
B(~ever_active,:) = 0;
Intercept = muY-muX*B;

% ------------------------------------------
% Calculate Mean Prediction Squared Error
% ------------------------------------------

fits = Intercept + X*B;
residuals = Y - fits;
if ~observationWeights
    mspe = mean(residuals.^2);
else
    % This line relies on the weights having been normalized.
    mspe = weights * (residuals.^2);
end

end %-lassoFit

% ===================================================
%                 cdescentCycle() 
% ===================================================

function [b,active,r] = cdescentCycle(X0, wX0, Y0, ...
    b, active, totalweight, shrinkFactor, threshold,r)
    [b,active,r] = internal.stats.lassoCoordDescentCycle(X0,wX0,Y0,b,active',totalweight,shrinkFactor',threshold,r); 
active = active';
end %-cdescentCycle

% ===================================================
%                 thresholdScreen() 
% ===================================================

function [potentially_active,b_ts_old,r] = thresholdScreen(X0, wX0, Y0, ...
    b, active, threshold, b_ts_old, r)
    if ~isempty(b_ts_old) && ~isempty(r)
        % Update r without calculating the entire Y0-X0*b
        delta_b = b - b_ts_old;
        active = (delta_b~=0)&active'; %transposition of active is done here
        r = internal.stats.lassoYminusXb( X0, Y0, delta_b, active, r );
    else
        % First time calculating r. Calculate the entire Y0-X0*b
        r = internal.stats.lassoYminusXb( X0, Y0, b, active' );
    end
    b_ts_old = b;
    potentially_active = abs(r' *wX0) > threshold;
end %-thresholdScreen

% ===================================================
%                 standardizeXY() 
% ===================================================

function [X0,Y0,muX,sigmaX,muY] = standardizeXY(X, Y, weights, withConstant, standardize, constantPredictors)
    
    observationWeights = ~isempty(weights);
    
    if withConstant
        if standardize
            if observationWeights
                % withConstant, standardize, observationWeights
                muX=weights*X;
                muY=weights*Y;
                sigmaX=sqrt(weights*((X-muX).^2));
                sigmaX(constantPredictors) = 1;
            else % ~observationWeights
                % withConstant, standardize, ~observationWeights
                muX=mean(X);
                muY=mean(Y);
                sigmaX=std(X,1);
                sigmaX(constantPredictors) = 1;
            end
        else % ~standardize
            if observationWeights
                % withConstant, ~standardize, observationWeights
                muX=weights*X;
                muY=weights*Y;
                sigmaX=1;
            else % ~observationWeights
                % withConstant, ~standardize, ~observationWeights
                muX=mean(X);
                muY=mean(Y);
                sigmaX=1;
            end
        end
    else % ~withConstant
        muX=zeros(1,size(X,2));
        muY=0;
        sigmaX=1;
    end
    X0 = (X-muX)./sigmaX;
    Y0 = Y-muY;
    
end %-standardizeXY

% ===================================================
%                 standardizeW() 
% ===================================================
function weights0 = standardizeW(weights)

    % Normalize weights up front.
    weights0 = weights / sum(weights);
    
    % Below, the convention is that weights is a row vector.
    weights0 = weights0(:)';
    
end %-standardizeW


% ===================================================
%                 computeLambdaMax() 
% ===================================================
function [lambdaMax,nullMSE]=computeLambdaMax(X0,Y0,weights,alpha)
    
    observationWeights = ~isempty(weights);
    N = size(X0,1);

    % Calculate max lambda that permits non-zero coefficients
    if ~observationWeights
        dotp = abs(X0' * Y0);
        lambdaMax = max(dotp) / (N*alpha);
    else
        wX0 = X0 .* weights';
        dotp = abs(sum(wX0 .* Y0));
        lambdaMax = max(dotp) / alpha;
    end
    
    if ~observationWeights
        nullMSE = mean(Y0.^2);
    else
        % This works because weights are normalized and Y0 is already
        % weight-centered.
        nullMSE = weights * (Y0.^2);
    end
    
end %-computeLambdaMax
