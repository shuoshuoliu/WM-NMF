function [Beta_hat,Psi_hat,stats,b_hat] = nlmefit(X,Y,grp,V,f,Beta_0,varargin)
%NLMEFIT Nonlinear mixed-effects estimation.
%   BETA = NLMEFIT(X,Y,GROUP,V,MODELFUN,BETA0) fits a nonlinear
%   mixed-effects regression model and returns estimates of the fixed
%   effects in BETA. By default, NLMEFIT fits a model where each model
%   parameter is the sum of a corresponding fixed and random effect, and
%   the covariance matrix of the random effects is diagonal, i.e.,
%   uncorrelated random effects.
%
%   X is an N-by-H matrix of N observations on H predictor variables. Y is
%   an N-by-1 vector of responses.
%
%   GROUP is a grouping variable indicating which of M groups each
%   observation belongs to. Valid GROUP values include a categorical
%   vector, a numeric vector, a string array or a cell array of strings.
%   See "help groupingvariable" for more information.
%
%   V is an M-by-G matrix of G group-specific predictor variables for each
%   of the M groups in the data. These are predictor values that take on
%   the same value for all observations in a group. Rows of V are ordered
%   according to GRP2IDX(GROUP). Use an M-by-G cell array for V if any of
%   the group- specific predictor values vary in size across groups.
%   Specify [] for V if there are no group predictors.
%
%   MODELFUN is a handle to a function that accepts predictor values and
%   model parameters, and returns fitted values. MODELFUN has the form
%
%      YFIT = MODELFUN(PHI,XFUN,VFUN) 
%
%   with input arguments
%
%      PHI    A 1-by-P vector of model parameters.
%      XFUN   An L-by-H array of predictor variables where L is 1 if XFUN is a
%             single row of X, NI if XFUN contains the rows of X for a single
%             group of size NI, or N if XFUN contains all rows of X.
%      VFUN   Either a 1-by-G vector of group-specific predictors for a single
%             group, corresponding to a single row of V; or an N-by-G matrix,
%             where if the K-th observation is in group I, then the K-th row
%             of VFUN is V(I,:). If V is empty, NLMEFIT calls MODELFUN with
%             only two inputs.
%
%   and returning an L-by-1 vector of fitted values YFIT. When either PHI or
%   VFUN contains a single row, that one row corresponds to all rows in the
%   other two input arguments. Note: for improved performance, use the
%   'Vectorization' parameter name/value pair (described below) if MODELFUN
%   can compute YFIT for more than one vector of model parameters in one call.
%
%   BETA0 is a F-by-1 vector with initial estimates for the F fixed effects.
%   By default, F is equal to the number of model parameters P.
%
%   NLMEFIT fits the model by maximizing an approximation to the marginal
%   likelihood, i.e., with the random effects integrated out, and assumes
%   that:
%      a) the random effects are multivariate normally distributed, and
%         independent between groups, and
%      b) the observation errors are independent, identically normally
%         distributed, and independent of the random effects. (However,
%         this assumption is changed by the ErrorModel parameter.)
%
%   [BETA,PSI] = NLMEFIT(...) returns PSI, an R-by-R estimated covariance
%   matrix for the random effects. By default, R is equal to the number of
%   model parameters P.
%
%   [BETA,PSI,STATS] = NLMEFIT(...) returns STATS, a structure with fields:
%       logl        The maximized log-likelihood for the fitted model
%       rmse        The root mean squared residual (computed on the log
%                   scale for the 'exponential' error model)
%       errorparam  The estimated parameters of the error variance model
%       aic         The Akaike information criterion
%       bic         The Bayesian information criterion
%       sebeta      The standard errors for BETA
%       covb        The estimated covariance of the parameter estimates
%       dfe         The error degrees of freedom
%       pres        The population residuals (y - y_population), where
%                   y_population is the population predicted values
%       ires        The individual residuals (y - y_individual), where
%                   y_individual is the individual predicted values
%       pwres       The population weighted residuals
%       cwres       The conditional weighted residuals
%       iwres       The individual weighted residuals
%
%   [BETA,PSI,STATS,B] = NLMEFIT(...) returns B, an R-by-M matrix of estimated
%   random effects for the M groups. By default, R is equal to the number of
%   model parameters P.
%
%   [...] = NLMEFIT(X,Y,GROUP,V,FUN,BETA0,'param1',val1,...) specifies
%   additional parameter name/value pairs that allow you to define the model
%   and control the estimation algorithm, as described below.
%
%   By default, NLMEFIT fits a model where each model parameter is the sum of
%   a corresponding fixed and random effect. Use the following parameter
%   name/value pairs to fit a model with a different number of or dependence
%   on fixed or random effects. Use at most one parameter name with an 'FE'
%   prefix and one parameter name with an 'RE' prefix. Note that some choices
%   change the way NLMEFIT calls MODELFUN, as described further below.
%
%       'FEParamsSelect' A vector specifying which elements of the model
%                        parameter vector PHI include a fixed effect, as a
%                        numeric vector with elements in 1:P, or as a 1-by-P
%                        logical vector.  The model will include F fixed
%                        effects, where F is the specified number of elements.
%       'FEConstDesign'  A P-by-F design matrix ADESIGN, where ADESIGN*BETA
%                        are the fixed components of the P elements of PHI.
%       'FEGroupDesign'  A P-by-F-by-M array specifying a different P-by-F
%                        fixed effects design matrix for each of the M groups.
%       'FEObsDesign'    A P-by-F-by-N array specifying a different P-by-F
%                        fixed effects design matrix for each of the N
%                        observations.
%
%       'REParamsSelect' A vector specifying which elements of the model
%                        parameter vector PHI include a random effect, as a
%                        numeric vector with elements in 1:P, or as a 1-by-P
%                        logical vector.  The model will include R random
%                        effects, where R is the specified number of elements.
%       'REConstDesign'  A P-by-R design matrix BDESIGN, where BDESIGN*B are
%                        the random components of the P elements of PHI.
%       'REGroupDesign'  A P-by-R-by-M array specifying a different P-by-R
%                        random effects design matrix for each of M groups.
%       'REObsDesign'    A P-by-R-by-N array specifying a different P-by-R
%                        random effects design matrix for each of N
%                        observations.
%
%   The default model is equivalent to setting both 'FEConstDesign' and
%   'REConstDesign' to EYE(P), or to setting both 'FEParamsSelect' and
%   'REParamsSelect' to 1:P.
%
%   Additional optional parameter name/value pairs control the iterative
%   algorithm used to maximize the likelihood:
%
%       'ApproximationType' The method used to approximate the non-linear
%                           mixed effects model likelihood:
%              'LME'    Use the likelihood for the linear mixed-effects
%                       model at the current conditional estimates of BETA
%                       and B. This is the default.
%              'RELME'  Use the restricted likelihood for the linear
%                       mixed-effects model at the current conditional
%                       estimates of BETA and B.
%              'FO'     First order (Laplacian) approximation without
%                       random effects.
%              'FOCE'   First order (Laplacian) approximation at the
%                       conditional estimates of B.
%
%       'CovParameterization'  Specifies the parameterization used
%                              internally for the scaled covariance matrix
%                              (PSI/sigma^2). 'chol' for the Cholesky
%                              factorization, or 'logm' (the default) for
%                              the Cholesky factorization of the matrix
%                              logarithm.
%
%       'CovPattern'         Specifies an R-by-R logical or numeric matrix
%                            PAT that defines the pattern of the random
%                            effects covariance matrix PSI. NLMEFIT
%                            computes estimates for the variances along the
%                            diagonal of PSI as well as covariances that
%                            correspond to non-zeroes in the off-diagonal
%                            of PAT.  NLMEFIT constrains the remaining
%                            covariances, i.e., those corresponding to
%                            off-diagonal zeroes in PAT, to be zero. PAT
%                            must be a row-column permutation of a block
%                            diagonal matrix, and NLMEFIT adds non-zero
%                            elements to PAT as needed to produce such a
%                            pattern. The default value of PAT is EYE(R),
%                            corresponding to uncorrelated random effects.
%
%                            Alternatively, specify PAT as a 1-by-R vector
%                            containing values in 1:R. In this case,
%                            elements of PAT with equal values define
%                            groups of random effects, NLMEFIT estimates
%                            covariances only within groups, and constrains
%                            covariances across groups to be zero.
%
%       'OptimFun'  Either 'fminsearch' or 'fminunc', specifying the
%                   optimization function to be used in maximizing the
%                   likelihood.  Default is 'fminsearch'.  You may only
%                   specify 'fminunc' if Optimization Toolbox is available.
%
%       'Options'  A structure created by a call to STATSET. NLMEFIT uses the
%                  following STATSET parameters:
%            'TolX'         Termination tolerance on the estimated fixed
%                           and random effects. Default 1e-4.
%            'TolFun'       Termination tolerance on the log-likelihood
%                           function. Default 1e-4.
%            'MaxIter'      Maximum number of iterations allowed. Default 200.
%            'Display'      Level of display during estimation. 'off' (the
%                           default) displays no information, 'final'
%                           displays information after the final iteration
%                           of the estimation algorithm, 'iter' displays
%                           information at each iteration.
%            'DerivStep'   - Relative difference used in finite difference gradient
%                            calculation.  May be a scalar, or a vector
%                            with the number of elements equal to the
%                            number of model parameters P.  Defaults to
%                            EPS^(1/3).
%            'FunValCheck'  'on' (the default) to check for invalid values
%                           (such as NaN or Inf) from MODELFUN, or 'off' to
%                           skip this check.
%            'OutputFcn'    Function handle specified using @, a cell array
%                           with function handles or an empty array
%                           (default). NLMEFIT calls all output functions
%                           after each iteration. See NLMEFITOUTPUTFCN for
%                           an example of an output function.
%
%       'ParamTransform' A vector of P values specifying a transformation
%                        function f() for each of the P parameters:
%                            XB = ADESIGN*BETA + BDESIGN*B
%                            PHI = f(XB)
%                        Each element of the vector must be one of the
%                        following integer codes specifying the
%                        transformation for the corresponding value of PHI:
%                             0: PHI = XB  (default for all parameters)
%                             1: log(PHI) = XB
%                             2: probit(PHI) = XB
%                             3: logit(PHI) = XB
%
%       'ErrorModel' A string specifying the form of the error term.
%                    Default is 'constant'. Each model defines the error
%                    using a standard normal (Gaussian) variable e, the
%                    function value f, and one or two parameters a and b.
%                    Choices are:
%                       'constant'         y = f + a*e
%                       'proportional'     y = f + b*f*e
%                       'combined'         y = f + (a+b*abs(f))*e
%                       'exponential'      y = f*exp(a*e), or equivalently
%                                          log(y) = log(f) + a*e
%                    If this parameter is given, the output STATS.errorparam
%                    field has the value
%                        a       for 'constant' and 'exponential'
%                        b       for 'proportional'
%                        [a b]   for 'combined'
%
%       'RefineBeta0'  Determines whether NLMEFIT will make an initial
%                      refinement of BETA0 by fitting the model defined by
%                      MODELFUN without random effects. Default is 'on'.
%
%       'Vectorization'  Determines the possible sizes of the PHI, XFUN,
%                        and VFUN input arguments to MODELFUN.  Possible
%                        values are:
%             'SinglePhi'    MODELFUN is a function (such as an ODE solver)
%                            that can only compute YFIT for a single set of
%                            model parameters at a time, i.e., PHI must be
%                            a single row vector in each call. NLMEFIT
%                            calls MODELFUN in a loop if necessary using a
%                            single PHI vector and with XFUN containing
%                            rows for a single observation or group at a
%                            time. VFUN may be a single row that applies to
%                            all rows of XFUN, or a matrix with rows
%                            corresponding to rows in XFUN.
%            'SingleGroup'   MODELFUN can only accept inputs corresponding
%                            to a single group in the data, i.e., XFUN must
%                            contain rows of X from a single group in each
%                            call. Depending on the model, PHI is a single
%                            row that applies to the entire group, or a
%                            matrix with one row for each observation. VFUN
%                            is a single row.
%            'Full'          MODELFUN can accept inputs for multiple
%                            parameter vectors and multiple groups in the
%                            data. Either PHI or VFUN may be a single row
%                            that applies to all rows of XFUN, or a matrix
%                            with rows corresponding to rows in XFUN. Using
%                            this option can improve performance by
%                            reducing the number of calls to MODELFUN, but
%                            may require MODELFUN to perform singleton
%                            expansion on PHI or V.
%                        The default for 'Vectorization' is 'SinglePhi'. In
%                        all cases, if V is empty, NLMEFIT calls MODELFUN
%                        with only two inputs.
%
%   Example:
%      % Fit a model to data on concentrations of the drug indomethacin in
%      % the bloodstream of six subjects over eight hours
%      load indomethacin
%      model = @(phi,t)(phi(:,1).*exp(-phi(:,2).*t) + phi(:,3).*exp(-phi(:,4).*t));
%      phi0 = [1 1 1 1];
%      xform = [0 1 0 1]; % log transform for 2nd and 4th parameters
%      [beta,PSI,stats,br] = nlmefit(time,concentration,subject,[],...
%                                    model,phi0, 'ParamTransform',xform);
%
%      % Plot the data along with an overall "population" fit
%      phi = [beta(1), exp(beta(2)), beta(3), exp(beta(4))];
%      h = gscatter(time,concentration,subject);
%      xlabel('Time (hours)')
%      ylabel('Concentration (mcg/ml)')
%      title('{\bf Indomethacin Elimination}')
%      xx = linspace(0,8);
%      line(xx,model(phi,xx),'linewidth',2,'color','k')
%
%      % Plot individual curves based on random effect estimates
%      for j=1:6
%          phir = [beta(1)+br(1,j), exp(beta(2)+br(2,j)), ...
%                  beta(3)+br(3,j), exp(beta(4)+br(4,j))];
%          line(xx,model(phir,xx),'color',get(h(j),'color'))
%      end
%      legend(h)
%
%   See also NLINFIT, NLMEFITOUTPUTFCN, NLMEFITSA, GROUPINGVARIABLE.

%   Copyright 2008-2016 The MathWorks, Inc.

%   NLMEFIT calls MODELFUN with PHI, XFUN, and VFUN having the following sizes
%   under the following conditions:
%
%   1a. If 'Vectorization' is 'SinglePhi', and neither 'REObsDesign' nor
%       'FEObsDesign' is specified, NLMEFIT calls MODELFUN once for each group
%       of observations:
%           PHI is 1-by-P, XFUN is NI-by-H, VFUN is 1-by-G
%       where NI is the number of observations in the I-th group.
%
%   1b. If 'Vectorization' is 'SinglePhi', and 'REObsDesign' or 'FEObsDesign'
%       is specified, NLMEFIT calls MODELFUN separately for each observation:
%           PHI is 1-by-P, XFUN is 1-by-H, VFUN is 1-by-G
%
%   2.  If 'Vectorization' is 'SingleGroup', NMLMEFIT calls MODELFUN once for
%       each group of observations:
%           PHI is 1-by-P or NI-by-P, XFUN is NI-by-H, VFUN is 1-by-G
%       where NI is the number of observations in the I-th group, and PHI's
%       size depends on whether or not 'REObsDesign' or 'FEObsDesign' is
%       specified.
%
%   3.  If 'Vectorization' is 'Full', NLMEFIT may call MODELFUN in one of two
%       ways.  First, once for each group of observations:
%           PHI is 1-by-P or NI-by-P, XFUN is NI-by-H, VFUN is 1-by-G
%       where NI is the number of observations in the I-th group, and PHI's
%       size depends on whether or not 'REObsDesign' or 'FEObsDesign' is
%       specified.  Second, once for all observations:
%           PHI is N-by-P, XFUN is N-by-H, VFUN is N-by-G
%       where if the K-th observation is in group I, then the K-th row of PHI
%       is the vector of model parameters that contains random effects for
%       group I. Similarly, VFUN is an expanded version of V, so that if the
%       K-th observation is in group I, the K-th row of VFUN is V(I,:).
%
%   4.  If 'Vectorization' is 'SinglePhi' or 'Full', NLMEFIT also calls
%       MODELFUN under some circumstances with all observations:
%           PHI is 1-by-P, XFUN is N-by-H, VFUN is N-by-G
%       NLMEFIT calls MODELFUN this way when fitting without random effects in
%       the initial refinement of BETA0 when 'RefineBeta0' is 'on', or when
%       computing the first-order approximation of the likelihood, i.e., when
%       'ApproximationType' is 'FO'.  VFUN is an expanded version of V, so
%       that if the K-th observation is in group I, the K-th row of VFUN is
%       V(I,:).
%
%   In all cases, if V is empty, NMLEFIT calls MODELFUN with only two input
%   arguments.  MODELFUN must return YFIT with length corresponding to the
%   number of rows in XFUN.


%% Get and validate mandatory inputs
if nargin > 2
    grp = convertStringsToChars(grp);
end

if nargin > 6
    [varargin{:}] = convertStringsToChars(varargin{:});
end

[N,NP] = size(X);   % N -> number of observations,   NP -> number of predictors in X
if N<2 || NP==0
    error(message('stats:nlmefit:invalidX'))
end
X = double(X);
if (size(Y,1)~=N) || (size(Y,2)~=1)
    error(message('stats:nlmefit:invalidY'))
end
Y = double(Y);
if any(isnan(Y)) || any(isnan(X(:)))
    error(message('stats:nlmefit:noSupportNaNs'))
end
if isa(grp,'categorical')
    % remove levels that do not appear in data
    grp = removecats(grp);
end
[grp,GN] = grp2idx(grp);
M = numel(GN); % M = number of groups
if M<2
    error(message('stats:nlmefit:oneGROUP'))
end
if (size(grp,1)~=N) || (size(grp,2)~=1)
    error(message('stats:nlmefit:invalidGROUP'))
end
% A lot of searching for matching groups is required. To optimize
% performance, only do this search once.
grpMatch = bsxfun(@eq, grp, 1:M);

if isempty(V)
    V = zeros(M,0);
elseif size(V,1)~=M
    error(message('stats:nlmefit:invalidV'))
end
if ~isa(f,'function_handle')
    error(message('stats:nlmefit:invalidFUN'))
end
if ~isvector(Beta_0)
    error(message('stats:nlmefit:invalidBETA0'))
end
Beta_0 = Beta_0(:);
q = numel(Beta_0);    % q -> fixed effects (size of Beta)

%% Parse optional parameter name-value pair input arguments
[~,r,A,B,patternCov,dofCov,refineBeta0,ApproximationType,...
    modelVectorization,parType,Options,minMethod,outputFcn,...
    paramTransform,errorModel,errorModelParam] = ...
    parseInVarargin(q,M,N,varargin{:});
numParam = dofCov+q+numel(errorModelParam);
haveOutputFcn = ~isempty(outputFcn);
stop = false;

if strcmp(errorModel, 'exponential')
    % Exponential, y = f*exp(a*e), or log(y) = log(f) + a*e
    if ~all(Y>0)
        error(message('stats:nlmefit:PositiveYRequired'))
    else
        Y = log(max(Y,realmin));
    end
end


% The first iteration of the Alternating Algorithm (or the Laplace
% approximation algorithm in case alternatingMethod==false) assumes a
% constant variance error model:
weights = ones(N,1);

% Extracting algorithm parameters from STATSET structure
tolF = max(Options.TolFun,sqrt(eps(X(1))));
tolX = max(Options.TolX,sqrt(eps(Beta_0(1))));
fdiffstep = Options.DerivStep;
if isscalar(fdiffstep)
    % Expand to fdiffstep to be a row vector of the number of parameters,
    % which is equal to the number of rows in A.
    fdiffstep = repmat(fdiffstep, 1, size(A, 1));
end
% Ensure fdiffstep is a row vector.
fdiffstep = fdiffstep(:)';

maxIter = Options.MaxIter;
switch Options.Display
    case 'none',   verbose = 0;
    case 'off',    verbose = 0;
    case 'final',  verbose = 2;
    case 'iter',   verbose = 3;
    case 'on',     verbose = 3;
    otherwise,     verbose = 0;
end

switch ApproximationType
    case 'LME'
        alternatingMethod = true;
        restricted = false;
        laplacianMethod = false;
        foceMethod = false;
    case 'RELME'
        alternatingMethod = true;
        restricted = true;
        laplacianMethod = false;
        foceMethod = false;
    case 'FO'
        alternatingMethod = true;
        restricted = false;
        laplacianMethod = true;
        foceMethod = false;
    case 'FOCE'
        alternatingMethod = true;
        restricted = false;
        laplacianMethod = true;
        foceMethod = true;
end

% Options not documented and not supported
DeltaInitFrac = 0.375;
funValCheckInMin = 'off';

% Setup function handles
[computeFAllGroupsFcn, computeFByGroupFcn, getFCountFcn, ...
    computeJacobianAllGroupsFcn, computeJacobianByGroupFcn] ...
    = getFcns(X,V,A,B,f,grp,grpMatch,modelVectorization,paramTransform,fdiffstep,errorModel);

%% Refine the values of Beta without random effects
Beta_hat = Beta_0;
if refineBeta0
    if verbose>2
        disp(getString(message('stats:nlmefit:WithoutRandomEffects')))
        dispVal('Beta_0''',Beta_hat')
        t = cputime;
    end
    
    LMfitOptions = struct('TolFun',tolF,'TolX',tolX,'FunValCheck',Options.FunValCheck,'MaxIter',maxIter*q);
    % Create function handles as a function of Beta only
    bDummy = zeros(r,M);
    model = @(Beta) computeFAllGroupsFcn(Beta, bDummy);
    % First output from computeJacobianAllGroupsFcn is JBeta.
    jacobian = @(f0,Beta) computeJacobianAllGroupsFcn(f0, Beta, bDummy);
    Beta_hat = LMfit(Y,Beta_hat,LMfitOptions,model,jacobian);
    
    if verbose>2
        dispVal('Beta_hat''',Beta_hat')
        fprintf('%s',getString(message('stats:nlmefit:ModelEvaluationsInSeconds',getFCountFcn(),sprintf('%f',cputime-t))));
    end
else
    if verbose>2
        disp(getString(message('stats:nlmefit:WithoutRandomEffectsOFF')))
        disp(' ')
    end
end

%% Safe initialization of random effects
% Get safe initial value for Delta (Bates&Pinheiro 1998, page 12)
if verbose>2
    disp(getString(message('stats:nlmefit:InitializeRandomEffects')))
end
t = cputime;
b_hat = zeros(r,M);
f0 = computeFAllGroupsFcn(Beta_hat,b_hat);
[~, Zw] = computeJacobianAllGroupsFcn(f0, Beta_hat, b_hat);

Delta = diag(DeltaInitFrac*sqrt(sum(Zw.^2,1)/M));
theta_hat = Delta2theta(Delta,r,parType,patternCov);

if verbose>2
    dispVal('Delta_0',Delta)
    dispVal('theta_0''',theta_hat')
    fprintf('%s',getString(message('stats:nlmefit:ModelEvaluationsInSeconds',getFCountFcn(),sprintf('%f',cputime-t))))
end

%% Estimate by the Alternating algorithm (Pinheiro&Bates, page 313)

if alternatingMethod
    
    % take initial values from last estimation or input guesses
    Beta = Beta_hat;
    b = b_hat;
    theta = theta_hat;
    Delta = theta2Delta(theta,r,parType,patternCov);
    if verbose>2
        if restricted
            disp(getString(message('stats:nlmefit:AlternatingAlgorithmRestricted')))
        else
            disp(getString(message('stats:nlmefit:AlternatingAlgorithmNonrestricted')))
        end
        dispVal('Beta_0''',Beta')
        dispVal('Delta_0',Delta)
        dispVal('theta_0''',theta')
    end
    
    t = cputime;
    
    % Algorithm constants
    tolFac = 2;
    
    PNLSOpt = struct('TolFun',tolF*(10^tolFac),'TolX',tolX*(10^tolFac),...
        'FunValCheck',Options.FunValCheck,...
        'MaxIter',(r*M+q)*maxIter,'OutputFcn',[]);
    
    minMethodOpt = optimset(...
        'TolFun',tolF*(10^tolFac),...
        'TolX',tolX*(10^tolFac),...
        'MaxFunEvals',inf,...
        'MaxIter',(dofCov+q)*maxIter,...
        'FunValCheck',funValCheckInMin,...
        'Display','none',...
        'LargeScale','off',... % Used by fminunc
        'OutputFcn',[]);
    
    % Precalculate the constant part of the loglikelihood
    if restricted
        Klike = (log(N-q)-log(2*pi)-1)*(N-q)/2; % constant part of the loglikelihood
    else
        Klike = (log(N)-log(2*pi)-1)*N/2; % constant part of the loglikelihood
    end
    if strcmp(errorModel, 'exponential')
        Klike = Klike - sum(Y);
    end
    
    LMEnegloglikelihoodFcnFcn = getLMEnegloglikelihoodFcnFcn(r,q,M,N,grpMatch,restricted,parType,patternCov);
    
    nIter = 0;
    nlALT = inf;
    sigma2_hat = NaN;
    llike = NaN;
    
    initStop = false;
    if haveOutputFcn
        ofopt = struct('procedure','ALT','iteration',nIter,'inner',...
            struct('procedure','none','state','none','iteration',NaN),...
            'fval',llike,'Psi', Delta2Psi_sigma2(Delta,r)*sigma2_hat,...
            'theta',theta,'mse',sigma2_hat,'caller','nlmefit');
        initStop = callOutputFcns(outputFcn,Beta,ofopt,'init');
    end
    
    while nIter < maxIter % Alternating algorithm (Pinheiro&Bates, page 313)
        nIter = nIter+1;
        % PNLS step: The precision factor (Delta) is fixed, find the
        % conditional modes of the random effects (b) and the conditional
        % estimates of the fixed effects (Beta) by minimizing a penalized
        % non-linear least squares objective function:
        
        if haveOutputFcn
            ofopt = struct('procedure','ALT','iteration',nIter-1,'inner',struct,...
                'fval',llike,'Psi',Delta2Psi_sigma2(Delta,r)*sigma2_hat,...
                'theta',theta,'mse',sigma2_hat,'caller','nlmefit');
            PNLSOpt.OutputFcn = @(Beta,inneriter,innerstate) pnlsOutputFcn(outputFcn,Beta,ofopt,inneriter,innerstate);
        end
        
        [BetaNew,bNew,stop] = GNfitPNLS(Delta,Y,Beta,b,weights,PNLSOpt,grp,grpMatch,M,N,r,q,computeFAllGroupsFcn,computeJacobianAllGroupsFcn);
        
        MaxBetabDiff = max([abs(BetaNew - Beta); abs(bNew(:)-b(:))]);
        Beta = BetaNew;
        b = bNew;
        
        % LME step: Updates the estimate of the precision factor (Delta) by
        % linearizing around the conditional modes of the random effects (b)
        % and the conditional estimates of the fixed effects (Beta):
        
        % LINEARIZATION: (Eq.(7.11) Pinheiro&Bates, page 313) (all groups at
        % once), new linear model is: ww_i = Xw_i Beta + Zw_i b_i + e_i
        fBetab = computeFAllGroupsFcn(Beta,b);
        [Xw, Zw] = computeJacobianAllGroupsFcn(fBetab, Beta, b);
        ww = Y-fBetab+Xw*Beta+sum(Zw.*b(:,grp)',2);
        Xwww = [Xw ww];
        
        % The derivative matrices and working vector are weighted by the error
        % model:
        Xwww = bsxfun(@times,weights,Xwww);
        Zw = bsxfun(@times,weights,Zw);
        
        % "Burn-in" Zw and Xwww
        fcn = LMEnegloglikelihoodFcnFcn(Zw,Xwww);
        
        % Minimize the negative profiled loglikelihood of the LME model
        % (Eqs. (2.20 or 2.23) Pinheiro&Bates, pp 71-76):
        if haveOutputFcn
            if stop || initStop
                break;
            end
            ofopt = struct('procedure','ALT','iteration',nIter-1,'inner',struct,...
                'fval',llike,'Psi',Delta2Psi_sigma2(Delta,r)*sigma2_hat,...
                'theta',theta,'mse',sigma2_hat,'caller','nlmefit');
            minMethodOpt.OutputFcn = @(optimX,optimStruc,optimState) lmeOutputFcn(optimX,optimStruc,optimState,outputFcn,ofopt,Beta,Klike);
        end
        
        
        [theta_hat,nlLME,exitFlag] = robustMinimizer(fcn,theta,minMethod,minMethodOpt);
        
        if haveOutputFcn && exitFlag==-1
            stop = true;
            break;
        elseif exitFlag==0
            warning(message('stats:nlmefit:IterationLimitExceededLME'));
        elseif exitFlag==-3
            warning(message('stats:nlmefit:ObjectiveLimitExceededLME'));
        end
        
        % Update weights
        switch errorModel
            case 'constant'
                % weights = ones(N,1);
            case 'proportional'
                weights = 1./abs(fBetab);
                weights(isinf(weights)) = 1;
            case 'combined'
                ab = abs(fminsearch(@(ab) error_ab(ab,Y,fBetab),errorModelParam));
                weights = 1./abs(ab(1)+ab(2)*abs(fBetab));
                weights(isinf(weights)) = 1;
                errorModelParam = ab;
            case 'exponential'
                % weights = ones(N,1);
        end
        
        theta = theta_hat;
        Delta = theta2Delta(theta,r,parType,patternCov);
        
        nlALTDiff = abs(nlALT-nlLME);
        nlALT = nlLME;
        
        if haveOutputFcn % call outer level iteration outputFcn
            [nlLME,c_1]= fcn(theta);
            llike = Klike - nlLME; % constant part of the loglikelihood + lALT
            sigma2_hat = (c_1.^2)/N;
            ofopt = struct('procedure','ALT','iteration',nIter-1,'inner',...
                struct('procedure','none','state','none','iteration',NaN),...
                'fval',llike,'Psi', Delta2Psi_sigma2(Delta,r)*sigma2_hat,...
                'theta',theta,'mse',sigma2_hat,'caller','nlmefit');
            stop = callOutputFcns(outputFcn,Beta,ofopt,'iter');
            if stop, break; end
        end
        
        % Check output conditions for the Alternating algorithm
        if nlALTDiff <= tolF*(10.^tolFac) &&  MaxBetabDiff<= tolX*(10.^tolFac)
            if tolFac <= 0
                break
            else
                tolFac = tolFac - 1;
                LMEminMethodOpt.TolFun = tolF*(10.^tolFac);
                LMEminMethodOpt.TolX = tolX*(10.^tolFac);
                PNLSOpt.TolFun = tolF*(10.^tolFac);
                PNLSOpt.TolX = tolX*(10.^tolFac);
            end
        end
        
    end
    
    if nIter>=maxIter
        warning(message('stats:nlmefit:IterationLimitExceededAlt'));
    end
    
    [nlLME,c_1,R00]= fcn(theta);
    llike = Klike - nlLME; % constant part of the loglikelihood + lALT
    llike = llike + sum(log(weights));  % including error model weights into ll
    
    Beta_hat = Beta;                    % use Beta computed in the PNLS step
    b_hat = b;                          % use b computed in the PNLS step
    Psi_sigma2 = Delta2Psi_sigma2(Delta,r);  % use Delta computed in the LME step
    
    sigma2_hat = (c_1.^2)/N;            % use sigma computed in the LME step
    Psi_hat = Psi_sigma2 * sigma2_hat;
    theta_hat = theta;
    Delta_hat = Delta;
    if (verbose>2) || (~laplacianMethod&&verbose)
        dispVal('Beta_hat''',Beta_hat')
        dispVal('theta_hat''',theta_hat')
        dispVal('logLike',llike)
        dispVal('sigma2_hat',sigma2_hat)
        dispVal('Psi_hat',Psi_hat,'%9.6f')
        dispVal('Delta_hat',Delta_hat,'%9.6f')
        dispVal('Deg.Freedom',numParam,'%d')
        dispVal('AIC',-2*llike+2*numParam)
        dispVal('BIC',-2*llike+log(M-q*restricted)*numParam)
        fprintf('%s',getString(message('stats:nlmefit:IterationsAndModelSeconds',nIter,getFCountFcn(),sprintf('%f',cputime-t))))
        % Optional values that the user may want to display (when debugging)
        % dispVal('b_hat''',b_hat','%6.3f')
        % dispVal('log(sigma2_hat)',log(sigma2_hat))
        % dispVal('log(diag(L))',log(diag(chol(Psi_hat./sigma2_hat)))')
        % dispVal('triu(L,1)',nonzeros(triu(chol(Psi_hat./sigma2_hat),1))')
        % disp(' ')
    end
    
else
    nIter = 0;
    if verbose>2
        disp(getString(message('stats:nlmefit:AlternatingAlgorithmOFF')))
        disp(' ')
    end
end

ALTniter = nIter;

%% Estimate by the Laplace approximation algorithm

if laplacianMethod && ~(haveOutputFcn && stop)
    
    Beta = Beta_hat;
    b = b_hat;
    theta = theta_hat;
    Delta = theta2Delta(theta,r,parType,patternCov);
    
    % Algorithm constants
    tolFac = 2;
    
    if verbose>2
        disp(getString(message('stats:nlmefit:LaplacianApproximation')))
        dispVal('Beta_0''',Beta')
        dispVal('Delta_0',Delta)
        dispVal('theta_0''',theta')
    end
    t = cputime;
    
    GPNLSOpt = struct('TolFun',tolF*(10^tolFac),'TolX',tolX*(10^tolFac),...
        'FunValCheck',Options.FunValCheck,'MaxIter',r*maxIter,...
        'OutputFcn',[]);
    
    % Precalculate the constant part of the loglikelihood
    Klike = -N/2 * (1+log(2*pi));
    if strcmp(errorModel, 'exponential')
        Klike = Klike - sum(Y);
    end
    
    LAnegloglikelihoodFcnFcn = getLAnegloglikelihoodFcnFcn(Y,grpMatch,q,r,dofCov,M,N,parType,patternCov,computeFAllGroupsFcn,computeJacobianAllGroupsFcn);
    
    nIter = 0;
    nlLA = inf;
    sigma2_hat = NaN;
    llike = Klike - nlLA;
    
    if ~stop
        while nIter < maxIter % Laplacian approximation (Pinheiro&Bates, page 315)
            nIter = nIter+1;
            
            % PNLS step: The precision factor (Delta) and the conditional
            % estimates of the fixed effects (Beta) are fixed, find the conditional
            % modes of the random effects (b) by minimizing M penalized
            % non-linear least squares objective functions:
            
            bNew = zeros(r,M);
            if foceMethod
                
                if haveOutputFcn % call outer level iteration outputFcn
                    ofopt = struct('procedure','LAP','iteration',ALTniter+nIter-1,'inner',...
                        struct('procedure','PNLS','state','init','iteration',0),...
                        'fval',llike,'Psi', Delta2Psi_sigma2(Delta,r)*sigma2_hat,...
                        'theta',theta,'mse',sigma2_hat,'caller','nlmefit');
                    stop = callOutputFcns(outputFcn,Beta,ofopt,'iter');
                    if stop, break; end
                end
                
                for i = 1:M
                    bNew(:,i) = GNfitGPNLS(Delta,Y(grpMatch(:,i)),weights(grpMatch(:,i)),Beta,b(:,i),GPNLSOpt,i,r,computeFByGroupFcn,computeJacobianByGroupFcn);
                end
                
                if haveOutputFcn % call outer level iteration outputFcn
                    ofopt = struct('procedure','LAP','iteration',ALTniter+nIter-1,'inner',...
                        struct('procedure','PNLS','state','done','iteration',M),...
                        'fval',llike,'Psi', Delta2Psi_sigma2(Delta,r)*sigma2_hat,...
                        'theta',theta,'mse',sigma2_hat,'caller','nlmefit');
                    stop = callOutputFcns(outputFcn,Beta,ofopt,'iter');
                    if stop, break; end
                end
            end
            MaxbDiff = max(abs(bNew(:)-b(:)));
            b = bNew;
            
            % Maximize Likelihood Approximation step: The random effects (b) are
            % fixed. Find the estimates of the fixed effects (Beta) and the
            % precision factor (Delta) by minimizing the negative profiled
            % loglikelihood of the Laplacian approximation profiled on sigma2 given
            % by Eq.(7.19) in Pinheiro&Bates (page 319).
            fcn = LAnegloglikelihoodFcnFcn(b,weights);
            
            params0 = [Beta;theta];
            
            
            if haveOutputFcn
                ofopt = struct('procedure','LAP','iteration',ALTniter+nIter-1,'inner',struct,...
                    'fval',llike,'Psi',Delta2Psi_sigma2(Delta,r)*sigma2_hat,...
                    'theta',theta,'mse',sigma2_hat,'caller','nlmefit');
                minMethodOpt.OutputFcn = @(optimX,optimStruc,optimState) plmOutputFcn(optimX,optimStruc,optimState,outputFcn,ofopt,Klike,q);
            end

            [params_hat,nlLAnew,exitFlag] = robustMinimizer(fcn,params0,minMethod,minMethodOpt);
           
            if haveOutputFcn && exitFlag==-1
                stop = true;
                break;
            elseif exitFlag==0
                warning(message('stats:nlmefit:IterationLimitExceededLA'));
            elseif exitFlag==-3
                warning(message('stats:nlmefit:ObjectiveLimitExceededLA'));
            end
            
            MaxBetaDiff = max(abs(Beta-params_hat(1:q)));
            MaxBetabDiff = max([MaxBetaDiff MaxbDiff]);
            Beta = params_hat(1:q);
            theta = params_hat(q+1:q+dofCov);
            Delta = theta2Delta(theta,r,parType,patternCov);
            
            nlLADiff = abs(nlLAnew-nlLA);
            nlLA = nlLAnew;
            
            % Update weights
            switch errorModel
                case 'constant'
                    % weights = ones(N,1);
                case 'proportional'
                    fBetab = computeFAllGroupsFcn(Beta,b);
                    weights = 1./abs(fBetab);
                    weights(isinf(weights)) = 1;
                case 'combined'
                    fBetab = computeFAllGroupsFcn(Beta,b);
                    ab = abs(fminsearch(@(ab) error_ab(ab,Y,fBetab),errorModelParam));
                    weights = 1./abs(ab(1)+ab(2)*abs(fBetab));
                    weights(isinf(weights)) = 1;
                    errorModelParam = ab;
                case 'exponential'
                    % weights = ones(N,1);
            end
            
            if haveOutputFcn % call outer level iteration outputFcn
                [nlLA,sigma2_hat] = fcn([Beta;theta]);
                llike = Klike - nlLA; % constant part of the loglikelihood + llLA
                ofopt = struct('procedure','LAP','iteration',ALTniter+nIter-1,'inner',...
                    struct('procedure','none','state','none','iteration',NaN),...
                    'fval',llike,'Psi', Delta2Psi_sigma2(Delta,r)*sigma2_hat,...
                    'theta',theta,'mse',sigma2_hat,'caller','nlmefit');
                stop = callOutputFcns(outputFcn,Beta,ofopt,'iter');
                if stop, break; end
            end
            
            % Check output conditions for the iterative algorithm
            if nlLADiff <= tolF*(10.^tolFac) &&  MaxBetabDiff<= tolX*(10.^tolFac)
                if tolFac <= 0
                    break
                else
                    tolFac = tolFac - 0.2;
                    minMethodOpt.TolFun = tolF*(10.^tolFac);
                    minMethodOpt.TolX = tolX*(10.^tolFac);
                    GPNLSOpt.TolFun = tolF*(10.^tolFac);
                    GPNLSOpt.TolX = tolX*(10.^tolFac);
                end
            end
        end
    end
    
    
    if nIter>=maxIter
        warning(message('stats:nlmefit:IterationLimitExceededLaplace'));
    end
    
    Beta_hat = Beta;
    theta_hat = theta;
    Delta_hat = Delta;
    
    for i = 1:M
        b(:,i) = GNfitGPNLS(Delta,Y(grpMatch(:,i)),weights(grpMatch(:,i)),Beta,b(:,i),GPNLSOpt,i,r,computeFByGroupFcn,computeJacobianByGroupFcn);
    end
    b_hat = b;
    
    % Calculate remaining output values
    fcn = LAnegloglikelihoodFcnFcn(b,weights);
    [nlLA,sigma2_hat] = fcn([Beta;theta]);
    llike = -N/2 * (1+log(2*pi)) - nlLA;
    Psi_sigma2 = Delta2Psi_sigma2(Delta,r);
    Psi_hat = Psi_sigma2 * sigma2_hat;
    
    if verbose
        dispVal('Beta_hat''',Beta_hat')
        dispVal('theta_hat''',theta_hat')
        dispVal('logLike',llike)
        dispVal('sigma2_hat',sigma2_hat)
        dispVal('Psi_hat',Psi_hat,'%9.6f')
        dispVal('Delta_hat',Delta_hat,'%9.6f')
        dispVal('Deg.Freedom',N-numParam,'%d')
        dispVal('AIC',-2*llike+2*numParam)
        dispVal('BIC',-2*llike+log(M)*numParam)
        fprintf('%s',getString(message('stats:nlmefit:IterationsAndModelSeconds',nIter,getFCountFcn(),sprintf('%f',cputime-t))))
        %    Optional values that the user may want to display (when debugging)
        %    dispVal('b_hat''',b_hat','%6.3f')
        %    dispVal('log(sigma2_hat)',log(sigma2_hat))
        %    dispVal('log(diag(L))',log(diag(chol(Psi_hat./sigma2_hat)))')
        %    dispVal('triu(L,1)',nonzeros(triu(chol(Psi_hat./sigma2_hat),1))')
        %    disp(' ')
    end
    
elseif verbose>2
    disp(getString(message('stats:nlmefit:LaplacianApproximationOFF')))
    disp(' ')
end

if haveOutputFcn
    if stop
        warning(message('stats:nlmefit:AlgorithmTerminatedByOutputFcn'));
    else
        if laplacianMethod
            ofopt = struct('procedure','LAP','iteration',ALTniter+nIter-1,'inner',...
                struct('procedure','none','state','none','iteration',NaN),...
                'fval',llike,'Psi', Delta2Psi_sigma2(Delta,r)*sigma2_hat,...
                'theta',theta,'mse',sigma2_hat,'caller','nlmefit');
        else
            ofopt = struct('procedure','ALT','iteration',ALTniter-1,'inner',...
                struct('procedure','none','state','none','iteration',NaN),...
                'fval',llike,'Psi', Delta2Psi_sigma2(Delta,r)*sigma2_hat,...
                'theta',theta,'mse',sigma2_hat,'caller','nlmefit');
        end
        callOutputFcns(outputFcn,Beta,ofopt,'done');
    end
end

%% Prepare output structure for stats
%       logl   The maximized log-likelihood for the fitted model
%       mse    The estimated error variance for the fitted model
%                 Undocummented once rmse was introduced to be consistent
%                 with nlmefitsa, left in the output structure to be
%                 backwards compatible.
%       rmse   The root mean squared residual
%       errorparam  The estimated parameters of the error variance model
%       aic    The Akaike information criterion for the fitted model
%       bic    The Bayesian information criterion for the fitted model
%       sebeta The standard errors for BETA
%       dfe    The error degrees of freedom for the model

Yhat = computeFAllGroupsFcn(Beta,b);
res = Y - Yhat; % no weights here
stats.dfe = N-numParam;
stats.logl = llike;
stats.mse = sigma2_hat;
stats.rmse = sqrt(sum(abs(res).^2) / stats.dfe);
switch errorModel
    case 'combined'
        stats.errorparam = abs(fminsearch(@(ab) error_ab(ab,Y,Yhat),errorModelParam));
    otherwise % case {'constant', 'proportional', 'exponential'}
        % Note that this is the maximum likelihood estimate. The
        % restricted likelihood estimate is given by
        %     sqrt( (c_1.^2)/dfe )
        % See Pinheiro & Bates Eqs. 2.12 and 2.19, and discussion on p. 76
        % above Eq. 2.23.
        stats.errorparam = sqrt(sigma2_hat);
end
stats.aic = -2*llike+2*numParam;
stats.bic = -2*llike+log(M-q*restricted)*numParam;
R00inv = R00\eye(q);
stats.covb = (R00inv*R00inv')*sigma2_hat;
stats.sebeta = sqrt(diag(stats.covb))';

% Calculate residuals and return them in the stats structure.
regFcn = @(b) computeFAllGroupsFcn(Beta, b);
% Jb is the second output; use helper function to easily get this output.
jacFcn = @(f0, b) fevalSecondOutput(computeJacobianAllGroupsFcn, f0, Beta, b);

warnState = warning('off', 'stats:nlmeres:ErrorVarianceCloseToZero');
cleanupObj = onCleanup(@() warning(warnState));
residuals = nlmeres(regFcn, jacFcn, Y, grpMatch, Psi_hat, errorModel, stats.errorparam, b);

stats.ires  = residuals.ires;
stats.pres  = residuals.pres;
stats.iwres = residuals.iwres;
stats.pwres = residuals.pwres;
stats.cwres = residuals.cwres;
[~,id]    = lastwarn;
if strcmpi(id, 'stats:nlmeres:ErrorVarianceCloseToZero')
    warning(message('stats:nlmefit:ErrorVarianceCloseToZero'));
end
end % function nlmefit


function fcn = getLAnegloglikelihoodFcnFcn(Y,grpMatch,q,r,dofCov,M,N,parType,patternCov,model,computeJacobianAllGroupsFcn)
%% Burn in problem-specific variables for LAnegloglikelihood
fcn = @LAnegloglikelihoodFcnFcn;

    function fcn = LAnegloglikelihoodFcnFcn(b, weights)
    %% Burn in b, weights, and related parameters for LAnegloglikelihood
    
    diagInvLambda = weights.*weights;
    halfsumlogLambda = sum(log(weights));
    
    fcn = @LAnegloglikelihood;
    
        function [negloglike,sigma2] = LAnegloglikelihood(params)
        %% Compute the LAnegloglikelihood
        Beta = params(1:q);
        Delta = theta2Delta(params(q+1:q+dofCov),r,parType,patternCov);
        if ~all(isfinite(Delta))
            negloglike = 1/eps;
            return
        end
        DeltaTDelta = Delta'*Delta;
        
        % Compute the approximation to the Jacobian of g (G) (Eq.(7.18) in
        % Pinheiro&Bates, page 317), and in Pinheiro&Bates (page 331) for
        % non-constant error models. We calculate sum(log(|G|)) and sigma2 in the
        % same loop:
        fb = model(Beta,b);
        [~, dfdb] = computeJacobianAllGroupsFcn(fb, Beta, b);
        
        sigma2 = (sum(sum((Delta*b).^2)) + sum(bsxfun(@times,weights,(Y-fb)).^2))/N;
        detG = zeros(M,1);
        for i = 1:M
            G = dfdb(grpMatch(:,i),:)' * diag(diagInvLambda(grpMatch(:,i))) * dfdb(grpMatch(:,i),:) + DeltaTDelta;
            detG(i) = det(G);
        end
        
        % Varying part of the negative loglikelihood Laplacian approximation
        % profiled on sigma2 given by  Eq.(7.19) in Pinheiro&Bates (page 319) and
        % in Pinheiro&Bates (page 331) for non-constant error models.
        negloglike =  - M*log(det(Delta)) + sum(log(detG))/2 + N*log(sigma2)/2 - halfsumlogLambda;
        end % function LAnegloglikelihood
    
    end % function LAnegloglikelihoodFcnFcn

end % function getLAnegloglikelihoodFcnFcn

function fcn = getLMEnegloglikelihoodFcnFcn(r,q,M,N,grpMatch,restricted,parType,patternCov)
%% Burn in problem-specific variables for LMEnegloglikelihood
R00c0 = zeros(N,q+1);
detR11 = zeros(M,1);

fcn= @LMEnegloglikelihoodFcnFcn;

    function fcn = LMEnegloglikelihoodFcnFcn(Z,Xy)
    %% Burn in Z and Xy for LMEnegloglikelihood
    fcn = @LMEnegloglikelihood;
    
        function [negloglike,c_1,R00] = LMEnegloglikelihood(theta)
        %% Compute the LMEnegloglikelihood
        Delta = theta2Delta(theta,r,parType,patternCov);
        if ~all(isfinite(Delta))
            negloglike = 1/eps;
            return
        end
        
        for i = 1:M
            % QR decomposition for LME estimation (Pinheiro&Bates, page 68)
            [Q_i,Rx1_i] = qr([Z(grpMatch(:,i),:);Delta]);
            if r == 1
                detR11(i) = abs(Rx1_i(1));
            else
                detR11(i) = prod(diag(Rx1_i));
            end
            % Building up Eq.(2.17) Pinheiro&Bates, page 70
            Rx0cx_i = Q_i'*[Xy(grpMatch(:,i),:);zeros(r,q+1)];
            R00c0(grpMatch(:,i),:) = Rx0cx_i(r+1:end,:);
        end
        R00cx = qr(R00c0,0);
        c_1 = abs(R00cx(q+1,q+1));
        if restricted  % Varying part of Eq. (2.23)
            negloglike = (N-q)*log(c_1) + log(abs(prod(diag(R00cx(1:q,1:q))))) - sum(log(abs(det(Delta)./detR11)));
        else   % Varying part of Eq. (2.21)
            negloglike = N*log(c_1) - sum(log(abs(det(Delta)./detR11)));
        end
        if nargout>2
            R00 = triu(R00cx(1:q,1:q));
        end
        end % function LMEnegloglikelihood
    
    end % function LMEnegloglikelihoodFcnFcn

end % function getLMEnegloglikelihoodFcnFcn

function b = GNfitGPNLS(Delta,Y,weights,Beta,b,options,i,r,computeFByGroupFcn,computeJacobianByGroupFcn)
%% Gauss-Newton algorithm for nonlinear regression of the random effects for the i-th group.

ptol = options.TolX;
rtol = options.TolFun;
funValCheck = strcmp(options.FunValCheck, 'on');
maxIter = options.MaxIter;

hD = diag(Delta) > min(diag(Delta)).*1e8;

res = (Y - computeFByGroupFcn(Beta,b,i)).*weights;
sse = res'*res + sum((Delta*b).^2);

iter = 0;
while iter < maxIter
    % Taylor series approximation about the current estimates:
    %  1) compute df/db | b_current with finite differences
    f0 = computeFByGroupFcn(Beta,b,i);
    [~, Zw] = computeJacobianByGroupFcn(f0, Beta, b, i);
    %  2) Build the Delta-independent part of the augmented matrix
    ww = Y-f0+Zw*b;
    
    % The derivative matrices and working vector are weighted by the error
    % model:
    ww = bsxfun(@times,weights,ww);
    Zw = bsxfun(@times,weights,Zw);
    
    if funValCheck, checkFunVals(ww(:)); end
    
    %  3) Find the least squares estimates
    R = qr([Zw ww;Delta zeros(r,1)],0);
    R11 = triu(R(1:r,1:r));
    c1 = R(1:r,r+1);
    %c0 = R(r+1,r+1);
    %b_new = R11\c1;
    
    if any(hD)
        b_new = zeros(r,1);
        b_new(~hD) = linsolve(R11(~hD,~hD),c1(~hD),struct('UT',true));
    else
        b_new = linsolve(R11,c1,struct('UT',true));
    end
    
    % Evaluate the fitted values at the new coefficients and
    % compute the residuals and the SSE.
    res = (Y - computeFByGroupFcn(Beta,b_new,i)).*weights;
    sse_new = res'*res + sum((Delta*b_new).^2);
    if funValCheck && ~isfinite(sse_new), checkFunVals(sse_new); end
    
    step = b_new - b;
    
    while sse_new > sse
        if max(abs(step))<eps % check for stall
            warning(message('stats:nlmefit:UnableToDecreaseSSEinPNLS'));
            sse_new = sse;
            break
        end
        step = step/2;
        b_new = b +step;
        res = (Y - computeFByGroupFcn(Beta,b_new,i)).*weights;
        sse_new = res'*res + sum((Delta*b_new).^2);
        if funValCheck && ~isfinite(sse_new), checkFunVals(sse_new); end
    end
    
    sseDiff = abs(sse_new - sse);
    maxParDiff = max(abs(step));
    
    b = b_new;
    sse = sse_new;
    
    % Check output conditions
    if sseDiff <= rtol && maxParDiff <= ptol
        break
    end
    iter = iter+1;
end
if iter>=maxIter
    warning(message('stats:nlmefit:IterationLimitExceededPNLS'));
end
end % function GNfitGPNLS

function  [Beta,b,stop] = GNfitPNLS(Delta,Y,Beta,b,weights,options,grp,grpMatch,M,N,r,q,computeFAllGroupsFcn,computeJacobianAllGroupsFcn)
%% Gauss-Newton algorithm for coupled nonlinear regression of the fixed and random effects.

stop = false;
outputFcn = options.OutputFcn;
ptol = options.TolX;
rtol = options.TolFun;
funValCheck = strcmp(options.FunValCheck, 'on');
maxIter = options.MaxIter;

hD = diag(Delta) > min(diag(Delta)).*1e8;

res = (Y - computeFAllGroupsFcn(Beta,b)).*weights;
sse = res'*res + sum(sum((Delta*b).^2));

R11 = cell(M,1);
R00c0 = zeros(N,q+1);
R10c1 = zeros(r*M,q+1);
b_new = zeros(r,M);

iter = 0;
if ~isempty(outputFcn)
    stop = outputFcn(Beta,iter,'init');
end
if ~stop
    while iter < maxIter
        % LINEARIZATION step, (Eq.(7.11) Pinheiro&Bates, page 313) (all groups at
        % once)
        f0 = computeFAllGroupsFcn(Beta,b);
        [Xw, Zw] = computeJacobianAllGroupsFcn(f0, Beta, b);
        
        % Build the Delta-independent part of the augmented matrix (ZwXwww)
        % that is required to solve LME (Pinheiro&Bates, page 69) (all groups at
        % once)
        ww = Y - f0 + Xw*Beta + sum(Zw.*b(:,grp)',2);
        
        if funValCheck, checkFunVals(ww(:)); end
        Xwww = [Xw ww];
        
        % The derivative matrices and working vector are weighted by the error
        % model:
        Xwww = bsxfun(@times,weights,Xwww);
        Zw = bsxfun(@times,weights,Zw);
        
        for i = 1:M
            % QR decomposition for Linear Mixed Effects estimation (Pinheiro&Bates, page 69)
            [Q,R]=qr([Zw(grpMatch(:,i),:);Delta]);
            R11{i} = R(1:r,1:r);
            % Building up Eq.(2.17) Pinheiro&Bates, page 70
            Rx0cx = Q'*[Xwww(grpMatch(:,i),:);zeros(r,q+1)];
            R00c0(grpMatch(:,i),:) = Rx0cx(r+1:end,:);
            R10c1(r*(i-1)+(1:r),:) = Rx0cx(1:r,:); % used later in the next loop
        end
        % Find the residual vector for the penalized least-squares fit:
        R00c0 = qr(R00c0,0);
        %Beta_new = triu(R00c0(1:q,1:q)) \ R00c0(1:q,q+1);
        Beta_new = linsolve(triu(R00c0(1:q,1:q)),R00c0(1:q,q+1),struct('UT',true));
        if any(hD)
            for i = 1:M  % Compute BLUPs (Eq.(2.22) Pinheiro&Bates, page 71)
                c1mR10Beta = (R10c1(r*(i-1)+(1:r),q+1)-R10c1(r*(i-1)+(1:r),1:q)*Beta_new);
                b_new(~hD,i) = linsolve(R11{i}(~hD,~hD),c1mR10Beta(~hD),struct('UT',true));
            end
        else
            for i = 1:M  % Compute BLUPs (Eq.(2.22) Pinheiro&Bates, page 71)
                %b_new(:,i) = R11{i}\(R10c1(r*(i-1)+(1:r),q+1)-R10c1(r*(i-1)+(1:r),1:q)*Beta_new);
                b_new(:,i) = linsolve(R11{i},(R10c1(r*(i-1)+(1:r),q+1)-R10c1(r*(i-1)+(1:r),1:q)*Beta_new),struct('UT',true));
            end
        end
        
        % Evaluate the fitted values at the new coefficients and
        % compute the residuals and the SSE.
        res = (Y - computeFAllGroupsFcn(Beta_new,b_new)).*weights;
        sse_new = res'*res + sum(sum((Delta*b_new).^2));
        if funValCheck && ~isfinite(sse_new), checkFunVals(sse_new); end
        
        step = [Beta_new;b_new(:)] - [Beta;b(:)];
        
        while sse_new > sse
            if max(abs(step))<eps % check for stall
                warning(message('stats:nlmefit:UnableToDecreaseSSEinPNLS'));
                sse_new = sse;
                break
            end
            step = step/2;
            Beta_new = Beta + step(1:q);
            b_new(:) = b(:) + step(q+1:end);
            res = (Y - computeFAllGroupsFcn(Beta_new,b_new)).*weights;
            sse_new = res'*res + sum(sum((Delta*b_new).^2));
            if funValCheck && ~isfinite(sse_new), checkFunVals(sse_new); end
        end
        
        sseDiff = abs(sse_new - sse);
        maxParDiff = max(abs(step));
        
        Beta = Beta_new;
        b = b_new;
        sse = sse_new;
        
        if ~isempty(outputFcn)
            stop = outputFcn(Beta,iter,'iter');
            if stop
                break
            end
        end
        % Check output conditions
        if sseDiff <= rtol && maxParDiff <= ptol
            if ~isempty(outputFcn)
                stop = outputFcn(Beta,iter,'done');
            end
            break
        end
        iter = iter+1;
    end
end
if iter>=maxIter
    warning(message('stats:nlmefit:IterationLimitExceededPNLS'));
end
end % function GNfitPNLS

function [computeFAllGroupsFcn, computeFByGroupFcn, getFCountFcn, computeJacobianAllGroupsFcn, computeJacobianByGroupFcn] ...
    = getFcns(X,V,A,B,model,grp,grpMatch,modelVectorization,paramTransform,fdiffstep,errorModel)
%% Generate function handles to call the model and estimate the Jacobian

% Get the phi functions and the functions to convert JPhi to JBeta & JB.
% the chain rule to the Jacobian with respect to phi.
[phiAllGroupsFcn, phiByGroupFcn, JBetaFcn, JBFcn, JBetaByGroupFcn, JBByGroupFcn] = getPhiFcns(X,V,A,B,grp,grpMatch);

model = applyModelTransformations(model, errorModel, paramTransform, V);

% Construct function evaluation counter. count is shared with the various
% computeF functions.
count = 0;
    function c = getCount()
    c = count;
    end
getFCountFcn = @getCount;

% Store sizes of the design matrices.
Am = size(A,3);
Bm = size(B,3);
M = size(V,1);
N = size(X,1);

%% Helper functions for vectorization of the model.
% The name indicates the type of call inside the for loop:
% * vectorizeSinglePhiSingleX calls the model with a single phi vector and
% a single observation vector;
% * vectorizeSinglePhiSingleXByGroup calls the model (for a single group)
% with a single phi vector and a single observation vector;
% * vectorizeGroupPhiGroupX calls the model with a group's phi matrix and a
% group's observation matrix; and
% * vectorizeSinglePhiGroupX calls the model with a single phi vector and a
% group's observation matrix.

    function fval = vectorizeSinglePhiSingleX(phi,X,V)
    fval = zeros(N,1);
    for i=1:N
        fval(i) = model(phi(i,:),X(i,:),V(grp(i),:));
    end
    end % function vectorizeSinglePhiSingleX

    function fval = vectorizeSinglePhiSingleXByGroup(phi,X,V)
    nX = size(X,1); % number of observations in this group
    fval = zeros(nX,1);
    for i=1:nX
        fval(i) = model(phi(i,:),X(i,:),V);
    end
    end % function vectorizeSinglePhiSingleXByGroup

    function fval = vectorizeGroupPhiGroupX(phi,X,V)
    fval = zeros(N,1);
    for i=1:M
        k = grpMatch(:,i);
        fval(k) = model(phi(k,:),X(k,:),V(i,:));
    end
    end % function vectorizeGroupPhiGroupX

    function fval = vectorizeSinglePhiGroupX(phi,X,V)
    fval = zeros(N,1);
    for i=1:M
        k = grpMatch(:,i);
        fval(k) = model(phi(i,:),X(k,:),V(i,:));
    end
    end % function vectorizeSinglePhiGroupX

%% Vectorize the model function if necessary.
% Also set deltaCountAllGroups, the number of function evaluations per call
% to the vectorized function handle; set deltaCountByGroupFcn, the number
% of function evaluations per call to the group-based function handle. Note
% that this increment depends on the number of observations in a group, so
% this variable is a function handle that is a function of the logical
% vector comparing the group number to the variable GRP (the vector of
% group identifiers for each observation).

% If modelVectorization is singlephi or singlegroup, we call the model in
% one of three ways, depending on the design matrices and vectorization:
% 1) a single PHI vector  and a single X vector  (looping over observations)
% 2) a single PHI vector  and multiple X vectors (looping over groups)
% 3) multiple PHI vectors and multiple X vectors (looping over groups)

% We may also need to "vectorize" V, if the function is
% observation-specific and fully vectorized. By default we assume the
% function accepts group-specific covariates V
VVectorized = V;
if strcmp(modelVectorization, 'full')
    % Completely vectorized; one evaluation per call.
    deltaCountAllGroups = 1;
    modelVectorizedAllGroups = model;
    deltaCountByGroupFcn = @(i) 1;
    modelVectorizedByGroup = model;
    
    % If a design matrix is observation-specific, then the covariates V
    % need to be expanded to be observation-specific, as well.
    if Am == N || Bm == N
        VVectorized = V(grp,:);
        % Observation-specific indexing in statjacobian
        jacobianRowIdx = 1:N;
        jacobianByGroupRowIdxFcn = @(f0) 1:numel(f0);
    else
        % Group-specific indexing in statjacobian
        jacobianRowIdx = grp;
        jacobianByGroupRowIdxFcn = @(f0) [];
    end
    
elseif Am == N || Bm == N
    % One or both design matrices are observation-specific.
    % Observation-specific indexing in statjacobian
    jacobianRowIdx = 1:N;
    jacobianByGroupRowIdxFcn = @(f0) 1:numel(f0);
    if strcmp(modelVectorization, 'singlephi')
        % Call the full model with a single observation-specific PHI and a
        % single observation X, looping over observations.
        deltaCountAllGroups = N;
        modelVectorizedAllGroups = @vectorizeSinglePhiSingleX;
        
        % The model must also be vectorized for group-specific evaluations.
        grpMatchCount = sum(grpMatch, 1);
        deltaCountByGroupFcn = @(i) grpMatchCount(i);
        modelVectorizedByGroup = @vectorizeSinglePhiSingleXByGroup;
    else % modelVectorization = singlegroup
        % Call the full model with a group of PHI and X values, looping
        % over groups.
        deltaCountAllGroups = M;
        modelVectorizedAllGroups = @vectorizeGroupPhiGroupX;
        
        % The model is already vectorized for group-specific evaluations.
        deltaCountByGroupFcn = @(i) 1;
        modelVectorizedByGroup = model;
    end
else
    % We have a group-specific PHI, so we can call the model with multiple
    % X values, looping over groups.

    % Group-specific indexing in statjacobian
    jacobianRowIdx = grp;
    jacobianByGroupRowIdxFcn = @(f0) [];
    deltaCountAllGroups = M;
    modelVectorizedAllGroups = @vectorizeSinglePhiGroupX;
    deltaCountByGroupFcn = @(i) 1;
    modelVectorizedByGroup = model;
end

    function f = computeFAllAndCount(phi)
    count = count + deltaCountAllGroups;
    f = modelVectorizedAllGroups(phi, X, VVectorized);
    end

    function f = computeFByGroupAndCount(phi,i)
    count = count + deltaCountByGroupFcn(i);
    f = modelVectorizedByGroup(phi, X(grpMatch(:,i),:), V(i,:));
    end

%% Create computeF function handles
computeFAllGroupsFcn = @(Beta,b) computeFAllAndCount(phiAllGroupsFcn(Beta,b));
computeFByGroupFcn = @(Beta,b,i) computeFByGroupAndCount(phiByGroupFcn(Beta,b,i),i);

%% Create Jacobian function handles
    function [JBeta, Jb] = jacobian(f0, Beta, b)
    % Compute the jacobian of the function with respect to Beta and/or b
    % Uses the chain rule to estimate JBeta nd Jb from Jphi.
    phi = phiAllGroupsFcn(Beta,b);
    Jphi = statjacobian(@computeFAllAndCount, phi, fdiffstep, f0, jacobianRowIdx);
    JBeta = JBetaFcn(Jphi);
    if nargout > 1
        Jb = JBFcn(Jphi);
    end
    end % function jacobian

    function [JBeta, Jb] = jacobianByGroup(f0, Beta, b, i)
    % Compute the jacobian of the function with respect to Beta and/or b
    % Uses the chain rule to estimate JBeta and Jb from Jphi.
    phi = phiByGroupFcn(Beta,b,i);
    fFcn = @(phi) computeFByGroupAndCount(phi, i);
    Jphi = statjacobian(fFcn, phi, fdiffstep, f0, jacobianByGroupRowIdxFcn(f0));
    JBeta = JBetaByGroupFcn(Jphi, i);
    if nargout > 1
        Jb = JBByGroupFcn(Jphi, i);
    end
    end % function jacobian

computeJacobianAllGroupsFcn = @jacobian;
computeJacobianByGroupFcn = @jacobianByGroup;

end % function getFcns

function model = applyModelTransformations(model, errorModel, paramTransform, V)
%% Add parameter transformations and optional covariates to the model function handle.
if strcmp(errorModel, 'exponential')
    % Exponential error model. Linearize the model as
    %   y = f*exp(a*e), or log(y) = log(f) + a*e
    model = @(phi,X,varargin) log(max(model(phi,X,varargin{:}),realmin));
end

if isempty(V);
    % ignore V input argument
    if any(paramTransform)
        model = @(phi,X,V) model(transphi(phi,paramTransform),X);
    else
        model = @(phi,X,V) model(phi,X);
    end
else
    if any(paramTransform)
        model = @(phi,X,V) model(transphi(phi,paramTransform),X,V);
    end
end
end % function applyModelTransformations

function z = svtimes(x,y)
% Do mtimes along each slice in x and column vector y.
n = size(x,3);
z = zeros(size(x,1), n);
for i = 1:n
    z(:,i) = x(:,:,i)*y;
end
end % function svtimes

function z = sctimes(x,y)
% Do mtimes along each slice in x and each column in matrix y.
n = size(x,3);
z = zeros(size(x,1), n);
for i = 1:n
    z(:,i) = x(:,:,i)*y(:,i);
end
end % function sctimes

function [phiAllGroupsFcn, phiByGroupFcn, ...
    JBetaFcn, JBFcn, JBetaByGroupFcn, JBByGroupFcn] = ...
    getPhiFcns(X,V,A,B,grp,grpMatch)
%% Create phi function handles
% * phiAllGroupsFcn is a function handle taking input Beta and
% returning a vector Y of model predictions for all groups.
% * phiByGroupFcn is a function handle taking inputs Beta, B, and group
% index i, returning a vector Y of model predictinos for group i.
% * JBetaFcn is a function handle taking input JPhi and returning JBeta.
% * JBFcn is a function handle taking input JPhi returning Jb.
M = size(V,1);
N = size(X,1);
Am = size(A,3);
Bm = size(B,3);

% Construct functions to multiply A and Beta
if Am == 1 % Constant fixed effect
    AMultFcn = @mtimes;
    AMultByGroupFcn = @(A,Beta,j) A*Beta;
elseif Am == M % Group-specific fixed effect
    AMultFcn = @svtimes;
    AMultByGroupFcn = @(A,Beta,j) A(:,:,j)*Beta;
elseif Am == N % Observation-specific fixed effect
    AMultFcn = @svtimes;
    AMultByGroupFcn = @(A,Beta,j) svtimes(A(:,:,grpMatch(:,j)),Beta);
end

% Construct functions to multiply B and b
if Bm == 1 % Constant random effect
    BMultFcn = @mtimes;
    BMultByGroupFcn = @(B,b,j) B*b;
elseif Bm == M % Group-specific random effect
    BMultFcn = @sctimes;
    BMultByGroupFcn = @(B,b,j) B(:,:,j)*b;
elseif Bm == N % Observation-specific random effect
    BMultFcn = @(B,b)sctimes(B,b(:,grp));
    BMultByGroupFcn = @(B,b,j) svtimes(B(:,:,grpMatch(:,j)),b);
end

% Create a helper function that expands columns from group-specific to
% observation-specific.
expandColumnsFromGroupToObs = @(X) X(:,grp);

% Construct functions to add A*Beta and B*b
if Am == M && Bm == N
    % Need to expand group-specific A*Beta to add to observation-specific B*b
    phiAllGroupsFcn = @(Beta,b) bsxfun(@plus, ...
        expandColumnsFromGroupToObs(AMultFcn(A,Beta)), BMultFcn(B,b))';
elseif Am == N && Bm ~= N
    % (B = 1 or B = M results in group-specific B*b.)
    % Need to expand group-specific B*b to add to observation-specific A*Beta
    phiAllGroupsFcn = @(Beta,b) bsxfun(@plus, AMultFcn(A,Beta), ...
        expandColumnsFromGroupToObs(BMultFcn(B,b)))';
else
    % We can use bsxfun to add the A*Beta and B*b
    phiAllGroupsFcn = @(Beta,b) bsxfun(@plus, AMultFcn(A,Beta), BMultFcn(B,b))';
end
phiByGroupFcn = @(Beta,b,j) bsxfun(@plus, AMultByGroupFcn(A,Beta,j), ...
    BMultByGroupFcn(B,b,j))';

% Construct the JBeta and JB function handles.
% Via the chain rule, JPhi can be converted to JBeta and JB. Since the
% design matrices are linear, this conversion is a form of multiplication.
[JBetaFcn, JBetaByGroupFcn] = getJChainRuleFcn(A, M, N, grpMatch);
[JBFcn, JBByGroupFcn] = getJChainRuleFcn(B, M, N, grpMatch);
end % function getPhiFcns

function [JChainRuleFcn, JChainRuleByGroupFcn] = getJChainRuleFcn(X, M, N, grpMatch)
% Use the chain rule to convert JPhi to JBeta or JB, using the design
% matrix X (where X is A for JBeta and B for JB).

% JChainRuleFcn(JPhi) operates on the matrix JPhi, where JPhi(i,j) is the
% Jacobian of the i'th observation (across all groups) with respect to
% Phi(j).

% JChainRuleByGroupFcn(JPhi_g, g) operates on the matrix JPhi_g, where
% JPhi_g(i,j) is the Jacobian of the i'th observation of group g with
% respect to Phi(j).

sizeX2 = size(X, 2);
sizeX3 = size(X, 3);

if sizeX3 == 1
    JChainRuleFcn = @(JPhi) JPhi*X;
    JChainRuleByGroupFcn = @(JPhi_g, g) JPhi_g*X;
elseif sizeX3 == N
    % Observation-specific design matrix; multiply each observation
    % separately.
    JChainRuleFcn = @obsTimes;
    JChainRuleByGroupFcn = @obsTimesByGroup;
else % sizeX3 == M
    % Group-specific design matrix; multiply each group separately.
    JChainRuleFcn = @grpTimes;
    JChainRuleByGroupFcn = @grpTimesByGroup;
end

    function JX = obsTimes(JPhi)
    JX = zeros(N, sizeX2);
    for i = 1:N
        JX(i,:) = JPhi(i,:)*X(:,:,i);
    end
    end % function obsTimes

    function JX_g = obsTimesByGroup(JPhi_g, g)
    X_g = X(:,:,grpMatch(:,g));
    nObs = size(X_g, 3);
    JX_g = zeros(nObs, sizeX2);
    for i = 1:nObs
        JX_g = JPhi_g(i,:)*X_g(:,:,i);
    end
    end % function obsTimesByGroup
    
    function JX = grpTimes(JPhi)
    JX = zeros(N, sizeX2);
    for i = 1:M
        JX(grpMatch(:,i),:) = JPhi(grpMatch(:,i),:)*X(:,:,i);
    end
    end % function grpTimes
    
    function JX_g = grpTimesByGroup(JPhi_g, g)
    JX_g = JPhi_g*X(:,:,g);
    end % function grpTimesByGroup
end % function getJChainRuleFcn

%% Additional helper functions
function Psi_sigma2 = Delta2Psi_sigma2(Delta,r)
% Construct Psi_sigma2 given Delta, r is the number of random effects
ws = warning('off','MATLAB:nearlySingularMatrix');
Deltainv = linsolve(Delta,eye(r));
warning(ws)
Psi_sigma2 = Deltainv*Deltainv';
end % function Delta2Psi_sigma2

function theta = Delta2theta(Delta,r,parType,pat)
% Construct theta given Delta, sigma2 and a predefined parameterization
% method given by parType, r is the number of random effects
switch parType
    case 'logm' % Parameterizing logm of (Delta'*Delta)
        log_DeltaTDelta = logm(Delta'*Delta);
        theta = log_DeltaTDelta(pat);
    case 'logmcov'  % Parameterizing logm of the Scaled Covariance (Pinheiro&Bates, page 78)
        Deltainv = linsolve(Delta,eye(r));
        Psi_sigma2 = logm(Deltainv*Deltainv');
        theta = Psi_sigma2(pat);
    case 'chol' % Parameterizing the Precision Factors directly, i.e.
        % parameterizing the Cholesky Decomposition of (Delta'*Delta)
        theta = Delta(pat);
    case 'cholcov' % Parameterizing the Cholesky Decomposition of the Scaled Covariance
        L = linsolve(Delta,eye(r));
        theta = L(pat);
end
end % function Delta2theta

function Delta = theta2Delta(theta,r,parType,pat)
% Construct Delta given theta and a predefined parameterization method
% given by parType, r is the number of random effects
switch parType
    case 'chol' % Parameterizing the Precision Factors directly, i.e.
        % parameterizing the Cholesky Decomposition of (Delta'*Delta)
        Delta = zeros(r);
        Delta(pat) = theta;
    case 'logm' % Parameterizing logm of (Delta'*Delta)
        log_DeltaTDelta = zeros(r);
        log_DeltaTDelta(pat) = theta;
        log_DeltaTDelta = log_DeltaTDelta+triu(log_DeltaTDelta,1)';
        DeltaTDelta = expm(log_DeltaTDelta);
        try
            Delta = chol(DeltaTDelta);
        catch ME
            if isequal(ME.identifier,'MATLAB:posdef')
                Delta = NaN;
            else
                rethrow(ME)
            end
        end
    case 'logmcov' % Parameterizing logm of the Scaled Covariance (Pinheiro&Bates, page 78)
        log_Psi_sigma2 = zeros(r);
        log_Psi_sigma2(pat) = theta;
        log_Psi_sigma2 = log_Psi_sigma2+triu(log_Psi_sigma2,1)';
        Psi_sigma2 = expm(log_Psi_sigma2);
        if rank(Psi_sigma2)<r
            Delta = NaN;
            return
        end
        DeltaTDelta = linsolve(Psi_sigma2,eye(r),struct('SYM',true));
        try
            Delta = chol(DeltaTDelta);
        catch ME
            if isequal(ME.identifier,'MATLAB:posdef')
                Delta = NaN;
            else
                rethrow(ME)
            end
        end
    case 'cholcov' % Parameterizing the Cholesky Decomposition of the Scaled Covariance
        L = zeros(r);
        L(pat) = theta;
        if rank(L)<r
            Delta = NaN;
        else
            Delta = linsolve(L,eye(r));
        end
end
if ~all(diag(Delta)>=0)
    Delta = NaN;
end
end % function theta2Delta

function  beta = LMfit(y,beta,options,model,jacobian)
% Levenberg-Marquardt algorithm for nonlinear regression

lambda = 0.01; % 'lambda' is the initial weight for LM algorithm
betatol = options.TolX;
rtol = options.TolFun;
funValCheck = strcmp(options.FunValCheck, 'on');
maxiter = options.MaxIter;

% Set the iteration step
sqrteps = sqrt(eps(class(beta)));
p = numel(beta);

yfit = model(beta);
r = y - yfit;
sse = r'*r;

step = 0; % in case Jacobian is always 0, step is never defined
zerosp = zeros(p,1,class(r));
iter = 0;
breakOut = false;

while iter < maxiter
    iter = iter + 1;
    betaold = beta;
    sseold = sse;
    
    % Compute a finite difference approximation to the Jacobian
    J = jacobian(yfit, beta);
    
    % Levenberg-Marquardt step: inv(J'*J+lambda*D)*J'*r
    diagJtJ = sum(abs(J).^2, 1);
    if funValCheck && ~all(isfinite(diagJtJ)), checkFunVals(J(:)); end
    Jplus = [J; diag(sqrt(lambda*diagJtJ))];
    rplus = [r; zerosp];
    h = any(Jplus,1);
    if any(h)
        step = Jplus(:,h) \ rplus;
        beta(h) = beta(h) + step;
    end
    
    % Evaluate the fitted values at the new coefficients and
    % compute the residuals and the SSE.
    yfit = model(beta);
    r = y - yfit;
    sse = r'*r;
    if funValCheck && ~isfinite(sse), checkFunVals(r); end
    % If the LM step decreased the SSE, decrease lambda to downweight the
    % steepest descent direction.  Prevent underflowing to zero after many
    % successful steps; smaller than eps is effectively zero anyway.
    if sse < sseold
        lambda = max(0.1*lambda,eps);
        
        % If the LM step increased the SSE, repeatedly increase lambda to
        % upweight the steepest descent direction and decrease the step size
        % until we get a step that does decrease SSE.
    else
        while sse >= sseold
            lambda = 10*lambda;
            if lambda > 1e16
                breakOut = true;
                break
            end
            Jplus = [J; diag(sqrt(lambda*sum(J.^2,1)))];
            h = any(Jplus,1);
            if any(h)
                step = Jplus(:,h) \ rplus;
                beta(h) = betaold(h) + step;
            end
            yfit = model(beta);
            r = y - yfit;
            sse = r'*r;
            if funValCheck && ~isfinite(sse), checkFunVals(r); end
        end
    end
    % Check step size and change in SSE for convergence.
    if (norm(step) < betatol*(sqrteps+norm(beta))) &&...
            (abs(sse-sseold) <= rtol*sse)
        break
    elseif breakOut
        warning(message('stats:nlmefit:UnableToDecreaseSSEinLMalg'));
        break
    end
end
if (iter >= maxiter)
    warning(message('stats:nlmefit:IterationLimitExceededLMalg'));
end
% If the Jacobian is ill-conditioned, then two parameters are probably
% aliased and the estimates will be highly correlated.  Prediction at new x
% values not in the same column space is dubious. It may also be that the
% Jacobian has one or more columns of zeros, meaning model is constant with
% respect to one or more parameters.  This may be because those parameters
% are not even in the expression in the model function, or they are
% multiplied by another param that is estimated at exactly zero (or
% something similar), or because some part of the model function is
% underflowing, making it a constant zero. In the context of NLMEFIT it is
% preferable to error (instead of warning as in NLINFIT) since it is very
% likely that the NLME model will have the same issue.
[~,R] = qr(J,0);
if condest(R) > 1/(eps(class(beta)))^(1/2)
    if any(all(abs(J)<sqrt(eps(norm(J,1))),1),2) % one or more columns of zeros
        error(message('stats:nlmefit:ModelConstantWRTParamLMalg'))
    else  % no columns of zeros
        error(message('stats:nlmefit:IllConditionedJacobianLMalg'))
    end
end
end % function LMfit

function checkFunVals(v)
% Helper function to check if the function has the finite output
if any(~isfinite(v))
    error(message('stats:nlmefit:checkFunVals'));
end
end % function checkFunVals

function dispVal(str,x,f)
% Helper function to display values in the command window
if nargin<3, f = '%f'; end
if numel(x)==1
    disp([blanks(5) str blanks(12-numel(str)) '= ' sprintf(f,x)])
elseif size(x,1)==1
    disp([blanks(5) str blanks(12-numel(str)) '= [' sprintf([' ' f],x) ' ]'])
else
    disp([blanks(5) str blanks(12-numel(str)) '= [' sprintf([' ' f],x(1,:)) ])
    for i = 2:size(x,1)-1
        disp([blanks(20) sprintf([' ' f],x(i,:))])
    end
    disp([blanks(20) sprintf([' ' f],x(end,:)) ' ]'])
end
end % function dispVal

function [p,r,A,B,patternCov,thetaLen,RefineBeta0,Approximationtype,...
    modelVectorization,CovParameterization,Options,Optimfun,OutputFcn,...
    ParamTransform,ErrorModel,ErrorParameters] = parseInVarargin(q,M,N,varargin)
% Helper function for parsing and validate VARARGIN

% Set defaults
RESELECTgiven = false;
REDESIGNgiven = false;
REGROUPDESIGNgiven = false;
REOBSDESIGNgiven = false;
FESELECTgiven = false;
FEDESIGNgiven = false;
FEGROUPDESIGNgiven = false;
FEOBSDESIGNgiven = false;
COVPATTERN = NaN;
RE_p = NaN;
RE_r = NaN;
FE_p = NaN;
PSI_r = NaN;

RefineBeta0 = true;
Approximationtype = 'LME';
modelVectorization = 'SinglePhi';
CovParameterization = 'logm';
Optimfun = 'fminsearch';
OutputFcn = [];
ParamTransform = [];
ErrorModel = 'constant';
ErrorParameters = [];

Options = statset('nlmefit');
% Options.TolX = 1e-4;
% Options.DerivStep = eps^(1/3);
% Options.TolFun = 1e-4;
% Options.MaxIter = 200;
% Options.Display = 'off';
% Options.FunValCheck = 'on';
% Options.Jacobian = 'off';
% Options.OutputFcn = [];

% Process PVP
numvaragin = numel(varargin);
if numvaragin > 0
    if rem(numvaragin,2)
        error(message('stats:nlmefit:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'covpattern','reparamsselect','redesign','regroupdesign',...
        'reobsdesign','feparamsselect','fedesign','fegroupdesign',...
        'feobsdesign','refinebeta0','approximationtype','vectorization',...
        'covparameterization','options','optimfun','reconstdesign',...
        'feconstdesign','outputfcn','paramtransform','errormodel','errorparameters'};
    for j=1:2:numvaragin
        pname = varargin{j};
        if ~ischar(pname)
            error(message('stats:nlmefit:IllegalParamName'));
        end
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs, length(pname)));
        if isempty(k)
            error(message('stats:nlmefit:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('stats:nlmefit:AmbiguousParameterName', pname));
        end
        param = okargs{k};
        switch param
            case 'covpattern'
                if isvector(pval)
                    pval = grp2idx(pval);
                    COVPATTERN = false(numel(pval));
                    for i = 1:max(pval)
                        COVPATTERN(pval==i,pval==i) = true(sum(pval==i));
                    end
                elseif isnumeric(pval) && ndims(pval)==2 && ~diff(size(pval))
                    COVPATTERN = pval ~= 0;
                elseif ndims(pval)==2 && ~diff(size(pval))
                    COVPATTERN = pval;
                else
                    error(message('stats:nlmefit:InvalidPAT'))
                end
                PSI_r = size(COVPATTERN,1);
                if ~isequal(COVPATTERN,(COVPATTERN|COVPATTERN')|eye(PSI_r))
                    warning(message('stats:nlmefit:SymmetricPAT'))
                    COVPATTERN = (COVPATTERN|COVPATTERN')|eye(PSI_r);
                end
                if ~isequal(COVPATTERN,(COVPATTERN).^PSI_r)
                    warning(message('stats:nlmefit:MissingCovInPAT'))
                    COVPATTERN = COVPATTERN.^PSI_r;
                end
            case 'reparamsselect'
                if isvector(pval) && islogical(pval)
                    RE_p = numel(pval);
                    RE_r = sum(pval);
                    RESELECT = pval(:);
                elseif isvector(pval) && numel(unique(pval))==numel(pval)
                    RE_p = NaN;
                    RE_r = numel(pval);
                    RESELECT = pval(:);
                elseif isempty(pval)
                    RE_p = NaN;
                    RE_r = 0;
                    RESELECT = [];
                else
                    error(message('stats:nlmefit:InvalidREParamsSelect'))
                end
                RESELECTgiven = true;
            case {'redesign','reconstdesign'}
                if isnumeric(pval) && ndims(pval)==2
                    RE_p = size(pval,1);
                    RE_r = size(pval,2);
                    REDESIGN = pval;
                else
                    error(message('stats:nlmefit:InvalidREConstDesign'))
                end
                REDESIGNgiven = true;
            case 'regroupdesign'
                if isnumeric(pval) && size(pval,3)==M
                    RE_p = size(pval,1);
                    RE_r = size(pval,2);
                    REGROUPDESIGN = pval;
                else
                    error(message('stats:nlmefit:InvalidREGroupDesign'))
                end
                REGROUPDESIGNgiven = true;
            case 'reobsdesign'
                if isnumeric(pval) && size(pval,3)==N
                    RE_p = size(pval,1);
                    RE_r = size(pval,2);
                    REOBSDESIGN = pval;
                else
                    error(message('stats:nlmefit:InvalidREObsDesign'))
                end
                REOBSDESIGNgiven = true;
            case 'feparamsselect'
                if isvector(pval) && islogical(pval)
                    FE_p = numel(pval);
                    FESELECT = pval(:);
                elseif isvector(pval) &&  numel(unique(pval))==numel(pval)
                    FE_p = NaN;
                    FESELECT = pval(:);
                else
                    error(message('stats:nlmefit:InvalidFEParamsSelect'))
                end
                FESELECTgiven = true;
            case {'fedesign','feconstdesign'}
                if isnumeric(pval) && ndims(pval)==2
                    FE_p = size(pval,2);
                    FEDESIGN = pval;
                else
                    error(message('stats:nlmefit:InvalidFEConstDesign'))
                end
                FEDESIGNgiven = true;
            case 'fegroupdesign'
                if isnumeric(pval) && size(pval,3)==M
                    FE_p = size(pval,2);
                    FEGROUPDESIGN = pval;
                else
                    error(message('stats:nlmefit:InvalidFEGroupDesign'))
                end
                FEGROUPDESIGNgiven = true;
            case 'feobsdesign'
                if isnumeric(pval) && size(pval,3)==N
                    FE_p = size(pval,2);
                    FEOBSDESIGN = pval;
                else
                    error(message('stats:nlmefit:InvalidFEObsDesign'))
                end
                FEOBSDESIGNgiven = true;
            case 'refinebeta0'
                if islogical(pval) && isscalar(pval)
                    RefineBeta0 = pval;
                elseif ~ischar(pval)
                    error(message('stats:nlmefit:InvalidRefineBeta0'))
                else
                    ok = {'on','off'};
                    okv = find(strncmpi(pval,ok,numel(pval)));
                    if numel(okv)==1
                        RefineBeta0 = okv==1;
                    else
                        error(message('stats:nlmefit:InvalidRefineBeta0'))
                    end
                end
            case 'approximationtype'
                if ~ischar(pval)
                    error(message('stats:nlmefit:InvalidApproximationtype'))
                end
                ok = {'LME','RELME','FO','FOCE'};
                okv = find(strncmpi(pval,ok,numel(pval)));
                if numel(okv)==1
                    Approximationtype = ok{okv};
                elseif strncmpi(pval,'FO',2)
                    Approximationtype = 'FO';
                else
                    error(message('stats:nlmefit:InvalidApproximationtype'))
                end
            case 'vectorization'
                if ~ischar(pval)
                    error(message('stats:nlmefit:InvalidVectorization'))
                end
                ok = {'singlephi','singlegroup','full'};
                okv = find(strncmpi(pval,ok,numel(pval)));
                if numel(okv)==1
                    modelVectorization = ok{okv};
                else
                    error(message('stats:nlmefit:InvalidVectorization'))
                end
            case 'covparameterization'
                if ~ischar(pval)
                    error(message('stats:nlmefit:InvalidCovParameterization'))
                end
                ok = {'logm','logmcov','chol','cholcov'};
                okv = find(strncmpi(pval,ok,numel(pval)),1);
                % partial matches favor to logm and chol options which
                % are the ones documented for 9a
                if numel(okv)==1
                    CovParameterization = ok{okv};
                else
                    error(message('stats:nlmefit:InvalidCovParameterization'))
                end
            case 'options'
                try
                    Options = statset(Options,pval);
                catch ME
                    error(message('stats:nlmefit:InvalidOptions'))
                end
                if isa(Options.OutputFcn,'function_handle')
                    Options.OutputFcn = {Options.OutputFcn};
                end
            case 'optimfun'
                if ~ischar(pval)
                    error(message('stats:nlmefit:InvalidOptimFun'))
                end
                ok = {'fminsearch','fminunc'};
                okv = find(strncmpi(pval,ok,numel(pval)));
                if numel(okv)==1
                    Optimfun = ok{okv};
                else
                    error(message('stats:nlmefit:InvalidOptimFun'))
                end
                if strcmp(Optimfun,'fminunc') && isempty(ver('Optim'))
                    error(message('stats:nlmefit:NoOptim'))
                end
            case 'outputfcn'
                if iscell(pval) && all(cellfun(@(x) isa(x,'function_handle'),pval))
                    OutputFcn = pval;
                elseif isa(pval,'function_handle')
                    OutputFcn = {pval};
                elseif isempty(pval)
                    OutputFcn = {};
                else
                    error(message('stats:nlmefit:invalidOutputFcn'))
                end
            case 'paramtransform'
                ParamTransform = pval;
            case 'errormodel'
                if ~ischar(pval)
                    error(message('stats:nlmefit:InvalidErrorModel'))
                end
                ok = {'constant','proportional','combined','exponential'};
                okv = find(strncmpi(pval,ok,numel(pval)));
                if numel(okv)==1
                    ErrorModel = ok{okv};
                else
                    error(message('stats:nlmefit:InvalidErrorModel'))
                end
            case 'errorparameters'
                if numel(pval)<1 || numel(pval)>2 || ~isnumeric(pval)
                    error(message('stats:nlmefit:BadErrorParam'))
                end
                ErrorParameters = pval;
        end
    end
end

% Coalece OutputFcn and Options.OutputFcn
OutputFcn = [Options.OutputFcn(:);OutputFcn(:)];

% Check input parameters
if sum([RESELECTgiven,REDESIGNgiven,REGROUPDESIGNgiven,REOBSDESIGNgiven])>1
    opts = {'REParamsSelect' 'REConstDesign' 'REGroupDesign' 'REObsDesign'};
    opts = opts([RESELECTgiven,REDESIGNgiven,REGROUPDESIGNgiven,REOBSDESIGNgiven]);
    msg = sprintf('%s, ',opts{:});
    msg = msg(1:end-2);
    error(message('stats:nlmefit:MultipleREDesign', msg));
end

if sum([FESELECTgiven,FEDESIGNgiven,FEGROUPDESIGNgiven,FEOBSDESIGNgiven])>1
    opts = {'FEParamsSelect' 'FEConstDesign' 'FEGroupDesign' 'FEObsDesign'};
    opts = opts([FESELECTgiven,FEDESIGNgiven,FEGROUPDESIGNgiven,FEOBSDESIGNgiven]);
    msg = sprintf('%s, ',opts{:});
    msg = msg(1:end-2);
    error(message('stats:nlmefit:MultipleFEDesign', msg));
end

% the number of fixed effects (q) is always inferred from Beta0, which is
% mandatory, so we just check that FE specification complies
if FESELECTgiven
    if islogical(FESELECT)
        if sum(FESELECT)~=q
            error(message('stats:nlmefit:FinFEParamsSelectLogical'))
        end
        A = accumarray([find(FESELECT(:)),(1:q)',],1,[FE_p,q]);
    else
        if numel(FESELECT)~=q
            error(message('stats:nlmefit:FinFEParamsSelectIndex'))
        end
        A = accumarray([FESELECT(:),(1:q)'],1);
        if ~isempty(ParamTransform) && length(ParamTransform)>size(A,1)
            % We have guessed at the number of parameters, but if we have a
            % transform vector we can determine the correct number
            A(length(ParamTransform),1) = 0;
        end
    end
elseif FEDESIGNgiven
    if size(FEDESIGN,2) ~=q
        error(message('stats:nlmefit:conflictingFinFEConstDesign'))
    end
    A = FEDESIGN;
elseif FEGROUPDESIGNgiven
    if size(FEGROUPDESIGN,2) ~=q
        error(message('stats:nlmefit:conflictingFinFEGroupDesign'))
    end
    A = FEGROUPDESIGN;
elseif FEOBSDESIGNgiven
    if size(FEOBSDESIGN,2) ~=q
        error(message('stats:nlmefit:conflictingFinFEObsDesign'))
    end
    A = FEOBSDESIGN;
else % no FE specified, then use default
    A = eye(q);
    FE_p = q;
end

% The number of parameters (p) was inferred from A or the transformation
% vector. This is ambiguous when there is no transformation vector and
% FESELECT was a vector of indices, so check that B does not conflict.
p = size(A,1);
if RESELECTgiven
    if islogical(RESELECT)
        if RE_p~=p
            if isnan(FE_p)
                if RE_p<p
                    error(message('stats:nlmefit:PinFEParamsSelectIndex', RE_p));
                else
                    A(end+1:RE_p,:) = 0;
                    p = size(A,1);
                end
            else
                error(message('stats:nlmefit:PinREParamsSelectLogical', p))
            end
        end
        B = accumarray([find(RESELECT(:)),(1:RE_r)',],1,[p,RE_r]);
    else
        if ~all(ismember(RESELECT,1:p))
            if isnan(FE_p)
                p = max(RESELECT);
            else
                error(message('stats:nlmefit:PinREParamsSelectIndex', p))
            end
        end
        B = accumarray([RESELECT(:),(1:numel(RESELECT))'],1,[p,RE_r]);
    end
elseif REDESIGNgiven
    if RE_p ~=p
        if isnan(FE_p)
            if RE_p<p
                error(message('stats:nlmefit:PinFEParamsSelectIndex', RE_p));
            else
                A(end+1:RE_p,:) = 0;
                p = size(A,1);
            end
        else
            error(message('stats:nlmefit:PinREConstDesign', p))
        end
    end
    B = REDESIGN;
elseif REGROUPDESIGNgiven
    if RE_p ~=p
        if isnan(FE_p)
            if RE_p<p
                error(message('stats:nlmefit:PinFEParamsSelectIndex', RE_p));
            else
                A(end+1:RE_p,:) = 0;
                p = size(A,1);
            end
        else
            error(message('stats:nlmefit:PinREGroupDesign', p))
        end
    end
    B = REGROUPDESIGN;
elseif REOBSDESIGNgiven
    if RE_p ~=p
        if isnan(FE_p)
            if RE_p<p
                error(message('stats:nlmefit:PinFEParamsSelectIndex', RE_p));
            else
                A(end+1:RE_p,:) = 0;
                p = size(A,1);
            end
        else
            error(message('stats:nlmefit:PinREObsDesign', p))
        end
    end
    B = REOBSDESIGN;
else % no RE specified, then use default
    B = eye(p);
end

% the number of random effects (r) was inferred in all cases from B except
% when FESELECT was a vector of indices, check that B does not conflict
r = size(B,2);
if isnan(PSI_r) % default PSI
    COVPATTERN = eye(r);
elseif PSI_r ~= r
    error(message('stats:nlmefit:RinCovPattern', r))
end

% Change COVPATTERN to a linear index selector of the upper triangular part
patternCov = find(triu(COVPATTERN));
thetaLen = numel(patternCov);

if size(A,2)==0
    error(message('stats:nlmefit:nofixedeffects'))
end

if isempty(ParamTransform)
    ParamTransform = zeros(1,p);
elseif ~(isnumeric(ParamTransform) && isvector(ParamTransform) && ...
        numel(ParamTransform)==p  && all(ismember(ParamTransform,0:3)))
    error(message('stats:nlmefit:BadParamTransform', p));
end

% Ensure ErrorParameters is consistent with ErrorModel
% and apply default values if empty.
switch ErrorModel
    case 'combined'
        if isempty(ErrorParameters)
            ErrorParameters = [1 1];
        elseif numel(ErrorParameters)~=2
            error(message('stats:nlmefit:BadCombinedParam', ErrorModel));
        end
    case {'proportional', 'constant', 'exponential'}
        if isempty(ErrorParameters)
            ErrorParameters = 1;
        elseif numel(ErrorParameters)~=1
            error(message('stats:nlmefit:BadErrorParam1', ErrorModel))
        end
end
end % function parseInVarargin

function e = error_ab(ab,y,f)
g = abs(ab(1))+abs(ab(2))*abs(f);
e = sum( 0.5*((y-f)./g).^2 + log(g) );
end % function error_ab

function stop = pnlsOutputFcn(outputFcn,Beta,ofopt,inneriter,innerstate)
ofopt.inner = struct('procedure','PNLS','iteration',inneriter,'state',innerstate);
stop = false;
for i = 1:numel(outputFcn) % call each output function
    stop = stop | outputFcn{i}(Beta,ofopt,'iter');
end
end % function pnlsOutputFcn

function stop = lmeOutputFcn(optimX,optimStruc,optimState,outputFcn,ofopt,Beta,Klike)
stop = false;
if ~strcmp(optimState,'interrupt') % Here we only care about the states init, iter and done
    ofopt.inner = struct('procedure','LME','iteration',optimStruc.iteration,'state',optimState);
    ofopt.fval = Klike-optimStruc.fval;
    ofopt.theta = optimX;
    for i = 1:numel(outputFcn) % call each output function
        stop = stop | outputFcn{i}(Beta,ofopt,'iter');
    end
end
end % function lmeOutputFcn

function stop = plmOutputFcn(optimX,optimStruc,optimState,outputFcn,ofopt,Klike,q)
stop = false;
if ~strcmp(optimState,'interrupt') % Here we only care about the states init, iter and done
    ofopt.inner = struct('procedure','PLM','iteration',optimStruc.iteration,'state',optimState);
    ofopt.fval = Klike-optimStruc.fval;
    ofopt.theta = optimX(q+1:end);
    Beta = optimX(1:q);
    for i = 1:numel(outputFcn) % call each output function
        stop = stop | outputFcn{i}(Beta,ofopt,'iter');
    end
end
end % function plmOutputFcn

function stop = callOutputFcns(outputFcn,Beta,ofopt,state)
% call each output function
stop = false;
for i = 1:numel(outputFcn)
    stop = stop | outputFcn{i}(Beta,ofopt,state);
end
end % function callOutputFcns

function output2 = fevalSecondOutput(fcn, varargin)
% Helper function to get second output from a function handle.
[~, output2] = fcn(varargin{:});
end % function fevalSecondOutput

function [xOpt,fOpt,exitFlag] = robustMinimizer(fcn,x,minMethod,minMethodOpt)
% A minimizer that does nothing when the parameter vector x is empty.
if isempty(x)
    xOpt = x;
    fOpt = fcn(x);
    exitFlag = 1;
else
    [xOpt, fOpt, exitFlag] = feval(minMethod, fcn, x, minMethodOpt);
end
end % function robustMinimizer
