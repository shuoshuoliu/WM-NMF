function [ypred, delta, sch] = nlpredci(model,X,beta,resid,varargin)
%NLPREDCI Confidence intervals for predictions in nonlinear regression.
%   [YPRED, DELTA] = NLPREDCI(MODELFUN,X,BETA,RESID,'covar',SIGMA) returns
%   predictions (YPRED) and 95% confidence interval half-widths (DELTA)
%   for the nonlinear regression model defined by MODELFUN, at input values X.
%   MODELFUN is a function, specified using @, that accepts two arguments,
%   a coefficient vector and the array X, and returns a vector of fitted Y
%   values.  Before calling NLPREDCI, use NLINFIT to fit MODELFUN by
%   nonlinear least squares and get estimated coefficient values BETA,
%   residuals RESID, and estimated coefficient covariance matrix SIGMA.
%
%   [YPRED, DELTA] = NLPREDCI(MODELFUN,X,BETA,RESID,'jacobian',J) is an
%   alternative syntax that also computes 95% confidence intervals.  J is
%   the Jacobian computed by NLINFIT.  You should use the 'covar' input
%   rather than the 'jacobian' input if you use a robust option with
%   NLINFIT, because the SIGMA parameter is required to take the robust
%   fitting into account.
%
%   [YPRED, DELTA] = NLPREDCI(MODELFUN,X,BETA,RESID,...,'PARAM',val) allows
%   you to pass an optional name/value pair that specifies either the observation 
%   weights or an error variance model to be used for making predictions. You 
%   cannot use 'ErrorModelInfo' and 'Weights' at the same time.
%
%        Name              Value
%       'ErrorModelInfo'   A structure S returned by a previous call to nlinfit
%                          that specifies the error model to be used for making 
%                          predictions. S has the following fields: 
%
%                          ErrorModel           - Name of the error model previously chosen in
%                                                 nlinfit.
%                          ErrorParameters      - Error parameters estimated by nlinfit.
%                          ErrorVariance        - A function handle that accepts a n by p 
%                                                 matrix X and computes a n by 1 vector of 
%                                                 error variances at the n observations in X 
%                                                 for this error model.  
%                          MSE                  - Mean squared error.
%                          ScheffeSimPred       - Scheffe parameter for a simultaneous 
%                                                 prediction interval when using this error 
%                                                 model. 
%                          WeightFunction       - True if a custom weight function was used
%                                                 previously in nlinfit.
%                          FixedWeights         - True if fixed weights were used previously
%                                                 in nlinfit.
%                          RobustWeightFunction - True if a robust fitting was used previously
%                                                 in nlinfit.
%                          S can be obtained as an output by calling nlinfit with
%                          the desired error variance model. The specified 'ErrorModelInfo' 
%                          has an effect only when using 'predopt' as 'observation'. 
%                          When 'ErrorModelInfo' is not supplied, a 'constant' error 
%                          variance model will be assumed.
%
%       'Weights'          A vector W of real positive values with the same number
%                          of elements as the number of observations (or rows) in 
%                          X. The error variance at observation i is estimated as 
%                          MSE * (1/W(i)) where MSE is the mean squared error.               
%                          'Weights' can also be specified as a function handle that 
%                          accepts a vector of predicted response values and returns 
%                          a vector of real positive values as output. Default is no 
%                          weights.
%
%   [...] = NLPREDCI(...,'PARAM1',val1,'PARAM2',val2,...) allows you to
%   specify optional parameter name/value pairs as follows:
%
%       Name         Value
%      'alpha'       A value between 0 and 1 to specify the confidence level
%                    as 100(1-ALPHA)%.  Default is 0.05.
%      'mse'         The mean squared error returned by nlinfit.  This is
%                    required to predict new observations (see 'predopt')
%                    if you used a robust option with nlinfit; otherwise
%                    the mse is computed from the residuals and does not
%                    take the robust fitting into account.
%      'predopt'     Either 'curve' (the default) to compute confidence
%                    intervals for the estimated curve (function value) at
%                    X, or 'observation' for prediction intervals for a new
%                    observation at X.  If you specify 'observation' after
%                    using a robust option with nlinfit, you must also supply
%                    the 'mse' parameter to specify the robust estimate of
%                    the mean squared error.
%      'simopt'      Either 'on' for simultaneous bounds, or 'off' (the default)
%                    for non-simultaneous bounds.
%
%   NLPREDCI treats NaNs in RESID or J as missing values, and ignores the
%   corresponding observations.
%
%   When J does not have full column rank, some of the model parameters may
%   be non-identifiable. In this case, nlpredci will try to construct confidence
%   intervals for those predictions that are estimable and return NaN for those
%   that are not estimable. 
%
%   Example:
%      load reaction;
%      [beta,resid,J,Sigma] = nlinfit(reactants,rate,@hougen,beta);
%      newX = reactants(1:2,:);
%      [ypred, delta] = nlpredci(@hougen,newX,beta,resid,'Covar',Sigma);
%
%   See also NLINFIT, NLPARCI, NLINTOOL.

% Older syntax still supported:
%   [YPRED, DELTA] = NLPREDCI(FUN,X,BETA,RESID,J,ALPHA,SIMOPT,PREDOPT)

%   References:
%      [1] Seber, G.A.F, and Wild, C.J. (1989) Nonlinear Regression, Wiley.

%   To compute confidence intervals when the parameters or data are complex,
%   you will need to split the problem into its real and imaginary parts.
%   First, define your parameter vector BETA as the concatenation of the real
%   and imaginary parts of the original parameter vector.  Then concatenate the
%   real and imaginary parts of the response vector Y as a single vector.
%   Finally, modify your model function MODELFUN to accept X and the purely
%   real parameter vector, and return a concatenation of the real and
%   imaginary parts of the fitted values.  Given this formulation of the
%   problem, NLINFIT will compute purely real estimates, and confidence
%   intervals are feasible.

%   Copyright 1993-2011 The MathWorks, Inc.


% Set parameter defaults.
if nargin > 0
    model = convertStringsToChars(model);
end

if nargin > 4
    [varargin{:}] = convertStringsToChars(varargin{:});
end

J = []; Sigma = []; alpha = 0.05; mse = []; simopt = 'off'; predopt = 'curve';
errorModelInfo = []; weights = [];
    
% Parse input name/value pairs - respect old syntax.
if nargin < 5
    error(message('stats:nlpredci:TooFewInputs'));
end
if nargin>=5 && ischar(varargin{1})
    % Calling sequence with named arguments
    okargs =   {'jacobian','covariance','alpha','mse','simopt','predopt','ErrorModelInfo','Weights'};
    defaults = {[]        , []         ,  0.05 ,  [] ,   'off',  'curve',        []      , weights };
    [J,Sigma,alpha,mse,simopt,predopt,errorModelInfo,weights] = ...
                         internal.stats.parseArgs(okargs,defaults,varargin{:});
else
    % [YPRED, DELTA] = NLPREDCI(FUN,X,BETA,RESID,J,ALPHA,SIMOPT,PREDOPT)
    if nargin>=5, J = varargin{1}; end
    if nargin>=6, alpha = varargin{2}; end
    if nargin>=7, simopt = varargin{3}; end
    if nargin>=8, predopt = varargin{4}; end
end
    
% Validate alpha.
if isempty(alpha)
    % For a 95% confidence interval.
    alpha = 0.05;
elseif (~isscalar(alpha) || alpha<=0 || alpha >= 1)
    error(message('stats:nlpredci:BadAlpha'));
end

% Set default simopt and predopt.
if isempty(simopt)
    simopt = 'off'; 
end
if isempty(predopt)
    predopt = 'curve'; 
end

% Validate simopt and predopt.
switch(simopt)
 case 'on', dosim = 1;
 case 'off', dosim = 0;
 otherwise, error(message('stats:nlpredci:BadSimOpt'));
end
switch(predopt)
 case {'c' 'curve'}, newobs = 0;
 case {'o' 'observation'}, newobs = 1;
 otherwise, error(message('stats:nlpredci:BadPredOpt'));
end

% Make sure we have everything we need.
if isempty(resid) || (isempty(J) && isempty(Sigma))
   error(message('stats:nlpredci:TooFewInputsResid'));
end

% Make sure input is not complex.
if ~isreal(beta) || ~isreal(J) || ~isreal(resid) || ~isreal(Sigma)
    error(message('stats:nlpredci:ComplexParams'));
end

% Remove missing values from resid and get:
%   n = effective number of points used in nlinfit and
%   p = length of parameter vector
% Also, ensure that n > p.
resid = resid(:);
missing = isnan(resid);
resid(missing) = [];
n = numel(resid);
p = numel(beta);
if n <= p
   error(message('stats:nlpredci:NotEnoughData'));
end

% Odds are, an input of length n should be a column vector.
if (size(X,1)==1 && size(X,2)==n)
    X = X(:); 
end

% If required, remove missing values from J and check its size.
if ~isempty(J)
   J(missing,:) = [];
   if size(J,1)~=n || size(J,2)~=p
      error(message('stats:nlpredci:InputSizeMismatch'));
   end
   % Approximation when a column is zero vector
   temp = find(max(abs(J)) == 0);
   if ~isempty(temp)
      J(:,temp) = sqrt(eps(class(J)));
   end
end

% If mse is given, make sure it is sensible.
if ~isempty(mse) && ( ~isnumeric(mse) || ~isreal(mse) || ~isscalar(mse) || (mse < 0) )
    error(message('stats:nlpredci:BadMSE'));
end

% If required, make sure Sigma has the right size.
if ~isempty(Sigma)
   if size(Sigma,1) ~= p || size(Sigma,2) ~= p
       error(message('stats:nlpredci:SigmaSizeMismatch'));
   end
end

% Compute the predicted values at the new X.
ypred = feval(model,beta,X);
if ~isreal(ypred)
    error(message('stats:nlpredci:ComplexY'));
end

% Make sure ypred is a vector.
if ~isvector(ypred)
    error(message('stats:nlpredci:YPredNotVector'));
end

% Do we have an exponential error model? 
% If so, as done in nlinfit, transform it to a constant error model.
exponentialErrorModel = false;
if ( ~isempty(errorModelInfo) && strcmpi(errorModelInfo.ErrorModel,'exponential') )
    % Transform model.
    [ypred, model] = applyLogTransformation(ypred, model);
    exponentialErrorModel = true;
end

% If errorModelInfo is not empty, make sure it is sensible.
if ~isempty(errorModelInfo)
    errorModelInfo = checkErrorModelInfo(errorModelInfo,X,ypred);
end

% Check user supplied weights if required.
if ~isempty(weights)
    weights = checkWeights(weights, ypred);
end

% Enforce that both weights and errorModelInfo cannot be supplied.
if ~isempty(errorModelInfo) && ~isempty(weights)
    error(message('stats:nlpredci:WeightsErrorModelInfoConflict'));
end

% Approximate the Jacobian at the new X.
delta = zeros(length(ypred),numel(beta));
fdiffstep = eps(class(beta)).^(1/3);
for i = 1:length(beta)
   change = zeros(size(beta));
   if (beta(i) == 0)
      nb = sqrt(norm(beta));
      change(i) = fdiffstep * (nb + (nb==0));
   else
      change(i) = fdiffstep*beta(i);
   end
   predplus = feval(model, beta+change, X);
   delta(:,i) = (predplus - ypred)/change(i);
end

% Try to use J instead of Sigma whenever possible.
if ~isempty(J)    
    usingJ = true;% Use J.
else    
    usingJ = false;% Use Sigma.
end

% Mark if J or Sigma is ill-conditioned.
illConditioned = false;
if usingJ
    [~,R] = qr(J,0);
    if condest(R) > 1/sqrt(eps(class(beta)))
        illConditioned = true;    
    end
else
    if condest(Sigma) > 1/eps(class(beta))
        illConditioned = true;    
    end
end

% Get rankJ and mark estimable contrasts in rows of delta. If usingJ, also 
% get pinvJTJ for later use in calculating Sigma.
if illConditioned    
    % Set TolSVD for isEstimable.
    TolSVD = eps(class(beta));
    if usingJ
        % J is given.
        [estimable,rankJ,pinvJTJ] = internal.stats.isEstimable(delta, 'DesignMatrix', J,'TolSVD',TolSVD);   
    else
        % Sigma is given. Sigma = mse*pinv(J'*J). We will use Sigma as
        % pinvNormalMatrix instead of Sigma/mse as this should not affect
        % estimability or rank of J.
        [estimable,rankJ] = internal.stats.isEstimable(delta, 'pinvNormalMatrix', Sigma,'TolSVD',TolSVD);   
    end    
else
    % Use existing QR factorization, if required.
    if usingJ
        % J is given.
        rankJ = size(J,2);
        Rinv = R \ eye(size(R));
        pinvJTJ = (Rinv*Rinv');
    else
        % Sigma is given.
        rankJ = size(Sigma,2);
    end
    estimable = true(size(delta,1),1);
end

% Get MSE. The degrees of freedom when J is full rank is v = n-p and n-rank(J) otherwise.
if isempty(mse) && ~isempty(errorModelInfo) && ~isempty(errorModelInfo.MSE)
    % Get MSE from errorModelInfo if possible.
    mse = errorModelInfo.MSE;
elseif isempty(mse) && (isempty(Sigma) || newobs)
    % If we really need mse and errorModelInfo is not given, do this.
    mse = norm(resid)^2 / (n-rankJ);        
end

% Calculate Sigma if usingJ only if Sigma is currently empty. 
if usingJ && isempty(Sigma)
   Sigma = mse*pinvJTJ;
end

% Compute varpred.
varpred = sum((delta*Sigma) .* delta,2);
% Are we using an error model to compute errorVar?
if (newobs)
    if ~isempty(weights)
        if isa(weights,'function_handle')
            wVec = weights(ypred);
        else
            wVec = weights;
        end
        % Assume error variance is mse*(1./wVec).
        wVec = wVec(:);
        errorVar = mse./wVec;
    elseif ~isempty(errorModelInfo) && ~isempty(errorModelInfo.ScheffeSimPred) && ~isempty(errorModelInfo.ErrorVariance)
        % Use ErrorVariance function in errorModelInfo.
        errorVar = errorModelInfo.ErrorVariance(X);
        errorVar = errorVar(:);
    else
        % Assume a constant variance model if errorModelInfo and weights are 
        % not supplied.
        errorVar = mse * ones(size(delta,1),1);
    end    
    varpred = varpred + errorVar;
end

% Compute delta, the CI half-widths.
sch = [];
if (dosim)
   % Simultaneous CIs.
   if (newobs)
        % Simultaneous CI on new observations.
        if ~isempty(errorModelInfo) && ~isempty(errorModelInfo.ScheffeSimPred) && ~isempty(errorModelInfo.ErrorVariance)
            % Use the precomputed Scheffe parameter valid for using this error model.
            sch = errorModelInfo.ScheffeSimPred;
            % Make sure sch is either rankJ or (rankJ+1).
            if ( sch ~= rankJ && sch ~= (rankJ+1) )
                % Be conservative.
                sch = (rankJ+1);
            end            
            % Validate Scheffe parameter against Jacobian at prediction points.
            sch = validateScheffeParamUsingJacobianPred(sch, rankJ, delta, estimable, errorVar);                                    
        else
            % Be conservative.
            sch = (rankJ+1);
        end        
   else
       % Simultaneous CI on curve.
        sch = rankJ;
   end
   crit = sqrt(sch * finv(1-alpha, sch, n-rankJ));   
else
   % Pointwise CIs.
   crit = tinv(1-alpha/2,n-rankJ);
end

delta = sqrt(varpred) * crit;

% Set non-estimable contrasts in delta to NaN.
delta(~estimable) = NaN;

% delta should have the same size as ypred.
if size(ypred,1) == 1
    % ypred is a row vector, transpose delta.
    delta = delta';
end

% If we have an exponential error model, make predictions on the same
% scale as the original model function.
if exponentialErrorModel
    ypred = exp(ypred);
end

end % End of nlpredci.

%==== Subfunction applyLogTransformation ====
function [y, model] = applyLogTransformation(y, model)
    % Exponential error model. Linearize the model as
    %   y = f*exp(a*e), or log(y) = log(f) + a*e
    
    if ~all(y>0)
        error(message('stats:nlinfit:PositiveYRequired'));
    else
        y = log(max(y,realmin));
    end
    
    model = @(phi,X) log(max(model(phi,X),realmin));
end % End of applyLogTransformation.


%==== Subfunction validateScheffeParamUsingJacobianPred ====
function sch = validateScheffeParamUsingJacobianPred(sch, rankJ, delta, estimable, errorVar)
% Check if Jacobians used during model fitting and prediction have the same structure.

    if sch ~= (rankJ+1)                
        % Get only estimable rows in delta.
        delta_est = delta(estimable,:);
        % Since EVFit = ones(size(EVPred)), Unweighted Jacobian = Weighted Jacobian.
        EVPred = errorVar(estimable); EVFit = ones(size(EVPred));                
        [schEstDelta,rankEstDelta] = internal.stats.getscheffeparam('UnWeightedJacobian',delta_est,'Intopt','observation','ErrorVarianceFit',EVFit,'ErrorVariancePred',EVPred);
        if ( schEstDelta ~= rankEstDelta )
           % Prediction Jacobian has a different structure. If Jpred = Vr*Sigmar*Ur' 
           % then Vr*Vr'*dpred \neq dpred where dpred(i) = sqrt(Gpred_ii). Reset sch
           % to (rankJ+1).
           sch = (rankJ+1); 
        end        
    end
    
end % End of validateScheffeParamUsingJacobianPred.

%==== Subfunction checkWeights ====
function weights = checkWeights(weights, yfit)
% If weights is a function handle, make sure it works as required. Whether
% weights vector is given directly or through a function handle, make sure 
% it is numeric, real and positive and of the same length as the number of 
% predictions.

    if isa(weights, 'function_handle') || (isnumeric(weights) && isvector(weights))
        % If weights are set.

        if isa(weights, 'function_handle')
            % function handle.
            try
                wVec = weights(yfit);
            catch ME
                if isa(weights, 'inline')
                    m = message('stats:nlinfit:InlineWeightFunctionError');
                    throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
                elseif strcmp('MATLAB:UndefinedFunction', ME.identifier) ...
                        && ~isempty(strfind(ME.message, func2str(weights)))
                    error(message('stats:nlinfit:WeightFunctionNotFound',func2str(weights)));
                else
                    m = message('stats:nlinfit:WeightFunctionError',func2str(weights));
                    throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
                end
            end
        else
            % fixed weights.
            wVec = weights;
        end
        
        % wVec should be real positive vector of the same size as yfit.
        if ~isequal(size(wVec), size(yfit)) || ~isnumeric(wVec) || ~isreal(wVec) || ~isvector(wVec) || any( wVec(:) <= 0 ) 
            error(message('stats:nlinfit:InvalidWeights'));
        end

    else
       % Invalid weights.
       error(message('stats:nlinfit:InvalidWeights'));
    end

end % End of checkWeights.

%==== Subfunction checkErrorModelInfo ====
function errorModelInfo = checkErrorModelInfo(errorModelInfo,X,ypred)
%  Make sure the supplied errorModelInfo structure conforms to specifications.
%  errorModelInfo should have the following fields:
%
%  ErrorModel           - Name of the error model (char string).
%  ErrorParameters      - Parameters of the error model (numeric vector).
%  ErrorVariance        - Variance function for the error model (function handle).
%  MSE                  - Mean squared error (numeric scalar, > 0).  
%  ScheffeSimPred       - Scheffe parameter for simultaneous prediction intervals (numeric scalar, integer).
%  WeightFunction       - True if using 'Weights' with a weight function (logical scalar).
%  FixedWeights         - True if using 'Weights' with a fixed weight vector (logical scalar).
%  RobustWeightFunction - True if using Robust fitting with options.RobustWgtFun (logical scalar). 
    
    % Validate ErrorModel.
    allowedErrorModels = {'constant','combined','proportional','exponential'};
    if ~isfield(errorModelInfo,'ErrorModel') || ~ischar(errorModelInfo.ErrorModel) ... 
            || ~any(strcmpi(errorModelInfo.ErrorModel,allowedErrorModels))
        error(message('stats:nlpredci:InvalidErrorModel'));
    end
    
    % Validate ErrorParameters.
    if ~isfield(errorModelInfo,'ErrorParameters') || isempty(errorModelInfo.ErrorParameters) ... 
            || ~isnumeric(errorModelInfo.ErrorParameters) || ~isreal(errorModelInfo.ErrorParameters) ... 
            || ~isvector(errorModelInfo.ErrorParameters)
        error(message('stats:nlpredci:InvalidErrorParameters'));
    end
    
    % Validate ErrorVariance.
    if ~isfield(errorModelInfo,'ErrorVariance') || ~isa(errorModelInfo.ErrorVariance,'function_handle')
        error(message('stats:nlpredci:InvalidErrorVariance'));
    else
        % Ensure that ErrorVariance returns sensible results.
        errorVar = errorModelInfo.ErrorVariance(X);
        % Make errorVar and ypred into column vectors.
        errorVar = errorVar(:);
        ypred = ypred(:);
        if size(errorVar,1) ~= size(ypred,1)
            error(message('stats:nlpredci:InvalidErrorVariance'));
        end
    end
    
    % Validate MSE.
    if ~isfield(errorModelInfo,'MSE') || isempty(errorModelInfo.MSE) ... 
            || ~isnumeric(errorModelInfo.MSE) || ~isreal(errorModelInfo.MSE) ... 
            || ~isscalar(errorModelInfo.MSE) || (errorModelInfo.MSE <= 0)
    	error(message('stats:nlpredci:InvalidMSE'));
    end
    
    % Validate ScheffeSimPred.
    if ~isfield(errorModelInfo,'ScheffeSimPred') || isempty(errorModelInfo.ScheffeSimPred) ...
            || ~isnumeric(errorModelInfo.ScheffeSimPred) || ~isreal(errorModelInfo.ScheffeSimPred) ... 
            || ~isscalar(errorModelInfo.ScheffeSimPred) || errorModelInfo.ScheffeSimPred <= 0
        error(message('stats:nlpredci:InvalidScheffeSimPred'));
    else
        % Make sure ScheffeSimPred is an integer.
        errorModelInfo.ScheffeSimPred = floor(errorModelInfo.ScheffeSimPred);
    end
    
    % Validate WeightFunction.
    if ~isfield(errorModelInfo,'WeightFunction') || isempty(errorModelInfo.WeightFunction) ... 
            || ~islogical(errorModelInfo.WeightFunction) || ~isscalar(errorModelInfo.WeightFunction)
        error(message('stats:nlpredci:InvalidWeightFunctionFlag'));
    end
    
    % Validate FixedWeights.
    if ~isfield(errorModelInfo,'FixedWeights') || isempty(errorModelInfo.FixedWeights) ... 
            || ~islogical(errorModelInfo.FixedWeights) || ~isscalar(errorModelInfo.FixedWeights)
        error(message('stats:nlpredci:InvalidFixedWeightsFlag'));
    end
    
    % Validate RobustWeightFunction.
    if ~isfield(errorModelInfo,'RobustWeightFunction') || isempty(errorModelInfo.RobustWeightFunction) ... 
            || ~islogical(errorModelInfo.RobustWeightFunction) || ~isscalar(errorModelInfo.RobustWeightFunction)
        error(message('stats:nlpredci:InvalidRobustWeightFunctionFlag'));
    end
    
end % End of checkErrorModelInfo.
