function  out = nlmeres(regFcn, jacFcn, y, groupMatches, omega, errorModel, errorparam, randomEffects)
%NLMERES function to calculate residuals for nlmefit/nlmefitsa estimation.
%   OUT = NLMERES(REGFCN, JACFCN, Y, GROUPNUMBERS, OMEGA, ERRORMODEL, ERRORPARAM, B)
%   returns a structure OUT with the followign fields:
%       PRES        The population residuals (y - y_individual)
%       IRES        The individual residuals (y - y_population)
%       PWRES       The population weighted residuals
%       CWRES       The conditional weighted residuals
%       IWRES       The individual weighted residuals
%
% REGFCN is the regression function containing the structural model. It has
% the following signature:
%
%   YPRED = REGFCN(B)
%
% where B is a matrix containing the random effects for all groups.
%
% Y is a vector containing the observed values.
%
% JACFCN is the Jacobian of REGFCN with respect to B. It has the following
% signature:
%
%   JB = JACFCN(F0, B)
%
% where F0 is REGFCN(B) and B is a matrix of random effects for all groups.
% JB is matrix where element J(i,j) is the derivative of Y(i) with respect
% to B(j).
%
% GROUPMATCHES is a logical matrix, with element (i,j) indicating whether
% observation i belongs to group j.
%
% OMEGA is the covariance matrix for the random effects.
%
% ERRORMODEL is a scalar identifying the error model.
%
% ERRORPARAM is the estimated parameters of the error variance model.
%
% B is a r X n matrix of 'r' random effects for 'n' individuals.
%
%   See also NLMEFIT, NLMEFITSA.

%   References:
%      [1]  Hooker et al. Conditional weighted residuals (CWRES): a model
%      diagnostic for the FOCE method. Pharmaceutical research (2007) vol.
%      24 (12) pp. 2187-2197

%   Copyright 2010-2011 The MathWorks, Inc.
sigma            = 1;
nGroups = size(groupMatches, 2);
nB = size(omega,1);

out = struct('pres', [], 'ires', [], 'pwres', [], 'iwres', [], 'cwres', []);
out.pwres        = zeros(nGroups,1);
out.iwres        = zeros(nGroups,1);
out.cwres        = zeros(nGroups,1);

randomEffectZero = zeros(nB,nGroups);

y_pop = regFcn(randomEffectZero);
y_ind = regFcn(randomEffects);

out.pres = y - y_pop;
out.ires = y - y_ind;

% Calculate the jacobian of the strutural model with respect to random.
% effects.
df_db_0 = jacFcn(y_pop, randomEffectZero);
df_db_bhat = jacFcn(y_ind, randomEffects);

% Calculate dh/de, here h is the function for the error model.
dh_de = get_dh_de(errorModel, errorparam, y_pop);

for i = 1:nGroups
    idx = groupMatches(:,i);
    
    df_db_0_i = df_db_0(idx,:);
    df_db_bhat_i = df_db_bhat(idx,:);
    dh_de_i = dh_de(idx,:);
    
    % WRES =  (y - E(y))/R, where R is the cholesky factor of the covariance
    % matrix of y, covarianceMatrixOfY.
    % covarianceMatrixOfY = df_db_0*omega*df_db_0' + diag(diag(dh_de*sigma*dh_de')
    % where df_db_0 is the Jacobian of the model function at randomEffects =
    % 0 and dh_de is the Jacobian of the error function at e = 0. For CWRES
    % we use df_db_bhat, the Jacobian of the model function with respect to
    % random effects at individual estimates of the random effects. In the
    % following lines we compute the terms required to compute weighted
    % residuals.
    
    % Compute Covariance matrix. omega(r x r) is assumed to be positive
    % definite. df_db_0*omega*df_db_0'(m x m) may not positive definite (if
    % the number of time points for an individual is more than then number
    % of random effects). Adding diag(diag(dh_de*sigma*dh_de')), which is
    % the variance of error term, makes covarianceMatrixOfY Positive
    % definite. But if the variance of the error is too small,
    % covarianceMatrixOfY will not be positive definite.
    varianceOfError_i = diag(diag(dh_de_i*sigma*dh_de_i'));
    covarianceMatrixOfY_i = df_db_0_i*omega*df_db_0_i' + varianceOfError_i;
    
    % Compute PWRES and IWRES.
    [R, p]  = chol(covarianceMatrixOfY_i);
    if p>0
        out.pwres = [];
        out.cwres = [];
        out.iwres = [];
        warning(message('stats:nlmeres:ErrorVarianceCloseToZero'));
        return;
    end
    out.pwres(idx) = out.pres(idx)'*linsolve(R, eye(size(R)), struct('UT',true));
    out.iwres(idx) = (out.ires(idx)./dh_de_i)'; % dh_de_i is same as standard deviation of error.
    
    % Calculate CWRES.
    covarianceMatrixOfY_i = df_db_bhat_i*omega*df_db_bhat_i' + varianceOfError_i;
    [R, p] = chol(covarianceMatrixOfY_i);
    if p>0
        out.pwres = [];
        out.cwres = [];
        out.iwres = [];
        warning(message('stats:nlmeres:ErrorVarianceCloseToZero'));
        return;
    end
    
    % Theoretical expectation of y, based on FOCE approximation, is
    expected_y_i = y_ind(idx) - df_db_bhat_i*randomEffects(:,i);
    out.cwres(idx) = (y(idx) - expected_y_i)'*linsolve(R, eye(size(R)), struct('UT',true));
end
end

function dh_de = get_dh_de(errorModel, errorparam, y)

switch errorModel
    case 'constant'
        dh_de = repmat(errorparam, numel(y), 1);
        
    case 'proportional'
        dh_de = y*errorparam;
        
    case 'combined'
        dh_de = errorparam(1) + errorparam(2)*abs(y);
        
    case 'exponential'
        dh_de = repmat(errorparam, numel(y), 1);
        
    otherwise
        error(message('stats:internalerror:InvalidErrorModel'));
end
end
