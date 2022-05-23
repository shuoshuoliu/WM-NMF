function [nll,VarParam] = mvregresslike(X,Y,Param,Covar,EstMethod,VarType,VarFormat)
%MVREGRESSLIKE Negative log-likelihood for multivariate regression.
%   NLOGL=MVREGRESSLIKE(X,Y,BETA,SIGMA,ALG) computes the negative log-
%   likelihood function NLOGL for a multivariate regression of the multivariate
%   observations in the N-by-D matrix Y on the predictor variables in X,
%   evaluated for the P-by-1 column vector BETA of coefficient estimates
%   and the D-by-D matrix SIGMA specifying the covariance of a row of Y.
%   For any value of D>1 with X being an N-by-P matrix, the coefficient estimates
%   BETA should be a P-by-D matrix. 
%
%   ALG specifies the algorithm used in fitting the regression (see below).
%   NaN values in X or Y are taken to be missing.  Observations with
%   missing values in X are ignored.  Treatment of missing values in Y
%   depends on the algorithm.
%
%   Y is an N-by-D matrix of D-dimensional multivariate observations.  X
%   may be either a matrix or a cell array.  If D=1, X may be an N-by-P
%   design matrix of predictor variables.  For any value of D, X may also
%   be a cell array of length N, each cell containing a D-by-P design
%   matrix for one multivariate observation.  If all observations have the
%   same D-by-P design matrix, X may be a single cell. For any value of D, 
%   X may also be an N-by-P matrix. 
%
%   ALG should match the algorithm used by MVREGRESS to obtain the
%   coefficient estimates BETA, and must be one of the following values:
%         'ecm'    ECM algorithm
%         'cwls'   Least squares conditionally weighted by SIGMA
%         'mvn'    Multivariate normal estimates computed after omitting
%                  rows with any missing values in Y.
%
%   [NLOGL,VARPARAM]=MVREGRESSLIKE(...) also returns an estimated covariance
%   matrix of the parameter estimates BETA.
%
%   [NLOGL,VARPARAM]=MVREGRESSLIKE(...,VARTYPE,VARFORMAT) specifies the type
%   and format of VARPARAM.  VARTYPE is either 'hessian' (default) to use
%   the Hessian or observed information, or 'fisher' to use the Fisher or
%   expected information.  VARFORMAT is either 'beta' (default) to compute
%   VARPARAM for BETA only, or 'full' to compute VARPARAM for both BETA and
%   SIGMA. The 'hessian' method takes into account the increased
%   uncertainties due to missing data.  The 'fisher' method uses the
%   complete-data expected information, and does not include uncertainty
%   due to missing data.
%
%   See also MVREGRESS, REGSTATS, MANOVA1.

%    Copyright 2006-2013 The MathWorks, Inc.

  
if nargin > 4
    EstMethod = convertStringsToChars(EstMethod);
end

if nargin > 5
    VarType = convertStringsToChars(VarType);
end

if nargin > 6
    VarFormat = convertStringsToChars(VarFormat);
end

narginchk(4,7);

if isempty(X) || (iscell(X) && isempty(X{1}))
    error(message('stats:mvregresslike:EmptyDesignArray'));
end
if isempty(Y)
    error(message('stats:mvregresslike:EmptyDataArray'));
end
if isempty(Param)
    error(message('stats:mvregresslike:EmptyParam'));
end
if isempty(Covar)
    error(message('stats:mvregresslike:EmptyCovar'));
end

ny = size(Y,2);
if ~isequal(size(Covar),[ny ny])
    error(message('stats:mvregresslike:InconsistentDims', ny, ny));
elseif any(any(isnan(Covar)))
    error(message('stats:mvregresslike:NanCovar'));
end

[NumSamples, NumSeries] = size(Y);
if ~iscell(X)&&(NumSeries~=1)
    [n m] = size(X);
    if n ~= NumSamples
        error(message('stats:statcheckmvnr:InconsistentXRows', NumSamples));
    end
    if ~isequal(size(Param),[m NumSeries])
        error(message('stats:mvregresslike:BadBeta1', m, NumSeries));
    end
    Param = reshape(Param,[m*NumSeries 1]);
    c = cell(NumSamples,1);
    for j = 1:NumSamples
        c{j} = kron(eye(NumSeries),X(j,:));
    end
    X = c;
end

if ~isvector(Param)
    error(message('stats:mvregresslike:BadBeta'));
end
Param = Param(:);

if nargin<5
    EstMethod = '';  % use data-dependent default
end
if ~isempty(EstMethod)
    if ~ischar(EstMethod) || size(EstMethod,1)~=1
        EstMethod = [];
    else
        [~,EstMethod] = internal.stats.getParamVal(EstMethod,{'cwls' 'ecm' 'mvn'},'ALG');
    end
end
if ~(isempty(EstMethod) || ismember(EstMethod,1:3))
    error(message('stats:mvregresslike:BadAlgorithm'));
end

if nargin < 6 || isempty(VarType)
    VarType = 'hessian';
else
    % Check variance type
    okvals = {'hessian' 'fisher'};
    if ~ischar(VarType) || size(VarType,1)~=1
        VarType = [];
    else
        VarType = internal.stats.getParamVal(VarType,okvals,'VARTYPE');
    end
end

if nargin < 7 || isempty(VarFormat)
    VarFormat = 'paramonly';  % internal code for 'beta'
else
    % Check variance format
    okvals = {'beta' 'full'};
    if ~ischar(VarFormat) || size(VarFormat,1)~=1
        VarFormat = [];
    else
        [~,VarFormat] = internal.stats.getParamVal(VarFormat,okvals,'VARFORMAT');
    end
    okvals{1} = 'paramonly';
    VarFormat = okvals{VarFormat};
end

% Check inputs, ignoring NaN rows for mvn method (3)
if isempty(EstMethod) && ~any(isnan(Y(:)))
    EstMethod = 3;   % use faster method for cases where others are equivalent
end
[~,~,~,Y,X] = statcheckmvnr(Y, X, Param, Covar, isequal(EstMethod,3));

d = diag(Covar);
isdiagonal = isequal(Covar,diag(d));

if isdiagonal        % use faster formula for diagonal covariance
    CholState = sum(d<=0);
    CholCovar = diag(sqrt(d));
else
    [CholCovar, CholState] = chol(Covar);
end
if CholState > 0
    error(message('stats:mvregresslike:NonPosDefCov'));
end


if EstMethod==3                   % mvn method
    nll = -statmvnrobj(Y,X,Param,Covar,[],CholCovar,isdiagonal);
else                              % cwls or ecm method
    nll = -statecmobj(X,Y,Param,Covar,[],CholCovar,isdiagonal);
end

if nargout>=2
    Info = statecmmvnrfish(Y,X,Covar,VarType,VarFormat,CholCovar,isdiagonal);
    VarParam = inv(Info);
end
