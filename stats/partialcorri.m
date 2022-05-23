function [coef,pval] = partialcorri(Y,X,varargin)
%PARTIALCORRI Partial correlation coefficients with internal adjustments.
%   RHO = PARTIALCORRI(Y,X) returns the sample linear partial correlation
%   coefficients between pairs of variables between Y and X, adjusting
%   for the remaining variables in X.  Y is N-by-P matrix and X is N-by-Q
%   matrix, with rows corresponding to observations, and columns
%   corresponding to variables. RHO is a P-by-Q matrix, where RHO(I,J) 
%   is the sample linear partial correlation coefficient between Y(:,I) 
%   and X(:,J), controlling for all the columns of X except column J.
%
%   RHO = PARTIALCORRI(Y,X,Z) returns the sample linear partial correlation
%   coefficients between pairs of variables between Y and X, adjusting 
%   for the remaining variables in X, after first controlling both X and 
%   Y for the variables in Z.  Y is an N-by-P matrix, X is an N-by-Q 
%   matrix, and Z is an N-by-T matrix, with rows corresponding to observations, 
%   and columns corresponding to variables. RHO is a P-by-Q matrix.
%
%   [RHO,PVAL] = PARTIALCORRI(...) also returns PVAL, a matrix of p-values for
%   testing the hypothesis of no partial correlation against the alternative
%   that there is a non-zero partial correlation.  Each element of PVAL is the
%   p-value for the corresponding element of RHO.  If PVAL(i,j) is small, say
%   less than 0.05, then the partial correlation RHO(i,j) is significantly
%   different from zero.
%
%   [...] = PARTIALCORRI(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values.  Valid parameters are the
%   following:
%
%        Parameter  Value
%         'type'    'Pearson' (the default) to compute Pearson (linear)
%                   partial correlations or 'Spearman' to compute Spearman
%                   (rank) partial correlations.
%         'rows'    'all' (default) to use all rows regardless of missing
%                   values (NaNs), 'complete' to use only rows with no
%                   missing values.  Use 'pairwise' to use all available
%                   values in each column of Y when computing partial
%                   correlations coefficients and p-values corresponding
%                   to that column.  For each column of Y, rows will be
%                   dropped corresponding to missing values in X
%                   (and/or Z, if supplied).  However, remaining rows 
%                   with valid values in that column of Y will be used,
%                   even if there are missing values in other columns
%                   of Y.
%         'tail'    The alternative hypothesis against which to compute
%                   p-values for testing the hypothesis of no partial
%                   correlation.  Choices are:
%                      TAIL         Alternative Hypothesis
%                   ---------------------------------------------------
%                     'both'     correlation is not zero (the default)
%                     'right'    correlation is greater than zero
%                     'left'     correlation is less than zero
%
%   PARTIALCORRI computes p-values for linear and rank partial correlations
%   using a Student's t distribution for a transformation of the correlation.
%   This is exact for linear partial correlation when X and Z are normal, but
%   is a large-sample approximation otherwise.
%
% Examples:
%
% The "carsmall" dataset contains performance and design measurements
% on cars manufactured in 1970, 1976, and 1982.  Take MPG and 
% Acceleration as performance measures, and Displacement, Horsepower,
% and Weight as design variables.  Displacement, Horsepower, and
% Weight all have predictable qualitative relationships with MPG
% and Acceleration.  But, as an example, it is less clear whether
% Displacement retains a strong association with Acceleration, once the
% effect of Horsepower has been taken into account.  The partial
% correlation coefficient is one way to quantify this conditional
% association.
%  
% load carsmall
% 
% The variables MPG and Acceleration contain some missing values.
% Unless we remove them, the outputs from corr and partialcorri will 
% be all NaNs.  Use the 'rows' parameter to exclude observations with
% missing values.
%
% [coef,pval] = corr([MPG,Acceleration], ...
%     [Displacement,Horsepower,Weight], 'rows','complete');
% coef = array2table(coef, ...
%     'VariableNames',{'Displacement','Horsepower','Weight'}, ...
%     'RowNames',{'MPG','Acceleration'});
% 
% It may be confusing that "Acceleration" has a negative correlation
% with "Weight" and "Horsepower".  However, the "Acceleration" variable
% represents the time required to accelerate from a standing stop to
% 60 miles per hour.  Therefore, a high value for "Acceleration"
% corresponds to a vehicle with low acceleration.
% 
% pval = array2table(pval, ...
%     'VariableNames',{'Displacement','Horsepower','Weight'}, ...
%     'RowNames',{'MPG','Acceleration'});
% 
% disp('Correlation Coefficients')
% disp(coef)
% disp('P-values')
% disp(pval)
% 
% [coefi,pvali] = partialcorri([MPG,Acceleration], ...
%     [Displacement,Horsepower,Weight], ...
%     'rows','complete');
% 
% coefi = array2table(coefi, ...
%     'VariableNames',{'Displacement','Horsepower','Weight'}, ...
%     'RowNames',{'MPG','Acceleration'});
% 
% pvali = array2table(pvali, ...
%     'VariableNames',{'Displacement','Horsepower','Weight'}, ...
%     'RowNames',{'MPG','Acceleration'});
% 
% disp('Partial Correlation Coefficients and P-values')
% disp(coefi)
% disp('P-values')
% disp(pvali)
% 
% To illustrate use of the three matrix form of input, imagine that
% you have an additional measurement, "Headwind", the average headwind
% along the route on which the performance measurements were made. (The 
% headwind variable is fictitious and it is contrived just for the
% purpose of illustration.)  You know that the strength of a headwind
% or tailwind will affect MPG and zero-to-60 times.  Potentailly, a headwind 
% variable might also interact with the design variables: for example,
% tests with high displacement engines might occur on days with 
% disproportionately strong headwinds.  You might regard "Headwind" as
% a nuisance variable: you don't care about headwind, per se, but you
% want to condition the other measurements with respect to any confounding
% effects that environmental headwinds may have had.
% 
% Headwind = (10:-0.2:-9.8)' + 5*randn(100,1);
% [coefih,pvalih] = partialcorri([MPG,Acceleration], ...
%     [Displacement,Horsepower,Weight], ...
%     Headwind, ...
%     'rows','complete');
% 
% coefih = array2table(coefih, ...
%     'VariableNames',{'Displacement','Horsepower','Weight'}, ...
%     'RowNames',{'MPG','Acceleration'});
% 
% pvalih = array2table(pvalih, ...
%     'VariableNames',{'Displacement','Horsepower','Weight'}, ...
%     'RowNames',{'MPG','Acceleration'});
% 
% disp('Partial Correlation Coefficients and P-values, Accounting for Headwind')
% disp(coefih)
% disp('P-values')
% disp(pvalih)
%
%   See also PARTIALCORR, CORR, CORRCOEF, TIEDRANK.

%   References:  
%   (1) Kendall's Advanced Theory of Statistics, Stuart, Ord and Arnold, 6th Ed, 
%       Volume 2A, 2nd Ed, 2004 (Chapter 28).
%   (2) The distribution of the partial correlation coefficient, 
%       Fisher, R.A., Metron, Vol. 3, p. 329-332.

%   Copyright 1984-2014 The MathWorks, Inc.

% twoStage is true if the command line had partialcorri(X,Y,Z,...).
% An empty value for Z is created if not.
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

[X,Y,tail,twoStage,Z,pairwise] = processCommandLine(Y,X,varargin{:});

wantPvals = (nargout > 1);

alreadyWarnedPerfectlyConditionedXY = false;

[~,P] = size(Y);
[~,Q] = size(X);

coef = NaN(P,Q);
if wantPvals
    pval = coef;
end

if any(any(isnan(X))) || any(any(isnan(Z)))
    % If 'rows' is 'all' (the default), NaNs have not been removed.
    % The computation for all return values would yield NaNs, 
    % so just return them now.
    return
end

if isempty(Y)
    % Row deletions for missing values may have left Y empty.
    % Since processCommandLine tests for equal row counts and then
    % removes all observations with missing in X and Z, if Y is
    % empty at this point, so are X and Z.
    warning(message('stats:partialcorr:EmptyAfterDeletions'));
    return
end

if ~(pairwise && any(any(isnan(Y))))
    [coef,pv,~] = ...
        computeAll(Y,X,Z,tail,twoStage,wantPvals,alreadyWarnedPerfectlyConditionedXY);
    if wantPvals
        pval = pv;
    end
else
    % pairwise row deletion and some NaNs in Y.
    % Treat each column of Y separately, 
    % removing rows with NaN specific to that column of Y.
    
    % We warn if missing value deletions remove all observations.
    % Flag so we only issue the warning once (the situation might occur
    % any column of Y).
    warnAllDelete = true;
    
    for i=1:P
        notnan = ~isnan(Y(:,i));
        if ~any(notnan)
            % Row deletions for missing values have left column i of Y
            % empty (also, the adjusted X and Z). We pre-filled coefficients 
            % and pvals with NaNs, so just go on to the next column.
            if warnAllDelete
                warning(message('stats:partialcorr:EmptyAfterDeletions'));
                warnAllDelete = false;
            end
            
        else
            
            % We have useful observations.
            % Protect against indexing an empty Z
            if twoStage
                [cf,pv,alreadyWarnedPerfectlyConditionedXY] = ...
                    computeAll(Y(notnan,i),X(notnan,:),Z(notnan,:),tail,...
                    twoStage,wantPvals,alreadyWarnedPerfectlyConditionedXY);
            else
                [cf,pv,alreadyWarnedPerfectlyConditionedXY] = ...
                    computeAll(Y(notnan,i),X(notnan,:),Z,tail,twoStage,...
                    wantPvals,alreadyWarnedPerfectlyConditionedXY);
            end
            coef(i,:) = cf;
            if wantPvals
                pval(i,:) = pv;
            end
        end
    end
end

end %-partialcorri

% SUBFUNCTIONS -------------------------------------------------------

% computeAll ---------------------------------------------------------

function [coef, pval, alreadyWarnedPerfectlyConditionedXY] = ...
    computeAll(Y,X,Z,tail,twoStage,wantPvals, alreadyWarnedPerfectlyConditionedXY)
X = bsxfun(@minus,X,mean(X));
Y = bsxfun(@minus,Y,mean(Y));
if twoStage
    Z = bsxfun(@minus,Z,mean(Z));
end
[coef,rankX,rankZ,allNaN,alreadyWarnedPerfectlyConditionedXY] = ...
    computeCoefficients(Y,X,Z,twoStage,wantPvals,alreadyWarnedPerfectlyConditionedXY);
if wantPvals
    if allNaN
        % coefs are all NaN because of missing values
        pval = coef;
    else
        pval = computePvals(coef,size(Y,1),rankX,rankZ,tail);
    end
else
    pval = [];
end
end %-computeAll

% computePvals -------------------------------------------------------

function pval = computePvals(coef,N,rankX,rankZ,tail)
df = max(N-1-rankX-rankZ,0);
t = sign(coef) .* Inf;
k = (abs(coef) < 1);
t(k) = coef(k) ./ sqrt(1-coef(k).^2);
t = sqrt(df).*t;
switch tail
    case 'b' % 'both or 'ne'
        pval = 2*tcdf(-abs(t),df);
    case 'r' % 'right' or 'gt'
        pval = tcdf(-t,df);
    case 'l' % 'left or 'lt'
        pval = tcdf(t,df);
end
end %-computePvals

% computeCoefficients ------------------------------------------------

function [coef,rankX,rankZ,allNaN,alreadyWarnedPerfectlyConditionedXY] = ...
    computeCoefficients(Y,X,Z,twoStage,wantPvals,alreadyWarnedPerfectlyConditionedXY)
allNaN = false;
rankZ = 0;
if twoStage
    if wantPvals
        rankZ = rank(Z);
    else
        % We won't need the actual value of rank(Z), just an approximation
        % sufficient for calculating tolerances for equating residuals
        % to zero.
        rankZ = size(Z,2);
    end
    % Regress Y and X on Z and replace them with the residuals.
    Y = flooredResiduals(Y,Z,rankZ);
    X = flooredResiduals(X,Z,rankZ);
    Y = bsxfun(@minus,Y,mean(Y));
    X = bsxfun(@minus,X,mean(X));
    if any(all([Y X]==0)) & ~alreadyWarnedPerfectlyConditionedXY
        % Some columns of Y and/or X are perfectly explained by Z.
        warning(message('stats:partialcorr:PerfectlyConditionedXY'));
        alreadyWarnedPerfectlyConditionedXY = true;
    end
end

rankX = rank(X);

[N,P] = size(Y);
Q = size(X,2);

if rankX == Q
    % Full column rank.  We can use the quick QR  / post-permute method.
    coef = orthopcorr(Y,X);
else
    % Perform separate regressions for each Y_i vs X_j pair.
    coef = zeros(P,Q);
    YX = [Y,zeros(N,1)];
    for j=1:Q
        Xnotj = X(:,~ismember(1:Q,j));
        YX(:,P+1) = X(:,j);
        resid = flooredResiduals(YX,Xnotj,rankX);
        coef(:,j) = correlateResiduals(resid);
    end
end
end %-computeCoefficients

% orthopcorr   -------------------------------------------------------

function rho = orthopcorr(Y,X)
%RHO = ORTHOPCORR(Y,X) Partial correlation, Y with X adjusted for X.
%
% Let C(X) denote the column space of X.
% Let Xi denote column i of X, and let X[i] denote X absent column i.
% For arbitrary Y, let P(X)Y denote the projection of Y onto C(X).
%
% Suppose X is n x p, n>p, X is full column rank, and
% [Q,R] = qr(X,0) is the economy QR decomposition of X.
% - Qp is orthogonal to X[p], thus colinear with the residuals from
%   regressing Xp on X[p].
% - For any vector Y and column i, P(Q[i]) == P(Q)Y - P(Qi)Y;
% - Therefore the residual of Y controlling for X[p] is:
%   r = Y - P(Q)Y + P(Qp)Y;
% - And the Pearson partial corr of Xp and Y controlling for X[p] is:
%   <r,Qp> / norm(r);
% - Subsequently,  utilize Q*R*U = A*U, where U is a permutation
%   matrix that places other columns of X in the rightmost position.

% Throughout, we will set very small residuals to zero.  We do this to
% avoid spurious partial correlation coefficients when residuals should
% be zero.
sz = size(X);
tol = residualTolerance(sz(1),sz(2)-1,Y);

[Q,R] = qr(X,0);
Qy = Q'*Y;

% Residuals from the projection of Y on all of X (equivalently, all of Q).
res = Y - Q*Qy;
rss = sum(res.^2);
rss(rss<tol) = 0;

P1 = size(Y,2);
P2 = size(X,2);

rho = zeros(P1,P2);

qi = Qy(P2,:);
qi(abs(qi)<tol) = 0;

rho(:,P2) = sign(R(end)) * qi ./ sqrt(rss + qi.^2);

for i = 1:(P2-1)
    U = eye(P2);
    U(:,[i P2]) = U(:,[P2 i]);
    [qq,rr] = qr(R*U);
    qi = qq(:,P2)' * Qy;
    qi(abs(qi)<tol) = 0;
    rho(:,i) = sign(rr(end)) * qi ./ sqrt(rss + qi.^2);
end
end %-orthopcorr

% processCommandLine -------------------------------------------------

function [X,Y,tail,twoStage,Z,pairwise] = processCommandLine(Y,X,varargin)

if ~ismatrix(Y) || ~isnumeric(Y) || ~ismatrix(X) || ~isnumeric(X)
    error(message('stats:partialcorr:InputsMustBeMatrices'));
end

[N,~] = size(Y);
[Nx,~] = size(X);

if Nx ~= N
    error(message('stats:partialcorr:InputSizeMismatch'));
end;

% "twoStage" indicates
% (1) There is a 3rd matrix input (aka "Z").
% (2) We will regress Y & X on Z to get residuals Yres and Xres
% (3) Then get partial correlations of Yres(:,i) and Xres(:,j),
%     conditioned on Xres(:,notj).
twoStage = false;
Z = [];

if ~isempty(varargin) && isnumeric(varargin{1})
    twoStage = true;
    Z = varargin{1};
    varargin(1) = [];
    if ~ismatrix(Z) || ~isnumeric(Z)
        error(message('stats:partialcorr:InputsMustBeMatrices'));
    end
    [Nz,~] = size(Z);
    if Nz ~= N
        error(message('stats:partialcorr:InputSizeMismatch'));
    end;
end

if isempty(Y)
    % Y and X (optionally, Z) have been checked for same number of rows.
    % If Y is empty, all of the matrix inputs are empty.
    error(message('stats:partialcorr:EmptyInput'));
end

pnames = {'type'  'rows' 'tail'};
dflts  = {'p'     'a'    'both'};
[type,rows,tail] = internal.stats.parseArgs(pnames,dflts,varargin{:});

% Validate the rows parameter.
rows = internal.stats.getParamVal(rows,{'all' 'complete' 'pairwise'},'''Rows''');
rows = rows(1);

% Validate the type parameter.
try
    type = internal.stats.getParamVal(type,{'pearson' 'kendall' 'spearman'},'''Type''');
catch
    error(message('stats:partialcorr:UnknownType'));
end
type = type(1);
if type == 'k'
    error(message('stats:partialcorr:Kendall'));
end

% Validate the tail parameter.
tailChoices = {'left','both','right'};
if ischar(tail) && (size(tail,1)==1)
    i = find(strncmpi(tail,tailChoices,length(tail)));
    if isempty(i)
        i = find(strncmpi(tail,{'lt','ne','gt'},length(tail)));
    end
    if isscalar(i)
        tail = tailChoices{i}(1);
    elseif isempty(i)
        error(message('stats:partialcorr:UnknownTail'));
    end
else
    error(message('stats:partialcorr:UnknownTail'));
end

if rows ~= 'a'
    % In all syntax combinations, each column of the matrix X will
    % be regressed on the other columns of X.  If any row of X has
    % a missing value, the corresponding residual will be undefined,
    % whether X(:,j) is NaN, or one of the "predictor" columns is NaN.
    % Therefore, X must be pruned of rows containing NaN values.  Residuals
    % from the matrix Y will be paired with those from X.  Therefore, rows
    % containing NaNs in X must be removed from Y (and in the "twoStage"
    % case, rows of Z containing NaN as well).
    missing = any(isnan(X),2);
    if twoStage
        missing = missing | any(isnan(Z),2);
        Z = Z(~missing,:);
    end
    X = X(~missing,:);
    Y = Y(~missing,:);
    if rows == 'c'
        % If there are additional rows with missing values in Y,
        % and if 'complete observations were requested, we delete
        % those rows across the board.
        missing = any(isnan(Y),2);
        Y = Y(~missing,:);
        X = X(~missing,:);
        if twoStage
            Z = Z(~missing,:);
        end
    end
end

pairwise = (rows == 'p');

if type == 's'
    X = tiedrank(X);
    % tiedrank handles NaNs in Y gracefully (possible if 'pairwise')
    Y = tiedrank(Y);
end

if twoStage
    if type == 's'
        Z = tiedrank(Z);
    end
end

end %-processCommandLine

% residualTolerance --------------------------------------------------

function tol = residualTolerance(N,rankX,Y)
tol = max(N,rankX) .* eps(class(Y)) .* sqrt(sum(Y.^2,1));
end

% flooredResiduals ---------------------------------------------------

function residuals = flooredResiduals(Y,X,rankX)
% Some Y(:,i) might be perfectly predictable from X,
% and the residuals should then be zero, but roundoff could throw
% that off slightly.  If a column of residuals is effectively zero
% relative to the original variable, then assume we've predicted
% exactly.  This prevents computing spuriously valid correlations
% when they really should be NaN.

% Turn off rank deficiency warning from backslash, since a basic solution is
% perfectly fine for calculation of residuals.
savedWarnState = warning('off','MATLAB:rankDeficientMatrix');
cleanupObj = onCleanup(@() warning(savedWarnState));

N = size(Y,1);
tol = residualTolerance(N,rankX,Y);
% Use linsolve instead of "residuals = Y - X * (X\Y)" because '\' can
% give NaNs or +/-Inf when X is square.
opt.RECT = true;
residuals = Y - X * linsolve(X,Y,opt);
residuals(:, sqrt(sum(abs(residuals).^2,1)) < tol ) = 0;
end

% correlateResiduals -------------------------------------------------

function coef = correlateResiduals(resid)
Xresid = resid(:,end);  % These are the X column residuals
resid(:,end) = [];      % These are the Y column(s) residuals

coef = resid' * Xresid / norm(Xresid);
coef = coef ./ sqrt(sum(resid.^2,1))';
end
