function [coef, pval] = corr(x, varargin)
%CORR Linear or rank correlation.
%   RHO = CORR(X) returns a P-by-P matrix containing the pairwise linear
%   correlation coefficient between each pair of columns in the N-by-P
%   matrix X.
%
%   RHO = CORR(X,Y,...) returns a P1-by-P2 matrix containing the pairwise
%   correlation coefficient between each pair of columns in the N-by-P1 and
%   N-by-P2 matrices X and Y.
%
%   [RHO,PVAL] = CORR(...) also returns PVAL, a matrix of p-values for
%   testing the hypothesis of no correlation against the alternative that
%   there is a non-zero correlation.  Each element of PVAL is the p-value
%   for the corresponding element of RHO.  If PVAL(i,j) is small, say less
%   than 0.05, then the correlation RHO(i,j) is significantly different
%   from zero.
%
%   [...] = CORR(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values.  Valid parameters are the following:
%
%        Parameter  Value
%         'type'    'Pearson' (the default) to compute Pearson's linear
%                   correlation coefficient, 'Kendall' to compute Kendall's
%                   tau, or 'Spearman' to compute Spearman's rho.
%         'rows'    'all' (default) to use all rows regardless of missing
%                   values (NaNs), 'complete' to use only rows with no
%                   missing values, or 'pairwise' to compute RHO(i,j) using
%                   rows with no missing values in column i or j.
%         'tail'    The alternative hypothesis against which to compute
%                   p-values for testing the hypothesis of no correlation.
%                   Choices are:
%                      TAIL         Alternative Hypothesis
%                   ---------------------------------------------------
%                     'both'     correlation is not zero (the default)
%                     'right'    correlation is greater than zero
%                     'left'     correlation is less than zero
%
%   The 'pairwise' option for the 'rows' parameter can produce RHO that is
%   not positive semi-definite.  The 'complete' option always produces a
%   positive semi-definite RHO, but when data are missing, the estimates
%   may be based on fewer observations.
%
%   CORR computes p-values for Pearson's correlation using a Student's t
%   distribution for a transformation of the correlation.  This is exact
%   when X and Y are normal.  CORR computes p-values for Kendall's tau and
%   Spearman's rho using either the exact permutation distributions (for
%   small sample sizes), or large-sample approximations.
%
%   CORR computes p-values for the two-tailed test by doubling the more
%   significant of the two one-tailed p-values.
%
%   See also CORRCOEF, PARTIALCORR, TIEDRANK.

%   References:
%      [1] Gibbons, J.D. (1985) Nonparametric Statistical Inference,
%          2nd ed., M. Dekker.
%      [2] Hollander, M. and D.A. Wolfe (1973) Nonparametric Statistical
%          Methods, Wiley.
%      [3] Kendall, M.G. (1970) Rank Correlation Methods, Griffin.
%      [4] Best, D.J. and D.E. Roberts (1975) "Algorithm AS 89: The Upper
%          Tail Probabilities of Spearman's rho", Applied Statistics,
%          24:377-379.

%   Copyright 1984-2019 The MathWorks, Inc.


%   Spearman's rho is equivalent to the linear (Pearson's) correlation
%   between ranks.  Kendall's tau is equivalent to the linear correlation
%   between the "concordances" sign(x(i)-x(j))*sign(y(i)-y(j)), i<j, with
%   an adjustment for ties.  This is often referred to as tau-b.
%
%   Kendall's tau-b is identical to the standard tau (or tau-a) when there
%   are no ties.  However, tau-b includes an adjustment for ties in the
%   normalizing constant.

%   Spearman's rho and Kendall's tau are discrete-valued statistics, and
%   their distributions have positive probability at 1 and -1.  For small
%   sample sizes, CORR uses the exact permutation distributions, and thus,
%   the on-diagonal p-values from CORR(X,X) in those cases.

%   When there are ties in the data, the null distribution of Spearman's
%   rho and Kendall's tau may not be symmetric.  Computing a two-tailed
%   p-value in such cases is not well-defined.  CORR computes p-values for
%   the two-tailed test by doubling the smaller of the one-tailed p-values.

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 1 || isempty(x)
    error(message('stats:corr:TooFewInputs'));
end
[n,p1] = size(x);

% Only x given, compute pairwise rank correlations
if (nargin < 2) || ischar(varargin{1})
    corrXX = true;
    y = x;
    p2 = p1;

% Both x and y given, compute the pairwise rank cross correlations
else
    y = varargin{1};
    varargin = varargin(2:end);
    if size(y,1) ~= n
        error(message('stats:corr:InputSizeMismatch'));
    end
    corrXX = false;
    p2 = size(y,2);
end

pnames = {'type'  'rows' 'tail'};
dflts  = {'p'     'a'    'both'};
[type,rows,tail] = internal.stats.parseArgs(pnames,dflts,varargin{:});

% Validate the rows parameter.
rowsChoices = {'all' 'complete' 'pairwise'};
if ischar(rows)
    i = find(strncmpi(rows,rowsChoices,length(rows)));
    if isscalar(i)
        rows = rowsChoices{i}(1);
    end
end
switch rows
case 'a' % 'all'
    % Missing values are left in, so all cols have the same length.
case 'c' % 'complete'
    % Complete rows with missing values are removed, so all cols have the
    % same length.
    ok = ~any(isnan(x),2);
    if ~corrXX
        ok = ok & ~any(isnan(y),2);
    end
    n = sum(ok);
    if n < size(x,1)
        x = x(ok,:);
        if ~corrXX
            y = y(ok,:);
        end
    end
case 'p' % 'pairwise'
    % Missing values are removed pairwise, so each column pair may have a
    % different length, we'll have to see.
otherwise
    error(message('stats:corr:UnknownRows'));
end

% Can't do anything without at least two observations.
if n < 2
    outType = internal.stats.dominantType(x,y);
    coef = NaN(p1,p2,"like",outType);
    pval = NaN(p1,p2,"like",outType);
    return
end

% Validate the type parameter, and set a handle to the correlation function.
typeChoices = {'pearson' 'kendall' 'spearman'};
if ischar(type)
    i = find(strncmpi(type,typeChoices,length(type)));
    if isscalar(i)
        type = typeChoices{i}(1);
    end
end
switch type
case 'p' % Pearson linear correlation
    corrFun = @internal.stats.corrPearson;
case 'k' % Kendall rank correlation
    corrFun = @internal.stats.corrKendall;
case 's' % Spearman rank correlation
    corrFun = @internal.stats.corrSpearman;
otherwise
    error(message('stats:corr:UnknownType'));
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
        error(message('stats:corr:UnknownTail'));
    end
else
    error(message('stats:corr:UnknownTail'));
end

% Some things we just can't do for complex data.
complexdata = ~(isreal(x) || all(imag(x(:))==0));
if ~complexdata && ~corrXX
    complexdata = ~(isreal(y) || all(imag(y(:))==0));
end
if complexdata
    if nargout > 1
        error(message('stats:corr:PValComplexData'));
    elseif type ~= 'p'
        error(message('stats:corr:RankCorComplexData'));
    end
end

if corrXX
    if nargout<2
        coef = corrFun(rows,tail,x);
    else
        [coef,pval] = corrFun(rows,tail,x);
    end
else
    if nargout<2
        coef = corrFun(rows,tail,x,y);
    else
        [coef,pval] = corrFun(rows,tail,x,y);
    end
end
