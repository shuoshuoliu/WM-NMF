function [coef,pval] = corrSpearman(rows,tail,x,y)
%CORRPEARSON Computation of Spearman correlation matrix
%
%   [COEF,PVAL] = CORRSPEARMAN(ROWS,TAIL,X,Y) calculates the Spearman
%   cross-correlation matrix COEF. ROWS must be 'all', 'complete', or
%   'pairwise'. TAIL must be 'both', 'left', or 'right'.
%
%   See also: CORR, INTERNAL.STATS.CORRKENDALL, INTERNAL.STATS.CORRPEARSON.

% Copyright 2019-2020 The MathWorks, inc.
[n,p1] = size(x);
corrXX = (nargin < 4);
if corrXX
    p2 = p1;
    y = x;
else
    p2 = size(y,2);
end
outType = internal.stats.dominantType(x,y);

% Do column-pairwise NaN removal only if there are missing values.
if rows == 'p'
    isnanX = isnan(x);
    pairwise = any(isnanX(:));
    if corrXX
        isnanY = isnanX;
    else
        isnanY = isnan(y);
        pairwise = pairwise || any(isnanY(:));
    end
else
    pairwise = false;
end

% Computation of rank correlation, with column-pairwise NaN removal.  Can't
% hand off to corrPearson as for 'all' or 'complete', because we need to
% compute the ranks separately for each pair of columns.
if pairwise
    % Preallocate the outputs.
    coef = zeros(p1,p2,"like",outType);
    if nargout > 1
        pval = zeros(p1,p2,"like",outType);
        Dstat = zeros(p1,p2,"like",outType); % save the obs. stat. for p-value computation
        if corrXX
            needPVal = tril(true(p1),-1); % lower triangle only, no diagonal
        else
            needPVal = true(p1,p2);
        end
    end
    
    for i = 1:p1
        % Loop over all columns of y for corr(x,y), or over the lower triangle
        % and the diagonal for corr(x).
        j1 = p2*(1-corrXX) + i*corrXX; % 1:p2, or 1:i
        for j = 1:j1
            ok = ~(isnanX(:,i) | isnanY(:,j));
            nij = sum(ok);
            anyRowsRemoved = (nij < n);
            if anyRowsRemoved
                [xranki, xadj] = tiedrank(x(ok,i),0);
                [yrankj, yadj] = tiedrank(y(ok,j),0);
            else
                [xranki, xadj] = tiedrank(x(:,i),0);
                [yrankj, yadj] = tiedrank(y(:,j),0);
            end
            n3const = (nij+1)*nij*(nij-1) ./ 3;
            ties = ((xadj>0) || (yadj>0));
            
            D = sum((xranki - yrankj).^2);
            meanD = (n3const - (xadj+yadj)./3) ./ 2;
            
            % When either of the data vectors is constant, stdD should be
            % zero, get that exactly.
            n3const2 = (nij+1)*nij*(nij-1)/2;
            if (xadj == n3const2) || (yadj == n3const2)
                stdD = 0;
            else
                stdD = sqrt((n3const./2 - xadj./3).*(n3const./2 - yadj./3)./(nij-1));
            end
            
            coef(i,j) = (meanD - D) ./ (sqrt(nij-1)*stdD);
            if (i == j) && corrXX
                % Put an exact one on the diagonal for autocorrelation.
                coef(i,i) = sign(coef(i,i)); % preserves NaNs
                
                % Compute on-diag p-values for autocorrelation later
                
            elseif nargout > 1
                % If there are ties, or if there has been pairwise removal
                % of missing data, we'll compute the p-value separately here.
                if (anyRowsRemoved || ties)
                    pval(i,j) = pvalSpearman(tail, D, meanD, stdD, xranki, yrankj);
                    needPVal(i,j) = false; % this one's done
                else
                    Dstat(i,j) = D;
                end
            end
        end
    end
    % Limit off-diag correlations to [-1,1].
    t = find(abs(coef) > 1); coef(t) = coef(t)./abs(coef(t)); % preserves NaNs
    
    % Calculate the remaining p-values, except not the on-diag elements for
    % autocorrelation.  All cases with no ties and no removed missing values
    % can be computed based on a single null distribution.
    if nargout > 1 && any(needPVal(:))
        n3const = (n+1)*n*(n-1) ./ 3;
        meanD = n3const./2;
        stdD = n3const ./ (2*sqrt(n-1));
        pval(needPVal) = pvalSpearman(tail,Dstat(needPVal),meanD,stdD,n);
    end
    
    % If this is autocorrelation, reflect the lower triangle into the upper.
    if corrXX
        coef = tril(coef) + tril(coef,-1)'; % leave the diagonal alone
        if nargout > 1
            pval = pval + pval';
            % The p-values along the diagonal are always exactly one (unless
            % they're NaN), because Pr{coef(i,i) as/more extreme than 1} == 1,
            % regardless of which tail(s) we're testing against.
            pval(1:p1+1:end) = sign(diag(coef)); % preserves NaNs on diag
        end
    end
    
    % Vectorized computation of rank correlation.  No NaN removal, or NaNs already
    % removed in complete rows
else
    [xrank, xadj] = tiedrank(x,0);
    ties = any(xadj>0);
    if corrXX
        yrank = xrank; yadj = xadj;
    else
        [yrank, yadj] = tiedrank(y,0);
        ties = ties || any(yadj>0);
    end
    n3const = (n+1)*n*(n-1) ./ 3;
    
    % Compute all elements at once, fastest when there are no ties, or
    % if we don't need p-values
    if nargout == 1 || ~ties
        if corrXX
            coef = internal.stats.corrPearson(rows,tail,xrank);
        else
            coef = internal.stats.corrPearson(rows,tail,xrank,yrank);
        end
        if nargout > 1 % p-values requested, but no ties
            meanD = n3const./2;
            stdD = n3const ./ (2*sqrt(n-1));
            D = meanD - round(coef*(sqrt(n-1)*stdD));
            if corrXX
                % Compute p-values for the lower triangle and reflect into the upper.
                pval = zeros(p1,"like",outType);
                ltri = (tril(ones(size(pval)),-1) > 0); % lower triangle, no diagonal
                pval(ltri) = pvalSpearman(tail,D(ltri),meanD,stdD,n);
                pval = pval + pval';
                % The p-values along the diagonal are always exactly one
                % (unless they're NaN), because Pr{coef(i,i) as/more extreme
                % than 1} == 1, regardless of which tail(s) we're testing
                % against.
                pval(1:p1+1:end) = sign(diag(coef)); % preserves NaNs on diag
            else
                pval = pvalSpearman(tail,D,meanD,stdD,n);
            end
        end
        
        % Compute one row at a time when we need p-values and there are ties
    else
        coef = zeros(p1,p2,"like",outType);
        pval = zeros(p1,p2,"like",outType);
        for i = 1:p1
            xranki = xrank(:,i);
            if corrXX
                j = 1:i; % lower triangle and diagonal
                yrankj = xrank(:,j);
                yadjj = yadj(j);
            else
                j = 1:p2;
                yrankj = yrank;
                yadjj = yadj;
            end
            D = sum((xranki-yrankj).^2); % sum((xranki - yrankj).^2);
            meanD = (n3const - (xadj(i)+yadjj)./3) ./ 2;
            stdD = sqrt((n3const./2 - xadj(i)./3)*(n3const./2 - yadjj./3)./(n-1));
            
            % When either of the data vectors is constant, stdD should be
            % zero, get that exactly.
            n3const2 = (n+1)*n*(n-1)/2;
            stdD((xadj(i) == n3const2) | (yadjj == n3const2)) = 0;
            
            coef(i,j) = (meanD - D) ./ (sqrt(n-1)*stdD);
            pval(i,j) = pvalSpearman(tail,D,meanD,stdD,xranki,yrankj);
        end
        
        % Limit off-diag correlations to [-1,1].
        t = find(abs(coef) > 1); coef(t) = coef(t)./abs(coef(t)); % preserves NaNs
        
        % If this is autocorrelation, reflect the lower triangle into the upper.
        if corrXX
            coef = tril(coef) + tril(coef,-1)'; % leave diagonal alone
            % The p-values along the diagonal are always exactly one (unless
            % they're NaN), because Pr{coef(i,i) as/more extreme than 1} == 1,
            % regardless of which tail(s) we're testing against.
            pval = pval + pval';
            pval(1:p1+1:end) = sign(diag(coef)); % preserve NaNs on diag
        end
    end
end


%--------------------------------------------------------------------------

function p = pvalSpearman(tail, D, meanD, stdD, arg1, arg2)
%PVALSPEARMAN Tail probability for Spearman's D statistic.

% Without ties, D is symmetric about (n^3-n)/6, taking on values
% 0:2:(n^3-n)/3.  With ties, it's still in (but not on) that range, but
% not symmetric, and can take on odd and half-integer values.

% D, meanD, and stdD may be vectors when n is given (i.e., when there were no
% ties). When xrank and yrank are given (i.e., when there were ties), yrank
% may be a matrix, but xrank must be a single column, and D, meanD, stdD must
% have the same lengths as yrank has columns.

if nargin < 6 % pvalSpearman(tail, D, meanD, stdD, n), no ties
    noties = true;
    n = arg1;
else % pvalSpearman(tail, D, meanD, stdD, xrank, yrank), ties in data
    noties = false;
    
    % If stdD is zero, at least one of the data vectors was constant.  The
    % correlation coef in these cases falls out of the calculations correctly
    % as NaN, but the exact p-value calculations below would regard Pr{D==D0}
    % as 1.  Return a NaN p-value instead.
    D(stdD == 0) = NaN;
    
    xrank = arg1;
    yrank = arg2;
    n = length(xrank);
end
exact = (n < 10);
if exact
    p = NaN(size(D),"like",D);
    if noties
        tailProb = spearmanExactSub(tail,n);
        t = ~isnan(D(:));
        p(t) = tailProb(2*D(t)+1); % bins at half integers, starting at zero
    else
        for j = 1:size(yrank,2)
            if ~isnan(D(j))
                tailProb = spearmanExactSub(tail,xrank,yrank(:,j));
                p(j) = tailProb(2*D(j)+1); % bins at half integers, starting at zero
            end
        end
    end
    
else
    if noties
        % Use AS89, an Edgeworth expansion for upper tail prob of D.
        n3const = (n^3 - n)/3;
        switch tail
            case 'b' % 'both or 'ne'
                p = AS89(max(D, n3const-D), n, n3const);
                p = 2*p; p(p>1) = 1; % Don't count continuity correction at center twice
            case 'r' % 'right' or 'gt'
                p = AS89(n3const - D, n, n3const);
            case 'l' % 'left' or 'lt'
                p = AS89(D, n, n3const);
        end
    else
        % Use a t approximation.
        r = (meanD - D) ./ (sqrt(n-1)*stdD);
        t = Inf*sign(r);
        ok = (abs(r) < 1);
        t(ok) = r(ok) .* sqrt((n-2)./(1-r(ok).^2));
        switch tail
            case 'b' % 'both or 'ne'
                p = 2*tcdf(-abs(t),n-2);
            case 'r' % 'right' or 'gt'
                p = tcdf(-t,n-2);
            case 'l' % 'left' or 'lt'
                p = tcdf(t,n-2);
        end
    end
end


%--------------------------------------------------------------------------

function tailProb = spearmanExactSub(tail,arg1,arg2)
if nargin < 3
    % No ties, take permutations of 1:n
    n = arg1;
    nfact = factorial(n);
    Dperm = sum((repmat(1:n,nfact,1) - perms(1:n)).^2, 2);
else
    % Ties, take permutations of the midranks.
    %
    % With ties, we could consider only distinguishable permutations
    % (those for which equal ranks (ties) are not simply interchanged),
    % but generating only those is a bit of work.  Generating all
    % permutations uses more memory, but gives the same result.
    xrank = arg1;
    yrank = arg2;
    n = length(xrank);
    nfact = factorial(n);
    Dperm = sum((repmat(xrank(:)',nfact,1) - perms(yrank(:)')).^2, 2);
end
n3const = (n^3 - n)/3;
freq = histcounts(Dperm,(-.25):.5:(n3const+.25));

% Get the tail probabilities.  Reflect as necessary to get the correct
% tail: the left tail of D corresponds to right tail of rho.
switch tail
case 'b' % 'both or 'ne'
    % Use twice the smaller of the tail area above and below the
    % observed value.
    tailProb = min(cumsum(freq), rcumsum(freq,nfact)) ./ nfact;
    tailProb = min(2*tailProb, 1); % don't count the center bin twice
case 'r' % 'right' or 'gt'
    tailProb = cumsum(freq) ./ nfact;
case 'l' % 'left' or 'lt'
    tailProb = rcumsum(freq,nfact) ./ nfact;
end


%--------------------------------------------------------------------------

function p = AS89(D, n, n3const)
%AS89 Upper tail probability for Spearman's D statistic.
%   Edgeworth expansion for the upper tail probability of D in the no ties
%   case, with continuity correction, adapted from Applied Statistics
%   algorithm AS89.
c = [.2274 .2531 .1745 .0758 .1033 .3932 .0879 .0151 .0072 .0831 .0131 .00046];
x = (2*(D-1)./n3const - 1) * sqrt(n - 1);
y = x .* x;
u = x .* (c(1) + (c(2) + c(3)/n)/n + ...
    y .* (-c(4) + (c(5) + c(6)/n)/n - ...
    y .* (c(7) + c(8)/n - y .* (c(9) - c(10)/n + ...
    y .* (c(11) - c(12) * y)/n))/n))/n;
p = u ./ exp(.5*y) + 0.5 * erfc(x./sqrt(2));
p(p>1) = 1; p(p<0) = 0; % don't ignore NaNs

%--------------------------------------------------------------------------

function y = rcumsum(x,sumx)
%RCUMSUM Cumulative sum in reverse direction. Note that this is the same as
% CUMSUM(X,'reverse') if SUMX==SUM(X), but can be used for cases where we
% want to count down from some other total too.
y = repmat(sumx,size(x));
y(2:end) = sumx - cumsum(x(1:end-1));