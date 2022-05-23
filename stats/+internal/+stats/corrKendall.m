function [coef,pval] = corrKendall(rows,tail,x,y)
%CORRKENDALL Computation of Kendall correlation matrix.
%
%   [COEF,PVAL] = CORRKENDALL(ROWS,TAIL,X,Y) calculates the Kendall
%   cross-correlation matrix COEF. ROWS must be 'all', 'complete', or
%   'pairwise'. TAIL must be 'both', 'left', or 'right'.
%
%   See also: CORR, INTERNAL.STATS.CORRPEARSON, INTERNAL.STATS.CORRSPEARMAN.

% Copyright 2019 The MathWorks, inc.

[n,p1] = size(x);
corrXX = (nargin < 4);
if corrXX
    p2 = p1;
    y = x;
else
    p2 = size(y,2);
end

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

% Compute the ranks once when not doing column-pairwise NaN removal.
if ~pairwise
    if rows == 'a'
        % For 'all', it's quicker to recognize that a column has a NaN in it
        % than to do the computation explicitly and end up with a NaN anyway.
        nancolX = any(isnan(x),1)';
        if corrXX
            nancol = nancolX | nancolX';
        else
            nancolY = any(isnan(y),1);
            nancol = nancolX | nancolY;
        end
        
    else % rows == 'c' | rows  == 'p'
        % For 'complete' or 'pairwise', there are no NaNs at this point
        nancol = false(p1,p2);
    end
    
    [xrank, xadj] = tiedrank(x,1);
    if corrXX
        yrank = xrank; yadj = xadj;
    else
        [yrank, yadj] = tiedrank(y,1);
    end
end

% Preallocate the outputs.
outClass = superiorfloat(x,y);
coef = zeros(p1,p2,outClass);
if nargout > 1
    pval = zeros(p1,p2,outClass);
    Kstat = zeros(p1,p2,outClass); % save the obs. stat. for p-value computation
    if corrXX
        needPVal = tril(true(p1),-1); % lower triangle only, no diagonal
    else
        needPVal = true(p1,p2);
    end
end

for i = 1:p1
    if ~pairwise
        xranki = xrank(:,i); xadji = xadj(:,i);
    end
    
    % Loop over all columns of y for corr(x,y), or over the lower triangle and
    % the diagonal for corr(x).
    j1 = p2*(1-corrXX) + i*corrXX; % 1:p2, or 1:i
    for j = 1:j1
        if pairwise % column-pairwise NaN removal
            ok = ~(isnanX(:,i) | isnanY(:,j));
            nij = sum(ok);
            anyRowsRemoved = (nij < n);
            if anyRowsRemoved
                [xranki, xadji] = tiedrank(x(ok,i),1);
                [yrankj, yadjj] = tiedrank(y(ok,j),1);
            else
                [xranki, xadji] = tiedrank(x(:,i),1);
                [yrankj, yadjj] = tiedrank(y(:,j),1);
            end
        else % no NaN removal, or NaNs already removed in complete rows
            % Quicker to check and bail out rather than compute the NaN below
            if nancol(i,j)
                coef(i,j) = NaN;
                if nargout > 1
                    pval(i,j) = NaN;
                    needPVal(i,j) = false;
                end
                continue
            end
            nij = n;
            yrankj = yrank(:,j); yadjj = yadj(:,j);
            anyRowsRemoved = false;
        end
        n2const = nij*(nij-1) / 2;
        ties = ((xadji(1)>0) || (yadjj(1)>0));
        
        K = 0;
        for k = 1:nij-1
            K = K + sum(sign(xranki(k)-xranki(k+1:nij)).*sign(yrankj(k)-yrankj(k+1:nij)));
        end
        coef(i,j) = K ./ sqrt((n2const - xadji(1)).*(n2const - yadjj(1)));
        
        % Clean up the diagonal for autocorrelation
        if (i == j) && corrXX
            if ~ties
                % Put an exact one on the diagonal only when there's no ties.
                % Kendall's tau may not be exactly one along the diagonal when
                % there are ties.
                coef(i,i) = sign(coef(i,i)); % preserves NaN
            end
            
            % Compute on-diag p-values for autocorrelation later
            
        elseif nargout > 1
            % If there are ties, or if there has been pairwise removal of
            % missing data, compute the p-value separately here.
            if ties
                % The tied case is sufficiently slower that we don't
                % want to do it if not necessary.
                %
                % When either of the data vectors is constant, stdK should be
                % zero, get that exactly.
                if (xadji(1) == n2const) || (yadjj(1) == n2const)
                    stdK = 0;
                else
                    stdK = sqrt(n2const*(2*nij+5)./9 ...
                        + xadji(1)*yadjj(1)./n2const ...
                        + xadji(2)*yadjj(2)./(18*n2const*(nij-2)) ...
                        - (xadji(3) + yadjj(3))./18);
                end
                pval(i,j) = pvalKendall(tail, K, stdK, xranki, yrankj);
                needPVal(i,j) = false; % this one's done
            elseif anyRowsRemoved
                stdK = sqrt(n2const*(2*nij+5)./9);
                pval(i,j) = pvalKendall(tail, K, stdK, nij);
                needPVal(i,j) = false; % this one's done
            else
                Kstat(i,j) = K;
            end
        end
    end
end
% Limit off-diag correlations to [-1,1].
t = find(abs(coef) > 1); coef(t) = coef(t)./abs(coef(t)); % preserves NaNs

% Calculate the remaining p-values, except not the on-diag elements for
% autocorrelation.  All cases with no ties and no removed missing values can
% be computed based on a single null distribution.
if nargout > 1 && any(needPVal(:))
    n2const = n*(n-1) ./ 2;
    stdK = sqrt(n2const * (2*n+5) ./ 9);
    pval(needPVal) = pvalKendall(tail,Kstat(needPVal),stdK,n);
end

% If this is autocorrelation, reflect the lower triangle into the upper.
if corrXX
    coef = tril(coef) + tril(coef,-1)'; % leave the diagonal alone
    if nargout > 1
        pval = pval + pval';
        % The p-values along the diagonal are always exactly one (unless
        % they're NaN) conditional on the pattern of ties, because
        % Pr{coef(i,i) as/more extreme than observed value} == 1, regardless
        % of which tail(s) we're testing against.
        pval(1:p1+1:end) = sign(diag(coef)); % preserves NaNs on diag
    end
end


%--------------------------------------------------------------------------

function p = pvalKendall(tail, K, stdK, arg1, arg2)
%PVALKENDALL Tail probability for Kendall's K statistic.

% Without ties, K is symmetric about zero, taking on values in
% -n(n-1)/2:2:n(n-1)/2.  With ties, it's still in that range, but not
% symmetric, and can take on adjacent integer values.

% K and stdK may be vectors when n is given (i.e., when there were no ties).
% When xrank and yrank are given (i.e., when there were ties), K and stdK must
% be scalars and xrank and yrank must both be a single column.

if nargin < 5 % pvalKendall(tail, K, stdK, n), no ties
    noties = true;
    n = arg1;
    exact = (n < 50);
else % pvalKendall(tail, K, stdK, xrank, yrank), ties in data
    noties = false;
    
    % If stdK is zero, at least one of the data vectors was constant.  The
    % correlation coef in these cases falls out of the calculations correctly
    % as NaN, but the exact p-value calculations below would regard Pr{K==0}
    % as 1.  Return a NaN p-value instead.
    K(stdK == 0) = NaN;
    
    xrank = arg1;
    yrank = arg2;
    n = length(xrank);
    exact = (n < 10);
end
nfact = factorial(n);
n2const = n*(n-1)/2;

if exact
    if noties
        % No ties, use recursion to get the cumulative distribution of
        % the number, C, of positive (xi-xj)*(yi-yj), i<j.
        %
        % K = #pos-#neg = C-Q, and C+Q = n(n-1)/2 => C = (K + n(n-1)/2)/2
        freq = [1 1];
        for i = 3:n
            freq = conv(freq,ones(1,i));
        end
        freq = [freq; zeros(1,n2const+1)]; freq = freq(1:end-1)';
        
    else
        % Ties, take permutations of the midranks.
        %
        % With ties, we could consider only distinguishable permutations
        % (those for which equal ranks (ties) are not simply interchanged),
        % but generating only those is a bit of work.  Generating all
        % permutations uses more memory, but gives the same result.
        xrank = repmat(xrank(:)',nfact,1);
        yrank = perms(yrank(:)');
        Kperm = zeros(nfact,1);
        for k = 1:n-1
            U = sign(repmat(xrank(:,k),1,n-k)-xrank(:,k+1:n));
            V = sign(repmat(yrank(:,k),1,n-k)-yrank(:,k+1:n));
            Kperm = Kperm + sum(U .* V, 2);
        end
        freq = histcounts(Kperm,-(n2const+.5):(n2const+.5)); freq = freq(1:end-1);
    end
    
    % Get the tail probabilities.  Reflect as necessary to get the correct
    % tail.
    switch tail
        case 'b' % 'both or 'ne'
            % Use twice the smaller of the tail area above and below the
            % observed value.
            tailProb = min(cumsum(freq), rcumsum(freq,nfact)) ./ nfact;
            tailProb = min(2*tailProb, 1); % don't count the center bin twice
        case 'r' % 'right' or 'gt'
            tailProb = rcumsum(freq,nfact) ./ nfact;
        case 'l' % 'left' or 'lt'
            tailProb = cumsum(freq) ./ nfact;
    end
    p = NaN(size(K),class(K));
    t = ~isnan(K(:));
    p(t) = tailProb(K(t) + n2const+1); % bins at integers, starting at -n2const
    
else
    switch tail
        case 'b' % 'both or 'ne'
            p = normcdf(-(abs(K)-1) ./ stdK);
            p = 2*p; p(p>1) = 1; % Don't count continuity correction at center twice
        case 'r' % 'right' or 'gt'
            p = normcdf(-(K-1) ./ stdK);
        case 'l' % 'left' or 'lt'
            p = normcdf((K+1) ./ stdK);
    end
end

%--------------------------------------------------------------------------
function y = rcumsum(x,sumx)
%RCUMSUM Cumulative sum in reverse direction. Note that this is the same as
% CUMSUM(X,'reverse') if SUMX==SUM(X), but can be used for cases where we
% want to count down from some other total too.
y = repmat(sumx,size(x));
y(2:end) = sumx - cumsum(x(1:end-1));