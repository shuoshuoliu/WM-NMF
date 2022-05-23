function [coef,pval] = corrPearson(rows,tail,x,y)
%CORRPEARSON Computation of Pearson correlation matrix
%
%   [COEF,PVAL] = CORRPEARSON(ROWS,TAIL,X,Y) calculates the Pearson
%   cross-correlation matrix COEF. ROWS must be 'all', 'complete', or
%   'pairwise'. TAIL must be 'both', 'left', or 'right'.
%
%   See also: CORR, INTERNAL.STATS.CORRKENDALL, INTERNAL.STATS.CORRSPEARMAN.

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

% Computation of linear correlation with column-pairwise NaN removal
if pairwise
    % Preallocate the outputs.
    outType = internal.stats.dominantType(x,y);
    coef = zeros(p1,p2,"like",outType);
    nij = zeros(p1,p2,"like",outType);
    for i = 1:p1
        xi = x(:,i);
        isnanXi = isnanX(:,i);
        
        % Loop over all columns of y for corr(x,y), or over the lower triangle
        % without the diagonal for corr(x).
        j1 = p2*(1-corrXX) + (i-1)*corrXX; % 1:p2, or 1:(i-1)
        for j = 1:j1
            yj = y(:,j);
            ok = ~(isnanXi | isnanY(:,j));
            nij(i,j) = sum(ok);
            if nij(i,j) < n
                x0 = xi(ok); x0 = x0 - sum(x0)/nij(i,j);
                y0 = yj(ok); y0 = y0 - sum(y0)/nij(i,j);
            else
                x0 = xi - sum(xi)/n;
                y0 = yj - sum(yj)/n;
            end
            coef(i,j) = (x0'*y0) ./ (norm(x0) .* norm(y0));
        end
        % Compute diagonal elements for corr(x) to match the non-pairwise
        % case.  This is either 1, if the data in xi are not constant, or NaN
        % if they are.
        if corrXX
            ok = ~isnanXi;
            nij(i,i) = sum(ok);
            if nij(i,i) < n
                x0 = xi(ok); x0 = x0 - mean(x0);
            else
                x0 = xi - mean(xi);
            end
            coef(i,i) = (x0'*x0)./norm(x0)^2;
        end
    end
    
    % If this is autocorrelation, reflect the lower triangle into the upper.
    if corrXX
        coef = tril(coef) + tril(coef,-1)'; % leave the diagonal alone
    end
    
    n = nij;
    
    % Computation of linear correlation, all elements at once.  NaNs not removed,
    % or have already been removed in complete rows.
else
    x = x - sum(x,1)/n;  % Remove mean
    if corrXX
        coef = x' * x; % 1/(n-1) doesn't matter, renormalizing anyway
        d = sqrt(diag(coef)); % sqrt first to avoid under/overflow
        coef = coef./d; coef = coef./d'; % coef = coef ./ d*d';
    else
        y = y - sum(y,1)/n;  % Remove mean
        coef = x' * y; % 1/(n-1) doesn't matter, renormalizing anyway
        dx = vecnorm(x,2,1);
        dy = vecnorm(y,2,1);
        coef = coef./dx'; coef = coef./dy; % coef = coef ./ dx'*dy;
    end
end

% Limit off-diag elements to [-1,1], and put exact ones on the diagonal for
% autocorrelation.
t = find(abs(coef) > 1); coef(t) = coef(t)./abs(coef(t)); % preserves NaNs
if corrXX
    coef(1:p1+1:end) = sign(diag(coef)); % preserves NaNs
end

if nargout > 1
    if corrXX
        pval = zeros(p1,"like",coef);
        ltri = (tril(ones(size(pval)),-1) > 0); % lower triangle only, no diagonal
        if pairwise
            pval(ltri) = pvalPearson(tail, coef(ltri), n(ltri));
        else
            pval(ltri) = pvalPearson(tail, coef(ltri), n);
        end
        % Reflect the p-values from the lower triangle into the upper.
        pval = pval + pval';
        % The p-values along the diagonal are always exactly one (unless
        % they're NaN), because Pr{coef(i,i) as/more extreme than 1} == 1,
        % regardless of which tail(s) we're testing against.
        pval(1:p1+1:end) = sign(diag(coef)); % preserves NaNs on diag
    else
        pval = pvalPearson(tail, coef, n);
    end
end


%--------------------------------------------------------------------------

function p = pvalPearson(tail, rho, n)
%PVALPEARSON Tail probability for Pearson's linear correlation.
t = rho.*sqrt((n-2)./(1-rho.^2)); % +/- Inf where rho == 1
switch tail
    case 'b' % 'both or 'ne'
        p = 2*tcdf(-abs(t),n-2);
    case 'r' % 'right' or 'gt'
        p = tcdf(-t,n-2);
    case 'l' % 'left' or 'lt'
        p = tcdf(t,n-2);
end
