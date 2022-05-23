function ess = computeESS(xsmpl,maxlag,historysize)
%computeESS - Effective sample size (ESS) using Geyer's monotone estimator.
%   ESS = computeESS(XSMPL,MAXLAG) takes a P-by-N matrix of MCMC samples
%   where P is the number of parameters and N is the number of samples and
%   computes a P-by-1 vector ESS containing the effective sample size for
%   each parameter. MAXLAG is an integer specifying the maximum number of
%   lags to include in the estimation. ESS that exceeds N is recorded as N.
%
%   ESS = computeESS(XSMPL,MAXLAG,HSIZE) also accepts a scalar HSIZE such
%   that XSMPL(:,1:HSIZE+1:end) forms chain 1, XSMPL(:,2:HSIZE+1:end) forms
%   chain 2 and XSMPL(:,HSIZE+1:HSIZE+1:end) forms the last chain. In this
%   case, ESS is calculated separately for each chain and the result is the
%   sum over all chains.

    if nargin < 3
        historysize = 0;
    end

    if ( historysize == 0 )
        ess = computeGeyerMonotoneESS(xsmpl,maxlag);
    else
        P         = size(xsmpl,1);
        xsmplcell = cell(historysize+1,1);
        essqn     = zeros(P,historysize+1);
        
        for i = 1:(historysize+1)
            xsmplcell{i} = xsmpl(:,i:historysize+1:end);
            essqn(:,i)   = computeGeyerMonotoneESS(xsmplcell{i},maxlag);
        end
        ess = sum(essqn,2);
    end
end

function ess = computeGeyerMonotoneESS(xsmpl,maxlag)
    % 1. How many parameters and how many observations?
    [P,N] = size(xsmpl);
    
    % 2. Compute the raw auto-correlations. acf(:,i) corresponds to
    % parameter i.
    lag = max(0,min(maxlag,N-1));    
    
    acf = zeros(lag+1,P);    
    for i = 1:P
        acf(:,i) = stats.mcmc.utils.computeAutoCorr(xsmpl(i,:),lag);
    end
    
    % 3. Sum consecutive values of acf for each parameter.
    m   = floor(size(acf,1)/2);
    phi = zeros(m,P);
    for i = 1:P        
        for j = 1:m           
            phi(j,i) = acf(2*j-1,i) + acf(2*j,i);            
        end        
    end
    
    % 4. Compute minphi so minphi(j,i) is the smallest among phi(1:j,i).
    minphi = cummin(phi,1);
    
    % 5. Denominator for effective sample size. We could compute denom(i)
    % like this:
    %
    %   denom(i) = max(0,-1 + 2*sum(minphi(posidx,i)));
    %
    % This would allow ess(i) to be greater than N for some patterns of
    % alternating positive and negative autocorrelations. But it's possible
    % that -1 + 2*sum(minphi(posidx,i)) is < 0 for some posidx and i. In
    % that case, denom(i) will be 0 making ess(i) equal to Inf. To avoid
    % that, we require denom(i) >= 1. In other words, ess(i) that is
    % greater than N is reported as N.
    denom = zeros(P,1);
    for i = 1:P
        posidx   = minphi(:,i) > 0;
        denom(i) = max(1,-1 + 2*sum(minphi(posidx,i)));
    end
    
    % 6. Compute effective sample size.
    ess = N./denom;
end
