function ci = statnormci(parmhat,cv,alpha,x,censoring,freq) %#codegen %#internal
% Code generation compatible statnormci, a helper for distribution fitting functions
% STATNORMCI Confidence intervals for normal distribution

%   Copyright 2019 The MathWorks, Inc.


if isvector(parmhat)
    parmhat = parmhat(:);
end
muhat = parmhat(1,:);
sigmahat = parmhat(2,:);

if isempty(freq) || isequal(freq,1)
    if isvector(x)
        n = length(x);
    else
        n = size(x,1);
    end
else
    n = sum(freq);
end

if isempty(freq)
    if isempty(censoring)
        ncen = 0;
    else
        ncen = sum(censoring);
    end
    
else
    
    n = cast(sum(freq), 'like',x);
    if isempty(censoring)
        ncen = 0;
    else
        ncen = sum(freq.*censoring);
    end
    
end
nunc = n - ncen;

if any(censoring)
    % x will be a vector in this case
    if n == 0 || nunc == 0 || ~isfinite(parmhat(1))
        muci = coder.internal.nan(2,1);
        sigmaci = coder.internal.nan(2,1);
        ci = cast(cat(3,muci,sigmaci),'like',x);
        return
    end
end

% When all uncensored observations are equal and greater than all the
% censored observations, the likelihood surface becomes infinite at the
% boundary sigma==0.  Return something reasonable anyway.
if any(censoring)
    % x will be a vector in this case
    xunc = x(censoring==0);
    rangexUnc = stats.coder.internal.rangeWithCensoring(xunc,[]);
    if rangexUnc < realmin(class(x))
        if xunc(1) == max(x)
            if nunc > 1
                muci = [muhat;muhat];
                sigmaci = [0;0];
            else
                muci = [-coder.internal.inf;coder.internal.inf];
                sigmaci = [0;coder.internal.inf];
            end
            ci = cast(cat(3,muci,sigmaci),'like',x);
            return
        end
    end
end

% Get confidence intervals for each parameter
if (isempty(censoring) || ~any(censoring(:))) && ~isequal(cv,zeros(2,2))
    % Use exact formulas
    tcrit = tinv([alpha/2 1-alpha/2],n-1);
    muci = [muhat + tcrit(1)*sigmahat/sqrt(n); ...
        muhat + tcrit(2)*sigmahat/sqrt(n)];
    chi2crit = chi2inv([alpha/2 1-alpha/2],n-1);
    sigmaci = [sigmahat*sqrt((n-1)./chi2crit(2)); ...
        sigmahat*sqrt((n-1)./chi2crit(1))];
else
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(cv))';
    z = norminv(probs);
    
    % Compute the CI for mu using a normal distribution for muhat.
    muci = muhat + se(1).*z;
    
    % Compute the CI for sigma using a normal approximation for
    % log(sigmahat), and transform back to the original scale.
    % se(log(sigmahat)) is se(sigmahat) / sigmahat.
    logsigci = log(sigmahat) + (se(2)./sigmahat) .* z;
    sigmaci = exp(logsigci);
end

% Return as a single array
ci = cat(3,muci,sigmaci);
