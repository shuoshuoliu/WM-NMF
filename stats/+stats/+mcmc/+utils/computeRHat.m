function rhat = computeRHat(chains)
%computeRHat - Computes RHAT diagnostic for accessing MCMC convergence.
%   RHAT = computeRHat(CHAINS) computes a diagnostic RHAT using the cell
%   array CHAINS. CHAINS{i} is a P-by-Ni matrix of Ni samples in P
%   dimensions. RHAT is a P-by-1 vector containing the potential scale
%   reduction statistic for each dimension.
%
%   RHAT = computeRHat(CHAINS) computes a scalar diagnostic RHAT using the
%   M-by-N matrix CHAINS. M is the number of chains and N is the number of
%   samples per chain.

    if iscell(chains)
        % 1. Make chains a row cell array.
        chains = chains(:)';

        % 2. Get the number of chains M.
        M = length(chains);

        % 3. means and vars are P-by-M matrices and Nvals is a 1-by-M vector.
        means = cell2mat(cellfun(@(x) mean(x,2) ,chains,'UniformOutput',false));
        vars  = cell2mat(cellfun(@(x) var(x,0,2),chains,'UniformOutput',false));
        Nvals = cell2mat(cellfun(@(x) size(x,2) ,chains,'UniformOutput',false));

        % 4. Compute the mean of means and vars across chains for each
        % dimension. Variables below are of size P-by-1.
        meanmeans = mean(means,2);
        W         = mean(vars,2);

        % 5. Numerator term for rhat of size P-by-1.
        pooledvar = mean(bsxfun(@times,(Nvals-1)./Nvals,vars),2) + sum(bsxfun(@minus,means,meanmeans).^2,2)/(M-1);

        % 6. Compute rhat.
        rhat = sqrt(pooledvar./W);
    else
        [m,n] = size(chains);

        theta_idot   = mean(chains,2);
        theta_dotdot = mean(theta_idot);

        B = n*sum((theta_idot - theta_dotdot).^2)/(m-1);

        Si2 = var(chains,0,2);
        W   = mean(Si2);

        rhat = sqrt(1 + (1/n)*(B/W - 1));
    end
    rhat = max(1,rhat);
end