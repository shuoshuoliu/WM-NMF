function [L, R, mu, D, final] = alsmf(Y, k, varargin)
%Alternating Least Squares (ALS) low-rank matrix factorization.
%   [L,R] = alsmf(Y,K) finds the best rank-K approximation by factoring a
%   N-by-P matrix Y into a N-by-K left factor matrix L and a P-by-K right
%   factor matrix R. The factorization uses an iterative method starting
%   with random initial values for L and R.
%
%   [L,R,MU] = alsmf(Y,K) also returns MU, the estimated mean.
%
%   [L,R,MU,D] = alsmf(Y,K) also returns D, the K largest eigenvalues of
%   the covariance matrix of Y. It returns NaNs if 'Orthonormal' is not
%   specified.
%
%   [L,R,MU,D,S] = alsmf(Y,K) also stores final results at convergence
%   in the struct FINAL:
%
%       S.rmsResid   - Root mean squared residuals.
%       S.NumIter    - Number of iterations.
%
%   [...] = alsmf(Y,K,L0,R0,'PARAM1',val1,'PARAM2',val2,...) specifies one
%   or more of the following parameter name/value pairs:
%
%       'L0'         - Initial value for L, a N-by-K matrix. The default is
%                      a random matrix.
%
%       'R0'         - Initial value for R, a P-by-K matrix. The default is
%                      a random matrix.
%
%       'Orthonormal'- Indicate if rows of returned matrix R are
%                      orthonormal (in classic PCA basis).
%
%         'Weights'  - Observation weights, a vector of length N containing
%                      all positive elements.
%
%  'VariableWeights' - Variable weights. A vector of length P containing
%                      all positive elements.
%
%          'Options' - An options structure as created by the STATSET
%                      function. ALSMF uses the following fields:
%             'Display' - Level of display output.  Choices are 'off' (the
%                         default), 'final', and 'iter'.
%             'MaxIter' - Maximum number of steps allowed. The default is
%                         1000. Unlike in optimization settings, reaching
%                         MaxIter is regarded as convergence.
%              'TolFun' - Positive number giving the termination tolerance
%                         for the cost function.  The default is 1e-6.
%                'TolX' - Positive number giving the convergence threshold
%                         for relative change in the elements of L and R.
%                         The default is 1e-6.
%
%
%   Example:
%      >> load fisherires;
%      >> [score,coeff,mu] = alsmf(meas,4,'Orth',true);
%      >> [score,coeff,mu] = alsmf(meas,3,'L0',rand(size(meas,1),3));
%
%   Reference: 
%   [1] Sam Roweis. 1998. EM algorithms for PCA and SPCA. In Proceedings
%       of the 1997 conference on Advances in neural information processing
%       systems 10 (NIPS '97), Michael I. Jordan, Michael J. Kearns, and
%       Sara A. Solla (Eds.). MIT Press, Cambridge, MA, USA, 626-632
%   [2] Alexander Ilin and Tapani Raiko. 2010. Practical Approaches to
%       Principal Component Analysis in the Presence of Missing Values. J.
%       Mach. Learn. Res. 11 (August 2010), 1957-2000.
%
%   Copyright 2012 The MathWorks, Inc.


[n, p ]= size(Y); % no error checking on Y and k

inpar = inputParser;
inpar.addParameter('L0',randn(n,k,'like',Y),@(x)isnumeric(x)&&isequal(size(x),[n,k]));
inpar.addParameter('R0',randn(p,k,'like',Y),@(x)isnumeric(x)&&isequal(size(x),[p,k]));

inpar.addParameter('Orthonormal',true,@islogical);
inpar.addParameter('Centered',true,@islogical);

chkweightsfunc = @(x,siz)isnumeric(x)&&length(x)==siz&&all(x>=0);
inpar.addParameter('Weights',ones(n,1,'like',Y),@(x)chkweightsfunc(x,n));
inpar.addParameter('VariableWeights',ones(1,p,'like',Y),@(x)chkweightsfunc(x,p));

dftOpt = statset('alsmf');
inpar.addParameter('Options',dftOpt);

inpar.parse(varargin{:});

MaxIter = statget(inpar.Results.Options,'MaxIter',dftOpt,'fast');
TolX = statget(inpar.Results.Options,'TolX',dftOpt,'fast');
TolFun = statget(inpar.Results.Options,'TolFun',dftOpt,'fast');
displayOn = statget(inpar.Results.Options,'Display',dftOpt,'fast');
flagOrthogonal = inpar.Results.Orthonormal;
flagCentered = inpar.Results.Centered;
Omega = inpar.Results.Weights;
Phi = inpar.Results.VariableWeights;

dispnum = strcmp(displayOn,{'off','iter','final'});

Omega = Omega(:);% Make it a column vector
Phi = Phi(:);% Make it a column vector

% Record Missing Values Locations
idxMiss = isnan(Y);

% Initialize
R = inpar.Results.R0;
L = inpar.Results.L0;
mu = zeros(1,p);

final.NumIter =0;
dnormOld = inf;

header = 'Iteration\t   |Delta X|\t   rms resid\n';
fmtstr = '   %3d   \t%12g\t%12g\n';
cvgCond = '';
if dispnum(2)
    fprintf(header);
end

for iter = 1: MaxIter
    % Update parameters    
    munew = localMu(L,R);
    Rnew = localR(L,munew);
    Lnew = localL(Rnew,munew);
    
    % Check convergence
    diff = bsxfun(@minus,Y-Lnew*Rnew',munew);
    diff(idxMiss) = 0;
    dnorm = norm(diff,'fro') / (sqrt(n*p-sum(idxMiss(:))));
    dL = max(max(abs(L-Lnew) / (sqrt(eps)+max(max(abs(Lnew))))));
    dR = max(max(abs(R-Rnew) / (sqrt(eps)+max(max(abs(Rnew))))));
    dmu = max(abs(mu-munew)) / (sqrt(eps)+max(abs(munew)));
    delta = max([dL,dR,dmu]);
    
    if (dnormOld - dnorm) < TolFun*dnorm
        cvgCond = getString(message('stats:pca:TolFunReached'));
        break;
    elseif delta < TolX
        cvgCond = getString(message('stats:pca:TolXReached'));
        break;
    elseif iter == MaxIter;
        warning(message('stats:pca:MaxIterReached', MaxIter));
        break;
    else
        % Record last iteration results if continue iteration
        dnormOld = dnorm;
        L = Lnew;
        R = Rnew;
        mu = munew;
        if dispnum(2)
        fprintf(fmtstr,iter,delta,dnorm);
        end
    end
end

D = zeros(k,1);
if flagOrthogonal   
    if flagCentered
        muL = classreg.learning.internal.wnanmean(L,Omega);
        L = bsxfun(@minus,L,muL);
        mu = mu+muL*R';
    end
    
    [U, Dr] = eig(bsxfun(@times,R,Phi)'*R);
    Dr = diag(Dr);
    Adj = bsxfun(@times,L*U,sqrt(Dr)');
    [V,D] = eig(bsxfun(@times,Adj,Omega)'*Adj);
    
    R = bsxfun(@times,R*U,1./sqrt(Dr)')*V;
    L = Adj*V;
    D = diag(D)/(size(L,1)-flagCentered);
    
    [D,idxD] = sort(D,'descend');
    L = L(:,idxD);
    R = R(:,idxD);
end

% Display finial output:
if ~dispnum(1)
    finalstr = getString(message('stats:pca:FinalDisplay',iter,cvgCond,...
        sprintf('%g',dnorm)));
    fprintf('\n%s\n',finalstr);
end

% Store final results at convergence
if nargout> 4
    final.rmsResid = dnorm;
    final.Recon = bsxfun(@plus,L*R',mu);
    final.NumIter = iter;
end


%--- Updating Functions ---
    function L = localL(R,mu)
        %LOCALL computes LS estimate of L assuming R is fixed.
        
        L = zeros(n,k); %initialization
        for i = 1:n
            idx = ~idxMiss(i,:);
            L(i,:) = locallscov(R(idx,:),(Y(i,idx)-mu(idx))',Phi(idx))';

        end
    end % End of localL function

    function R = localR(L,mu)
        %LOCALR computes LS estimate of R assuming L is fixed.
        
        R = zeros(p,k);
        for j = 1:p
            idx = ~idxMiss(:,j);
            R(j,:) = locallscov(L(idx,:),Y(idx,j)-mu(j),Omega(idx));
        end
    end % End of localR function

    function mu = localMu (L,R)
        %LOCALMU updates the means
        
        diff = Y-L*R';
        if ~flagCentered
            mu = zeros(1,p,'like',Y);
        else
            mu = classreg.learning.internal.wnanmean(diff,Omega);
        end
    end % localMu

end % End of main function

%---------Sub Function----------
function x = locallscov(A,B,w)
%Weighted Least Squares. w must be a column vector

C=(bsxfun(@times,A,w)'*A);
if condest(C) > 1/sqrt(eps(class(C)))
    x = pinv(C)*(bsxfun(@times,A,w)'*B);
else
    x = C\(bsxfun(@times,A,w)'*B);
    % x = lscov(A,B,w);
end
end% End of locallscov