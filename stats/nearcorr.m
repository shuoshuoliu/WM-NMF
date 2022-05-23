function Y = nearcorr(A, varargin)
% NEARCORR Nearest correlation matrix
%   Compute the nearest correlation matrix Y by minimizing the Frobenius
%   distance.
%
%   Y = NEARCORR(A) returns a N-by-N nearest correlation matrix for
%   specified N-by-N symmetric matrix A.
%
%   Y = NEARCORR(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values.
%
% Inputs:
%
%   A               - A is an N-by-N symmetric approximate correlation
%                     matrix with all elements in [-1 1] and unit diagonal.
%
% Optional Inputs (Parameter Name/Value Pairs):
%
%   Tolerance       - Positive numerical tolerance for the iteration. 
%                     The default is 1e-6.
%
%   MaxIterations   - Positive integer, maximum number of solver iterations.
%                     The default is 200.
%
%   Method          - Method applied to solve nearest correlation matrix
%                     problem. Possible values are "newton"(default) and
%                     "projection", which represent "Newton algorithm" and
%                     "Alternating Projections algorithm" respectively.
%
%   Weights         - Weights may be specified as one of two types:
%
%                     As a symmetric matrix W with all elements >=0 to do
%                     element-wise weighting and the nearest correlation
%                     matrix Y will be computed by minimizing the norm of
%                     (W.*(A-Y)). Larger weight values place greater
%                     importance on the corresponding elements in A.
%	              
%                     Or as an N x 1 vector w with positive numerical
%                     values. In this case, the nearest correlation
%                     matrix Y will be computed by minimizing the norm of
%                     (diag(w)^0.5*(A-Y)*diag(w)^0.5)
%
%                     When "newton" method is specified, the Weights can be
%                     either a symmetric matrix or an N x 1 vector.
%                     When "projection" method is specified, the Weights
%                     can be an N x 1 vector only.
%
% Output:
%
%   Y               - Nearest correlation matrix to the input A
%

% Copyright 2019 The MathWorks, Inc.

%
% References:
%   [1] Higham, N.J., "Computing the nearest correlation matrix - a problem
%       from finance", IMA Journal of Numerical Analysis, Vol. 22, Issue. 3,
%       2002. 
%   [2] Pang, J.S., Sun, D., Sun, J. "Semismooth homeomorphisms and strong
%       stability of semidefinite and Lorentz complementarity problems",
%       Mathematics of Operation Research, Vol 28, 2003.
%   [3] Qi, H., Sun, D., "A quadratically convergent Newton method for
%       computing the nearest correlation matrix", SIAM J. Matrix Anal.
%       Appl., Vol. 28, 2006. 
%   [4] Borsdorf, R., Higham, N.J., "A preconditioned Newton algorithm for
%       the nearest correlation matrix", IMA Journal of Numerical Analysis,
%       Vol. 30, 2010.
%   [5] Qi, H., Sun, D., "An augmented Lagrangian dual approach for the
%       H-weighted nearest correlation matrix problem", IMA Journal of
%       Numerical Analysis, Vol 31, Issue 2, 2011.


if nargin < 1
    error(message('stats:nearcorr:TooFewInputs'));
end

if ~isnumeric(A) || isempty(A) || ~isreal(A) || ~ismatrix(A) || ~issymmetric(A) || min(min(A)) < -1 || max(max(A)) > 1 || any(diag(A)-1)
    error(message('stats:nearcorr:InvalidInputA'));
end

p = inputParser;
p.addParameter('Tolerance',1e-6,@(x)validateattributes(x, {'double'}, {'scalar' 'positive'},'','Tolerance'));
p.addParameter('MaxIterations',200,@(x)validateattributes(x, {'double'}, {'scalar' 'positive' 'integer'},'','MaxIterations'));
p.addParameter('Method',"newton", @(x)(isscalar(string(x))&&ismember(lower(x), ["newton","projection"])));
p.addParameter('Weights',[], @(x)validateattributes(x, {'double'}, {'real'},'','Weights'));

try
    p.parse(varargin{:});
catch ME
    newMsg = message('stats:nearcorr:optionalInputError');
    newME = MException(newMsg.Identifier,getString(newMsg));
    newME = addCause(newME,ME);
    throw(newME)
end

tol = p.Results.Tolerance;
maxIter = p.Results.MaxIterations;
method = p.Results.Method;
weight = p.Results.Weights;

if ~issymmetric(A)
    A = (A + A') / 2;
end

[Ndata, ~] = size(A);

if strcmpi(method, "projection")
    if ~isempty(weight)
        % Alternating Projections method only accepts diagonal W-form
        % weighting or (W^0.5*(A-Y)*W^0.5)
        if ~isvector(weight) || length(weight) ~= Ndata || min(weight) <= 0
            error(message('stats:nearcorr:InvalidWeightsProjection', Ndata));
        end
    else
        % no weighting
        weight = 1;
    end
    
else % Newton method
    if ~isempty(weight)
        % Newton method accepts diagonal W-form weighting - (W^0.5*(A-Y)*W^0.5)
        % or H-form weighting - (W .* (A-Y))
        
        if isvector(weight) 
            % for diagonal W-form weighting
            if length(weight) ~= Ndata || min(weight) <= 0
                error(message('stats:nearcorr:InvalidWeightsNewton', Ndata));
            end
        else
            % for element-wise H-form weighting
            if ~ismatrix(weight) || ~issymmetric(weight) || size(weight,1) ~= Ndata || min(min(weight)) < 0
                error(message('stats:nearcorr:InvalidWeightsNewton', Ndata));
            end
        end
        
    else
        % no weighting
        weight = 1;
    end
    
end

maxIterReached = true;
Y = 0;

% Alternating Projections method
if strcmpi(method, "projection")
    % Convert diagonal W-form weighting to element-wise matrix multiplication
    weight = weight(:);
    weight = (weight * weight').^0.5;
    
    % Set initial conditions
    dS = 0;
    Yold = A;
    Xold = A;
    
    diagIdx = 1:Ndata+1:Ndata^2;
    
    for iter = 0 : maxIter
        % add Dykstra's correction
        R = Yold - dS;        
        
        % do spectral decomposition and projection onto positive
        % semi-definite subspace
        [Q,D] = eig(weight .* R);
        D = max(sparse(D),0);
        X = (Q * D * Q') ./ weight;
        X = (X + X')/2;
        dS = X - R;        
        Y = X;
        
        % projection onto unit diagonal subspace
        Y(diagIdx) = 1;
        
        % calculate the exit conditions
        normY = norm(Y, 'fro');
        cond1 = norm(Y - Yold, 'fro')/normY;
        cond2 = norm(X - Xold, 'fro')/norm(X, 'fro');
        cond3 = norm(Y - X, 'fro')/normY;
                
        if  max([cond1, cond2, cond3]) <= tol
            maxIterReached = false;
            break
        end
        
        Yold = Y;
        Xold = X;
                
    end
    
    % restore the unit diagonal
    scaleMat = diag(X);
    scaleMat = (scaleMat * scaleMat').^0.5;
    Y = X ./ scaleMat;
    Y = full(Y + Y')/2;
    
end

% Newton method
if strcmpi(method, "newton")
    if isscalar(weight) || isequal(weight, ones(Ndata))
        % No weighting
        weight = 1;
        [Y, maxIterReached] = nearcorr_newton(A, Ndata, weight, tol, maxIter);
        
    else
        if isvector(weight)
            % for diagonal W-form weighting - (W^0.5*(A-Y)*W^0.5)
            % Convert diagonal W-form weighting to element-wise matrix multiplication
            weight = weight(:);
            weight = (weight * weight').^0.5;
        end
        
        % for element-wise H-form weighting - (W .* (A-Y))
        
        % set parameters
        c = 10;
        k = 1.4;
        eta = 1e-2;
        mu = 1e-12;
        tao1 = 1e-2;
        tao2 = 1e1;
        tao3 = 1e4;
        rho = 0.5;
        weight = (weight + weight')/2;
                        
        % Get initial guess by solving basic problem without weighting
        [X0, maxIterReached, y0] = nearcorr_newton(A, Ndata, 1, tol, maxIter);
        Z0 = X0 - A - diag(y0); % compute Z0
        Z0 = (Z0 + Z0') / 2; % make sure Z is symmetric
        
        [Q,D] = eig(Z0-X0);
        D = max(sparse(D),0);
        KKT1 = weight .* weight .* (X0 - A) - diag(y0) - Z0;
        KKT2 = 1 - diag(X0); % b - diag(X), b==e
        KKT3 = Z0 - (Q * D * Q');
        Tolk = max([norm(KKT1,'fro'), norm(KKT2,'fro'), norm(KKT3,'fro')]); % first KKT distance
        
        Z = Z0;
        X = X0;
        y = y0;
        
        % calculate gradient
        grad = gradient_Hform(A, y, weight, X, Z, c);
        cMax = norm(grad,'fro');
        c = max(c, cMax/.5e2);
        
        if Tolk > tol 
            % not meet the specified tolerance, proceed with augmented
            % Lagrangian dual approach for the element-wise H-form
            % weighting or a semismooth Newton method
            maxIterReached = true;
            
            for iter = 0 : maxIter % user specified iteration#
                
                [grad, theta_XyZ, Q, D, ~] = gradient_Hform(A, y, weight, X, Z, c);
                
                for solveX = 1 : 40
                    
                    % Get the sizes of positive, zero and negative eigenvalues to build W
                    alpha = find(D > 0);
                    beta = find(D == 0);
                    gamma = find(D < 0);
                    
                    Nalpha = numel(alpha);
                    Nbeta = numel(beta);
                    Ngamma = numel(gamma);
                    
                    lamda = D(alpha,ones(Ngamma,1))./(D(alpha,ones(Ngamma,1))-D(gamma,ones(Nalpha,1))');
                    
                    [Nrow,~] = size(lamda);
                    if Nrow ~= Nalpha
                        lamda = reshape(lamda,Nalpha,Ngamma);
                    end
                    
                    normGrad = norm(grad,'fro');
                    sj = min(tao1, tao2*normGrad);                    
                    tolcg = min(eta, tao3*normGrad) * normGrad;
                    M = precondMat(lamda, Q, Nalpha, Nbeta, Ngamma, Ndata);
                    
                    % solve the Newton equations using preconditioned CG
                    % and approximate the search direction
                    % without explicitly computing and storing Hessian matrix
                    dX = cg(-grad, weight, c, M, lamda, Q, sj, Ndata, tolcg, Nalpha, Nbeta, Ngamma);
                    
                    
                    cond = sum(sum(grad .* dX));
                    maxSearch = 20;                    
                    for m = 0:maxSearch
                        
                        [grad_m, theta_XyZ_m, Q, D, ~] = gradient_Hform(A, y, weight, X + rho^m * dX, Z, c);
                        
                        a = rho^m;
                        if theta_XyZ_m - theta_XyZ <= mu * rho^m * cond                            
                            break
                        end
                        
                         if abs(theta_XyZ_m - theta_XyZ) < 1e4 * eps * (1 + abs(theta_XyZ_m) + abs(theta_XyZ))
                             if norm(grad_m, 'fro')/normGrad <= 1 - 1e-2
                                 a = 1;
                             else
                                 dX = - grad;
                                 a = 1;
                             end
                             
                             [grad_m, theta_XyZ_m, Q, D, ~] = gradient_Hform(A, y, weight, X + a * dX, Z, c);
                            break 
                         end
                    end
                    
                    X = X + a * dX; % update X(k+1)
                    X = (X + X') / 2; % make sure X is symmetric
                    
                    grad = grad_m;
                    theta_XyZ = theta_XyZ_m;
                    
                    if norm(grad_m,'fro') <= min(10, 0.5 * Tolk) % X(k+1) solved                                              
                        break
                    end
                    
                end
                
                % update y, Z in the dual problem
                y = y + c * (1 - diag(X)); % update y(k+1) = y(k) + c(k)*(b-diag(X(k+1))), b == e;
                [Q,D] = eig(Z - c * X);
                D = max(sparse(D),0);
                Zold = Z;
                Z = (Q * D * Q'); % update Z(k+1) = (Z(k) - c(k)*X(k+1))+
                Z = (Z + Z') / 2; % make sure Z is symmetric
                
                % update Tolerance(k+1)
                % three terms correspond to the three conditions in the KKT
                % system
                Tolkold = Tolk;
                Tolk = max([norm(grad,'fro'), ...
                            norm(1 - diag(X), 'fro'), ...            % KKT(2): b - diag(X), b == e
                            norm(Z - Zold, 'fro')/nthroot(c, 2)]);
                
                if Tolk <= tol
                    % the whole augmented Lagrangian algorithm ends
                    maxIterReached = false;                    
                    break
                end
                
                % update c(k+1)                
                if Tolk > 0.25 * Tolkold
                    c = min(k * c, cMax);
                end
                
            end
            
            % restore the unit diagonal
            scaleMat = diag(X);
            scaleMat = (scaleMat * scaleMat').^0.5;
            Y = X ./ scaleMat;
            Y = full(Y + Y')/2;
            
        else
            Y = X0;
        end
        
    end
    
end

% handle computational rounding error causing small negative eigenvalues
minE = min(eig(Y));
if minE < 0
    Y = Y ./ (1 - minE + eps);
    Y(1:Ndata+1:Ndata^2) = 1;
    
    minE = min(eig(Y));
    n = 10;
    while minE < 0
        Y = Y ./ (1 + n*eps);
        Y(1:Ndata+1:Ndata^2) = 1;
        
        minE = min(eig(Y));
        n = n * 10;
    end
end

if maxIterReached
   warning(message("stats:nearcorr:MaxIter"));
end

end

function B = jacobMat(x, Q, Wag, weight, Nalpha, Nbeta, Ngamma, Ndata)
% Calculate the system V*h = diag(Q * (W .* (Q' * diag(h) * Q)) * Q') 

% Computing the generalized Jacobian V of gradient explicitly which is
% positive definite is prohibitively expensive and not necessary. Instead,
% the minres method is used to require matrix and vector product V*h only

weight = diag(weight);

if Nalpha == 0
    B = zeros(Ndata,1);
    return
elseif Nalpha == Ndata
    B = x./weight;
    return
end

if Nalpha <= Ngamma
    ind1 = 1:Nalpha;
    ind2 = Nalpha+1:Ndata;
        
    Wag = [ones(Nalpha,Nbeta), Wag];
    
    Ql = Q(:, ind1);
    Qr = Q(:, ind2);
    
    FG = bsxfun(@times, Ql, x);
    H = FG' * Qr;
    WH = Wag .* H; % G = H' since x is symmetric and Q is orthogonal;
    WF = Ql' * FG;
    
    FH = WF * Ql' + WH * Qr';
    GJ = WH' * Ql';
    
    B = sum([FH; GJ] .* Q', 1)'./weight; % diag(Q * B)
    
else
    % utilize eigenvector sign ambiguity
    ind1 = 1:Nalpha+Nbeta;
    ind2 = Nalpha+Nbeta+1:Ndata;
    
    Wag =  [ones(Nbeta,Ngamma); 1-Wag];
    
    Ql = Q(:, ind1);
    Qr = Q(:, ind2);
    
    HJ = bsxfun(@times, Qr, x);
    G = HJ' * Ql;
    WG = Wag' .* G; % G = H' since x is symmetric and Q is orthogonal;
    WJ = Qr' * HJ;
    
    GJ = WJ * Qr' + WG * Ql';
    FH = WG' * Qr';
    
    B = (x - sum([FH; GJ] .* Q', 1)')./weight; % diag(Q * B)
end
        
end


function [grad, theta_y, Q, D, maxD] = gradient(A, y, weight)
% Compute nearest positive semidefinite matrix to A + diag(y) and the
% gradient

[Q,D] = eig(weight .* (A + diag(y)));
[D,ind] = sort(diag(D),'descend');
Q = Q(:,ind);

maxD = max(sparse(D),0);

weight = diag(weight);
grad = Q * diag(maxD);
grad = sum(Q .* grad, 2)./weight - 1; % gradient(y) = diag(specDecompPlus_y) - b with b == e

normFsq = sum((maxD./weight).^2); % Frobenius norm square == sum(eigenvalues^2)
theta_y = 0.5 * normFsq - sum(y, 1); % 0.5 * F-norm^2 - e'*y, the objective function

end

function [grad, theta_XyZ, Q, D, maxD] = gradient_Hform(A, y, weight, X, Z, c)
% Compute nearest positive semidefinite matrix to Z - c*X and the gradient

normZ = norm(Z, 'fro');
Z = Z - c * X;
Z = (Z + Z') / 2;

[Q,D] = eig(Z);
[D,ind] = sort(diag(D),'descend');
Q = Q(:,ind);

maxD = max(sparse(D),0);
S = (Q * diag(maxD) * Q');

dy = 1 - diag(X); % b-diag(X) with b == e

grad = weight .* weight .* (X-A) - diag(y + c*dy) - S; % gradient
theta_XyZ = 0.5 * norm(weight.*(X-A), 'fro')^2 + y'*dy + 0.5*c*sum(dy.^2) ...
            + 0.5/c * (norm(S, 'fro')^2 - normZ^2); % the objective function

end

function M = precondMat(Wag, Q, Nalpha, Nbeta, Ngamma, Ndata)
%Calculate the Jacobi preconditioner for solving the V*h system
%Preconditioner 1. Mtmp = W * (Q.^2)' 
%               2. M = diag(V) = sum((Q.^2)' .* Mtmp, 1)

if Nalpha == 0 || Nalpha == Ndata
    M = eye(Ndata);
    return
end

Q2 = (Q.^2)'; % Q square = Q'.*Q';

% utilize the block matrix multiplication and flip the sign of eigenvectors
% depending on the size of Nalpha*Nalpha since matrix W is a sparse matrix
% W = [ones(Nalpha,Nalpha) ones(Nalpha,Nbeta) Wag(Nalpha, Ngamma)
%      ones(Nbeta,Nalpha)  zeros ...          zeros
%      Wag'(Ngamma,Nalpha) zeros ...          zeros              ]

if Nalpha <= Ngamma
    ind1 = 1:Nalpha;
    ind2 = Nalpha+1:Ndata;
        
    Wag = [ones(Nalpha,Nbeta), Wag];
    
    Qu = Q2(ind1, :);
    Qb = Q2(ind2, :);
    
    % block matrix multiplication and implicit expansion
    FH = sum(Qu,1) + Wag*Qb; % W*Q2 = [F H; G J]
    GJ = Wag'*Qu;
    
    M = sum([Qu.*FH; Qb.*GJ], 1); % preconditioner Mtmp = W * (Q.^2)' , M = diag(V) = sum((Q.^2)' .* Mtmp, 1)
else
    % utilize sign ambiguity
    ind1 = 1:Nalpha+Nbeta;
    ind2 = Nalpha+Nbeta+1:Ndata;
    
    Wag =  [ones(Nbeta,Ngamma); 1-Wag];
    
    Qu = Q2(ind1, :);
    Qb = Q2(ind2, :);
    
    % block matrix multiplication and implicit expansion
    GJ = sum(Qb,1) + Wag'*Qu; % W*Q2 = [F H; G J]
    FH = Wag*Qb;
    
    M = 1- sum([Qu.*FH; Qb.*GJ], 1); % preconditioner Mtmp = W * (Q.^2)' , M = diag(V) = sum((Q.^2)' .* Mtmp, 1)
end

M(M <= 1e-5) = 1e-5; % M must be positive definite
M = sparse(diag(M));
end

function [Y, maxIterReached, y, iter] = nearcorr_newton(A, Ndata, weight, tol, maxIter)
% Compute the nearest correlation matrix with inexact smoothing Newton
% method

maxIterReached = true;
y = zeros(Ndata,1);
eta = .5;
rho = 0.5;
sigma = 2e-4;

% Compute nearest positive semidefinite matrix to A + diag(y) and the
% gradient
[grad, theta_y, Q, D, maxD] = gradient(A, y, weight);

for iter = 0 : maxIter
    normGrad = norm(grad,'fro');
    
    if normGrad <= tol
        maxIterReached = false;
        break
    end
    
    % Get the sizes of positive, zero and negative eigenvalues to build W
    alpha = find(D > 0);
    beta = find(D == 0);
    gamma = find(D < 0);
    
    Nalpha = numel(alpha);
    Nbeta = numel(beta);
    Ngamma = numel(gamma);
    
    tao2 = D(alpha,ones(Ngamma,1))./(D(alpha,ones(Ngamma,1))-D(gamma,ones(Nalpha,1))');
    
    [Nrow,~] = size(tao2);
    if Nrow ~= Nalpha
        tao2 = reshape(tao2,Nalpha,Ngamma);
    end
    
    % preconditioning the Newton equation system
    M = precondMat(tao2, Q, Nalpha, Nbeta, Ngamma, Ndata);
    
    % A generalized Jacobian V_iter is positive semidefinite for all iters
    % and positive definite for sufficiently large iteration number iter
    % solve the Newton equations using minres and approximate the search
    % direction without explicitly computing and storing Jacobian matrix of
    % gradient
    [d,flag] = minres(@(x)jacobMat(x, Q, tao2, weight, Nalpha, Nbeta, Ngamma, Ndata), -grad, min(eta, normGrad), [], M);
    
    norm_d = norm(d,'fro');
    if ~(flag == 0 && (grad' * d <= -min(1e-6, normGrad) * norm_d^2))
        d = - grad;
    end
    
    % line search
    cond = grad' * d;
    maxSearch = 20;
    a = [];
    for m = 0:maxSearch
        [grad_y_m_d, theta_y_m_d, Q, D, maxD] = gradient(A, y + rho^m * d, weight);
        
        if theta_y_m_d - theta_y <= sigma * rho^m * cond
            a = rho^m;
            break
        end
        
        if abs(theta_y_m_d - theta_y) < 1e4 * eps * (1 + abs(theta_y_m_d) + abs(theta_y))
            if norm(grad_y_m_d, 'fro')/normGrad <= 1 - 1e-2
                a = 1;
            else
                d = - grad;
                a = 1;
            end
            
            [grad_y_m_d, theta_y_m_d, Q, D, maxD] = gradient(A, y + a * d, weight);
            break
        end
    end
    
    if isempty(a)
        d = - grad;
        a = 1;
        [grad_y_m_d, theta_y_m_d, Q, D, maxD] = gradient(A, y + a * d, weight);
    end
    
    % update y in the dual problem
    y = y + a * d;
    grad = grad_y_m_d;
    theta_y = theta_y_m_d;
    
end

% restore the unit diagonal
Y = Q * diag(maxD) * Q';
scaleMat = diag(Y);
scaleMat = scaleMat(:);
scaleMat = (scaleMat * scaleMat') .^ 0.5;
Y = Y ./ scaleMat;
Y = full(Y + Y') / 2;
end

function x = cg(b, weight, c, M, Wag, Q, sj, Ndata, tol, Nalpha, Nbeta, Ngamma)
% Preconditioned CG method
% The Newton system is (V + Sj*I)dX = -gradient(objectiveFun(X)). X is a 
% matrix from the element-wise weighting problem.
% V + Sj*I is positive definite so the CG method can always find dX to get
% a new search direction.

r = b; % zero matrix as the initial guess

M = diag(weight.*weight + c + sj + c*M); 
M = 1./M;
x = zeros(Ndata);

normr = norm(r,'fro');
if (normr <= tol)    
    return
end

for i = 1:40
    z = bsxfun(@times, r, M);  % z = M \ r;
    
    sumrz = sum(sum(r.*z));
    
    if i == 1        
        p = z;        
    else       
        beta =  sumrz ./  sumrzold;
        p = z + beta .* p;
    end
    
    p = (p + p')/2;
    % calculate (V + Sj*I)dX
    A = weight .* weight .* p + c * diag(diag(p)) + c * HessianMat(p, Q, Wag, Nalpha, Nbeta, Ngamma, Ndata) + sj * p;
        
    sumrzold = sumrz;
    
    pAp = sum(sum(p .* A));
    
    if pAp == 0
        x = p;
        break
    end
        
    alpha = sumrz ./ pAp;
    x = x + alpha .* p;
    r = r - alpha .* A;
        
    if norm(r, 'fro') < tol
        break;
    end    
end


    function B = HessianMat(x, Q, Wag, Nalpha, Nbeta, Ngamma, Ndata)
        % Computing the generalized Hessian product V * dX
        
        if Nalpha == 0
            B = 0;
            return
        elseif Nalpha == Ndata
            B = x;
            return
        end
        
        if Nalpha <= Ngamma
            ind1 = 1:Nalpha;
            ind2 = Nalpha+1:Ndata;
            
            Wag = [ones(Nalpha,Nbeta), Wag];
                        
            % utilize block matrix multiplication since W is in the form of
            % [Eaa Wa(b+g);W(b+g)a W(b+g)(b+g)] and is a sparse matrix
            % result matrix in the form of [F H; G J]
            Ql = Q(:, ind1);
            Qr = Q(:, ind2);
            
            FG = x * Ql;
            H = FG' * Qr;
            WH = Wag .* H; % G = H' since x is symmetric and Q is orthogonal;            
            WF = Ql' * FG;
            
            FH = WF * Ql' + WH * Qr';
            GJ = WH' * Ql';
            
            B = Ql * FH + Qr * GJ;
            
        else
            % utilize eigenvector sign ambiguity, rotate the input matrix by 180 degrees
            ind1 = 1:Nalpha+Nbeta;
            ind2 = Nalpha+Nbeta+1:Ndata;
            
            Wag =  [ones(Nbeta,Ngamma); 1-Wag];
            
            Ql = Q(:, ind1);
            Qr = Q(:, ind2);
            
            HJ = x * Qr;
            G = HJ' * Ql;
            WG = Wag' .* G; % G = H' since x is symmetric and Q is orthogonal;            
            WJ = Qr' * HJ;
            
            GJ = WJ * Qr' + WG * Ql';
            FH = WG' * Qr';

            B = x - (Ql * FH + Qr * GJ);
            
        end


    end
end

