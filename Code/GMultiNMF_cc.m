function [t,U,V,centroidV,wt] = GMultiNMF_cc(X, K, W, options)
%	Notation:
% 	X ... a cell containing all views for the normalized data
%   data. original data
% 	K ... number of hidden factors
% 	W ... weight matrix of the affinity graph 
% 	gnd...ground truth labels
%Note that columns are data vectors here
viewNum = length(X);
Rounds = options.rounds;
nSmp=size(X{1},2); % Number of data points
p=options.p;

U_ = [];
V_ = [];

U = cell(1, viewNum);
V = cell(1, viewNum);
%% initialize U,V， wt
% wt=cell(1,viewNum);
% for i=1:viewNum
%     wt{i}=diag(ones(1,nSmp)/viewNum);
% end
%tic
j = 0;
while j < 7
    j = j + 1;
    Goptions.alpha=options.Gaplpha;
    if j == 1
        %rand('twister',5489);
        [U{1}, V{1}] = GNMF(X{1}, K, W{1}, Goptions);
    else
        %rand('twister',5489);
        [U{1}, V{1}] = GNMF(X{1}, K, W{1}, options, U_, V{viewNum});      
    end
    for i = 2:viewNum
        %rand('twister',5489);
        [U{i}, V{i}] = GNMF(X{i}, K, W{i},Goptions, U_, V{i-1});
    end   
end
%toc
alphas=ones(1,viewNum)/viewNum; % initialize alphas in the paper
%% initialize the wt, V*
centroidV = alphas(1)^p * V{1};
for i = 2:viewNum
    centroidV = centroidV + alphas(i)^p * V{i};
end
centroidV = centroidV / sum(alphas.^p);

for i=1:viewNum
    wt{i}=diag(ones(1,nSmp)/viewNum);
end
    
%% ==================== outer loop ================================
jjj = 0;

l=zeros(1,Rounds);

tic
while jjj < Rounds  % ================ outer loop ===============
    jjj = jjj + 1;
    % ====================== inner loop =========================
    optionsForPerViewNMF = options; % used for per view updates
    for i = 1:viewNum
        optionsForPerViewNMF.alpha = alphas(i);
        %rand('twister',5489);
        [U{i}, V{i}] = PerViewNMF(X{i}, K, centroidV, W{i} , wt{i}, optionsForPerViewNMF, U{i}, V{i});
    end
    % ===================== update alpha ========================
    A=zeros(viewNum,1);
    for i=1:viewNum
        dV= V{i}-centroidV;
        A(i)=sum(sum(dV.^2));
    end
   
    dleft =sum((1./(A)).^(1/(p-1)));
    
    for i=1:viewNum
        dright=(A(i))^(1/(p-1));
        alphas(i)=1/dleft/dright;
    end
    % ================= update V* =======================
    centroidV = alphas(1)^p * V{1};
    for i = 2:viewNum
        centroidV = centroidV + alphas(i)^p * V{i};
    end
    centroidV = centroidV / sum(alphas.^p);
    % ================ update wt ========================
    wt=cell(1,viewNum);
    for i=1:viewNum
        Y{i}=X{i}-U{i}*V{i}'; 
    end
    
    for j=1:viewNum
        wt{j}=diag(ones(nSmp,1));
        for k=1:nSmp
            result=wtn(k,Y);
            Ysum=sum(Y{j}(:,k).^2);
            wt{j}(k,k)=1/result/Ysum;
        end
    end
    % =================== Update A ==================
    for i=1:viewNum
        W{i}=constructW_cai(wt{i}*X{i}',options); % matrix A in the paper  
    end
    % ====================================================
    logL = 0; % Loss for the round
    for i = 1:viewNum
        alpha=alphas(i);
        Wtemp = 10*alpha*W{i}; % Modify the weight matrix with the involved parameters
        %Wtemp=W{i};
        DCol = full(sum(Wtemp,2));
        D = spdiags(DCol,0,nSmp,nSmp);
        L = D - Wtemp; % Get matrix L
        if isfield(options,'NormW') && options.NormW
            D_mhalf = spdiags(DCol.^-.5,0,nSmp,nSmp) ;
            L = D_mhalf*L*D_mhalf;
        end
        % Compute the loss
        tmp1 = (X{i} - U{i}*V{i}')*wt{i};
        tmp2 = (V{i} - centroidV);
        logL = logL+sum(sum(tmp1.^2))+alpha^p* (sum(sum(tmp2.^2)))+options.beta*trace(V{i}'*L*V{i});%sum(sum((V{i}'*L).*V{i}')); 
    end

    l(jjj)=logL;
    %rand('twister',5489);
    if jjj==2
        break
    end
end
t=toc;
end



