clear all;
close all;
clc
%%
rng(2020,'twister');
addpath('tools/');
addpath('print/');
%addpath('stats/');
options = [];
options.nRepeat = 50;
options.minIter = 50;
options.meanFitRatio = 0.1; % used for perview calcuObject
options.kmeans = 2; % 0:max; 1:kmeans; 2:spectral; 3:kmedoids
%options.WeightMode='HeatKernel';
%options.WeightMode='Binary';
options.WeightMode='Cosine';
options.maxIter = 30; % used in PerViewNMF
options.error = 9.2e-8; 
options.Gaplpha=0;  % for initialization of U, V in GMultiNMF
options.beta=0.001; % beta used in the paper for beta_s; 
options.rounds = 30; 
options.p=5;
rep=20; % # of replications of results

% Comment: p=5, beta near 0 using nnmf(CNA,99) gives the best
% 

%% read dataset
CNA = readmatrix('LIHC_preprocessed_CNA.csv');
DNA = readmatrix('LIHC_preprocessed_DNAMeth.csv');
RNA = readmatrix('LIHC_preprocessed_RNASeq.csv');
%%
[cna_var,cna_i]=sort(var(CNA),'descend');
[dna_var,dna_i]=sort(var(DNA),'descend');
[rna_var,rna_i]=sort(var(RNA),'descend');
%%
CNA=CNA(:,cna_i(1:100));
DNA=DNA(:,dna_i(1:100));
RNA=RNA(:,rna_i(1:100));
% [CNA,~] = nnmf(CNA,99);
% [DNA,~] = nnmf(DNA,99);
% [RNA,~] = nnmf(RNA,99);
%%

data{1}=CNA';
data{2}=DNA';
data{3}=RNA';
gnd=readmatrix('LIHC_labels');
gnd=gnd+1;
K=2;

% tiledlayout(3,2)
% nexttile
% Y1 = tsne(data{1}');
% gscatter(Y1(:,1),Y1(:,2),gnd)
% title('View 1')
% nexttile
% Y2 = tsne(data{2}');
% gscatter(Y2(:,1),Y2(:,2),gnd)
% title('View 2')
% nexttile
% Y3 = tsne(data{3}');
% gscatter(Y3(:,1),Y3(:,2),gnd)
% title('View 3')

% data=horzcat(data{1}',data{2}',data{3}');
% %data=data{3}';
% 
% for i=1:rep
% % [W,H]=nnmf(data,K);
% % idx = kmeans(W,K);
% idx = kmeans(data,K);
% [acc(i), nmi(i), Pi(i), Ri(i), Fi(i), ARi(i)] = Measure(gnd, idx);
% end
% %%
% LIHC_Ck=[acc',nmi',Pi',Ri',Fi',ARi'];
% save('LIHC_Ck','LIHC_Ck')
%%
%gscatter(centroidV(:,1),centroidV(:,2),gnd)
%% normalize data matrix
nv=length(data);   % Number of views

for i = 1:nv 
    %W{i}=constructW_cai(data{i}',options); % matrix A in the paper             
    X{i} = data{i} / sum(sum(data{i})); % call the normalized as X
    W{i}=constructW_cai(X{i}',options); % matrix A in the paper 
end

%%
% U = cell(1,nv);
% V = cell(1,nv);
% wt = cell(1,nv);
for i = 1:rep
    i
   [t,U{i},V{i},centroidV,wt{i},alphas,acc(i), nmi(i), Pi(i), Ri(i), Fi(i), ARi(i),l]=GMultiNMF(X, K, W,gnd, options);
end

%% result
[mean(acc) std(acc)]
[mean(nmi) std(nmi)]
[mean(Pi) std(Pi)]
[mean(Ri) std(Ri)]
[mean(Fi) std(Fi)]
[mean(ARi) std(ARi)]
LIHC_result=[acc',nmi',Pi',Ri',Fi',ARi'];
%save LIHC500.mat


