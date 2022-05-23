clear all;
close all;
clc
%%
rng(2020,'twister');
addpath('tools/');
addpath('print/');
addpath('stats/');
options = [];
options.nRepeat = 10;
options.minIter = 30;
options.meanFitRatio = 0.1; % used for perview calcuObject
options.Gaplpha=1;  % for initialization of U, V in GMultiNMF

options.kmeans = 2; % 0:max; 1:kmeans; 2:spectral; 3:kmedoids
%options.WeightMode='HeatKernel';
%options.WeightMode='Binary';
options.WeightMode='Cosine';
options.maxIter = 30; % used in PerViewNMF
options.error = 9.2e-8;

options.beta=0.01; % beta used in the paper for graph.n; 
options.rounds = 2; 
options.p=5;
rep=20; % # of replications of results

%% read dataset
load handwritten.mat
data{1}=fourier';
data{2}=pixel';
data{3}=zer';
data{4}=profile';
K = 10;
gnd=gnd+1;


%% normalize data matrix
nv=length(data);   % Number of views

for i = 1:nv             
    XX{i} = data{i} / sum(sum(data{i})); % call the normalized as X
    W{i}=constructW_cai(XX{i}',options); % matrix A in the paper 
end

%%
for i = 1:rep
    i
   [t,U{i},V{i},centroidV,wt{i},alphas,acc(i), nmi(i), Pi(i), Ri(i), Fi(i), ARi(i),l]=GMultiNMF(XX, K, W,gnd, options);
end

%% result
acc_sum=[mean(acc) std(acc)];
nmi_sum=[mean(nmi) std(nmi)];
Pi_sum=[mean(Pi) std(Pi)];
Ri_sum=[mean(Ri) std(Ri)];
Fi_sum=[mean(Fi) std(Fi)];
ARi_sum=[mean(ARi) std(ARi)];
