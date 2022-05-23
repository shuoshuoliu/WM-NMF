clear all;
close all;
clc
%%
rng(2020,'twister');
addpath('simulation/');
addpath('tools/');
addpath('print/');
options = [];
options.nRepeat = 10;
options.minIter = 30;
options.meanFitRatio = 0.1; % used for perview calcuObject
options.Gaplpha=1;  % for initialization of U, V in GMultiNMF

options.kmeans = 1; % 0:max; 1:kmeans; 2:spectral; 3:kmedoids
options.WeightMode='HeatKernel';
%options.WeightMode='Binary';
%options.WeightMode='Cosine';
options.maxIter = 30; % used in PerViewNMF
options.error = 9.2e-8;
% Comment: choose a little bit larger beta if p is small (range:0.1-2)
%options.beta=0.01; % beta used in the paper for graph.n; 
options.rounds = 2; % set 2 for normal data and 3 for noisy data 
options.p=5;
rep=8; % # of replications of results

%%
load data6.mat
load gnd.mat
gnd=double(gnd);
K=10;
% 
% data=horzcat(data{1}',data{2}',data{3}',data{4}',data{5}',data{6}');
% %data=data{4}';
% 
% for i=1:20
%     i
% idx = spectralcluster(data,K);
% %idx = kmeans(data,K);
% [acc(i), nmi(i), Pi(i), Ri(i), Fi(i), ARi(i)] = Measure(gnd, idx);
% end
% %
% [mean(acc) std(acc)]
% [mean(nmi) std(nmi)]
% [mean(Pi) std(Pi)]
% [mean(Ri) std(Ri)]
% [mean(Fi) std(Fi)]
% [mean(ARi) std(ARi)]
%%
% tiledlayout(2,3)
% nexttile
% Y1 = tsne(data{1,1}');
% gscatter(Y1(:,1),Y1(:,2),gnd)
% title('View 1')
% nexttile
% Y2 = tsne(data{1,2}');
% gscatter(Y2(:,1),Y2(:,2),gnd)
% title('View 2')
% nexttile
% Y3 = tsne(data{1,3}');
% gscatter(Y3(:,1),Y3(:,2),gnd)
% title('View 3')
% nexttile
% Y4 = tsne(data{1,4}');
% gscatter(Y4(:,1),Y4(:,2),gnd)
% title('View 4')
% nexttile
% Y5 = tsne(data{1,5}');
% gscatter(Y5(:,1),Y5(:,2),gnd)
% title('View 5')
% nexttile
% Y6 = tsne(data{1,6}');
% gscatter(Y6(:,1),Y6(:,2),gnd)
% title('View 6')
%% normalize data matrix
nv=length(data);   % Number of views

for i = 1:nv             
    X{i} = data{i} / sum(sum(data{i})); % call the normalized as X
    W{i}=constructW_cai(X{i}',options); % matrix A in the paper 
end

%%
% U = cell(1,nv);
% V = cell(1,nv);
% wt = cell(1,nv);
res=zeros(5,6);
for j = 1:5
   options.beta=j*0.005;
   for i=1:rep
       [t,U{i},V{i},centroidV,wt{i},alphas,acc(i), nmi(i), Pi(i), Ri(i), Fi(i), ARi(i),l]=GMultiNMF(X, K, W,gnd, options);
   end
   acc_sum=mean(acc);
    nmi_sum=mean(nmi);
    Pi_sum=mean(Pi);
    Ri_sum=mean(Ri);
    Fi_sum=mean(Fi);
    ARi_sum=mean(ARi);
    res(j,:)=[acc_sum, nmi_sum, Pi_sum, Ri_sum, Fi_sum, ARi_sum]
end
save beta_syn.mat
% 
%% result
% [mean(acc) std(acc)]
% [mean(nmi) std(nmi)]
% [mean(Pi) std(Pi)]
% [mean(Ri) std(Ri)]
% [mean(Fi) std(Fi)]
% [mean(ARi) std(ARi)]
% %%
%save synthetic6_nono.mat

















