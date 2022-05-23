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
% Comment: choose a little bit larger beta if p is small (range:0.1-2)
options.beta=0.01; % beta used in the paper for graph.n; 
options.rounds = 2; % set 2 for normal data and 3 for noisy data 
options.p=5;
rep=3; % # of replications of results

%% read dataset
load handwritten.mat
data{1}=fourier';
data{2}=pixel';
data{3}=zer';
data{4}=profile';
K = 10;
gnd=gnd+1;

% load handwritten.mat
% data{1}=fourier';
% data{2} = pixel';
% aa = zer';
% aa(:,1:300)=normrnd(300,30,size(aa(:,1:300)));
% data{3}=aa;
% K = 10;
% gnd=gnd+1;

% load MSRC.mat;
% data{1}=X{1,1}';
% data{2}=X{1,2}';
% data{3}=X{1,3}';
% data{4}=X{1,4}';
% % data{5}=X{1,5}';
% K=7;
% gnd=Y;

% data=horzcat(data{1}',data{2}',data{3}',data{4}',data{5}');
% %data=data{5}';
% 
% for i=1:rep
% [W,H]=nnmf(data,K);
% idx = kmeans(W,K);
% % idx = kmeans(data,K);
% [acc(i), nmi(i), Pi(i), Ri(i), Fi(i), ARi(i)] = Measure(gnd, idx);
% end
% 
% [mean(acc) std(acc)]
% [mean(nmi) std(nmi)]
% [mean(Pi) std(Pi)]
% [mean(Ri) std(Ri)]
% [mean(Fi) std(Fi)]
% [mean(ARi) std(ARi)]
%% normalize data matrix
nv=length(data);   % Number of views

for i = 1:nv 
    %W{i}=constructW_cai(data{i}',options); % matrix A in the paper             
    XX{i} = data{i} / sum(sum(data{i})); % call the normalized as X
    W{i}=constructW_cai(XX{i}',options); % matrix A in the paper 
end

%%
% U = cell(1,nv);
% V = cell(1,nv);
% wt = cell(1,nv);
for i = 1:rep
    i
   [t,U{i},V{i},centroidV,wt{i},alphas,acc(i), nmi(i), Pi(i), Ri(i), Fi(i), ARi(i),l]=GMultiNMF(XX, K, W,gnd, options);
end
l;

%% result
acc_sum=[mean(acc) std(acc)];
nmi_sum=[mean(nmi) std(nmi)];
Pi_sum=[mean(Pi) std(Pi)];
Ri_sum=[mean(Ri) std(Ri)];
Fi_sum=[mean(Fi) std(Fi)];
ARi_sum=[mean(ARi) std(ARi)];
%save digit.mat
%%
%y = smooth(l) ;
% plot([2,3,4,5,6],time,'r','LineWidth',3)
% xlabel('Views #')
% ylabel('Time(s)')

%%
% k=X{1}*wt{1,1}{1,1};
% kk=X{2}*wt{1,1}{1,2};
% kkk=X{3}*wt{1,1}{1,3};
% kkkk=X{4}*wt{1,1}{1,4};
% %original data
% k=k(1:30,1:30);
% kk=kk(1:30,1:30);
% kkk=kkk(1:30,1:30);
% kkkk=kkkk(1:30,1:30);
% 
% %estimated data by our method
% z=U{1,1}{1,1}*V{1,1}{1,1}'*wt{1,1}{1,1};
% zz=U{1,1}{1,2}*V{1,1}{1,2}'*wt{1,1}{1,2};
% zzz=U{1,1}{1,3}*V{1,1}{1,3}'*wt{1,1}{1,3};
% zzzz=U{1,1}{1,4}*V{1,1}{1,4}'*wt{1,1}{1,4};
% z=z(1:30,1:30);
% zz=zz(1:30,1:30);
% zzz=zzz(1:30,1:30);
% zzzz=zzzz(1:30,1:30);
% 
% t = tiledlayout(1,4); % Requires R2019b or later
% nexttile
% % surf(k)
% % hold on
% imagesc(k)
% nexttile
% % surf(kk)
% % hold on
% imagesc(kk)
% nexttile
% % surf(kkk)
% % hold on
% imagesc(kkk)
% nexttile
% % surf(kkkk)
% % hold on
% imagesc(kkkk)
% 
% %%
% tt = tiledlayout(1,4);
% nexttile
% % surf(z)
% % hold on
% imagesc(z)
% nexttile
% % surf(zz)
% % hold on
% imagesc(zz)
% nexttile
% % surf(zzz)
% % hold on
% imagesc(zzz)
% nexttile
% % surf(zzzz)
% % hold on
% imagesc(zzzz)
% 
% %%

% %%
% a=diag(wt{1,1}{1,1});
% b=diag(wt{1,1}{1,2});
% c=diag(wt{1,1}{1,3});
% 
% g = tiledlayout(2,2);
% nexttile
% plot(1:2000,a);
% ylim([0 1])
% title('Weights for view 1')
% nexttile
% plot(1:2000,b);
% ylim([0 1])
% title('Weights for view 2')
% nexttile
% plot(1:300,c(1:300,:),'r')
% title('First 300 weights for view 3')
% nexttile
% plot(1:2000,c);
% ylim([0 1])
% title('Weights for view 3')
% 
% 
