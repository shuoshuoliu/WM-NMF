clear all;
close all;
clc
%%
rng(2022,'twister');
addpath('tools/');
addpath('print/');
addpath('stats/');

options = [];
% options.nRepeat = 10;
% options.meanFitRatio = 0.1; % used for perview calcuObject
% options.alpha=100;  % for initialization of U, V in GMultiNMF
% Goptions.alpha=100;

options.WeightMode='HeatKernel';
%options.WeightMode='Binary';
%options.WeightMode='Cosine';

%% digit dataset
% load handwritten.mat
% data{1}=fourier';
% data{2}=pixel';
% data{3}=zer';
% data{4}=profile';

% sync6 dataset
%load data6.mat
K = 20; % the range of possible K

% LIHC
CNA = readmatrix('LIHC_preprocessed_CNA.csv');
DNA = readmatrix('LIHC_preprocessed_DNAMeth.csv');
RNA = readmatrix('LIHC_preprocessed_RNASeq.csv');
%
[cna_var,cna_i]=sort(var(CNA),'descend');
[dna_var,dna_i]=sort(var(DNA),'descend');
[rna_var,rna_i]=sort(var(RNA),'descend');
CNA=CNA(:,cna_i(1:100));
DNA=DNA(:,dna_i(1:100));
RNA=RNA(:,rna_i(1:100));
data{1}=CNA';
data{2}=DNA';
data{3}=RNA';

%% normalize data matrix
%dt = data{3} / sum(sum(data{3})); % call the normalized as X 
dt = data{3}; % call the normalized as X 
%%
for i=2:K
    for j=1:K
        A=constructW_cai(dt,options); % matrix A in the paper
        [U, V] = GNMF(dt', i, A, options);
        %D = diag(sum(A));
        %L = D - A; % Get matrix L
        error(j)=norm(V*U'-dt);
        %error(j)=trace(V'*L*V);
    end
err3(i-1)=mean(error);
end

%%
figure
plot(2:K,smooth(err3))
title('LIHC Dataset')
set(gca,'XTick',2:20, 'YTick',[])
exportgraphics(gcf,'LIHC.png','Resolution',300)



