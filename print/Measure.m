function [acc, nmi, Pi, Ri, Fi, ARi] = Measure(truth, idx)

acc = 1-compute_CE(idx, truth); % clustering accuracy
[Fi,Pi,Ri] = compute_f(truth,idx); % F1, precision, recall
nmi = compute_nmi(truth,idx);
ARi = rand_index(truth,idx);


