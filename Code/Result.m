function [acc, nmi, Pi, Ri, Fi, ARi] = Result(X, gnd, K, kmeansFlag)

if kmeansFlag == 1
    label = litekmeans(X, K, 'Replicates',20);
elseif kmeansFlag == 2
    label=spectralcluster(X,K);
elseif kmeansFlag == 3
    label=kmedoids(X,K);
else
    [~, label] = max(X, [] ,2);
end

[acc, nmi, Pi, Ri, Fi, ARi] = Measure(gnd, label);
