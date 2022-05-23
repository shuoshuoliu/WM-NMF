function Q = calculateQ(U)
    %mFea = size(U,1);
    norms = sum(abs(U),1);
    norms = max(norms,1e-10);
    Q=diag(norms);
    %Q=repmat(norms,mFea,1);
end