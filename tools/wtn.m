function fea = wtn(Y)
% Y is an Ms*N matrix 
for i=1:viewNum
        Y{i}=X{i}-U{i}*V{i}'; 
        temp=Y{i};
        w_vec=zeros(nSmp,1);
        for j=1:nSmp % n-th observation
            Ysum=sum(temp(:,j).^2);
            w_vec(j)=;
            wt{i}=diag(w_vec);
        end
    end
end