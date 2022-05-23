function result = wtn(n, Y)
% Y is a cell of Ms*N matrix 
% returns the term in summation in the curly bracket on the denominator

viewNum = length(Y);
for i=1:viewNum
    temp=Y{i};
    sum_inv(i)=1/sum(temp(:,n).^2);
end
result=sum(sum_inv);
end