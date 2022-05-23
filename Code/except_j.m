function result=except_j(j,n,x,v,w)
%x,u,w is for single view only
for l=setdiff(1:n,j)
    
    A_a(l)=norm(v(l,:)-v(j,:),2)*(x(:,j)'*x(:,l))*w(l,l);
end
result=sum(A_a);
end