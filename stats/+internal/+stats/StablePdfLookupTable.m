function StablePdfLookupTable
%STABLEPDFLOOKUPTABLE Make density table for standard stable distribution.

ngrids = [51, 50, 500]; % Number of grids of BETA, ALPHA and X

a = linspace(0.4,2,ngrids(2));
b = linspace(-1,1,ngrids(1));
xgd = linspace(-0.999*pi/2,0.999*pi/2,ngrids(3));

p = zeros(ngrids);
xgrid = tan(xgd);
for k = 1:ngrids(1)
    for j = 1:ngrids(2)
        %fprintf('Calculating pdf for alpha = %4.8f and beta = %4.8f\n',a(j),b(k));
        p(k,j,:) = pdf('stable',xgrid,a(j),b(k),1,0);
    end
end

save('StablePdfTable','a','b','xgd','p');
end