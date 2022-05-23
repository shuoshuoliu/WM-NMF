function [ahat,bhat,aci,bci] = unifit(x,alpha)
%UNIFIT Parameter estimates for uniformly distributed data.
%   UNIFIT(X,ALPHA) Returns the maximum likelihood estimates (MLEs) of the  
%   parameters of the uniform distribution given the data in the vector, X.  
%
%   [AHAT,BHAT,ACI,BCI] = UNIFIT(X,ALPHA) gives MLEs and 
%   100(1-ALPHA) percent confidence intervals given the data.
%   ALPHA is optional. By default ALPHA = 0.05 which corresponds
%   to 95% confidence intervals. 
%
%   See also UNIFCDF, UNIFINV, UNIFPDF, UNIFRND, UNIFSTAT, MLE. 

%   Copyright 1993-2019 The MathWorks, Inc. 


if nargin < 2 
    alpha = 0.05;
else
    alpha = double(alpha);
end
ahat = min(x);
bhat = max(x);

% sample size per min/max value is the total number of data values 
% divided by the number of min/max values
ssz = numel(x)/numel(ahat);

tmp = (bhat - ahat)./alpha.^(1./ssz);
aci = [bhat - tmp; ahat];
bci = [bhat; ahat + tmp];
