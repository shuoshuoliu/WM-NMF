function kval = mlp_kernel(u,v,P1,P2,varargin)
%MLP_KERNEL Multilayer perceptron kernel for SVM functions

% Copyright 2004-2012 The MathWorks, Inc.


if nargin < 3 || isempty(P1)
    P1 = 1;
end

if nargin < 4 ||  isempty(P2)
    P2 = -1;
end

kval = tanh(P1*(u*v')+P2);
