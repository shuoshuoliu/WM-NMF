function K = quadratic_kernel(u,v,varargin)
%QUADRATIC_KERNEL Quadratic kernel for SVM functions

% Copyright 2004-2012 The MathWorks, Inc.


dotproduct = (u*v');
K = dotproduct.*(1 + dotproduct);
