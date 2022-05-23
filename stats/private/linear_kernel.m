function K = linear_kernel(u,v,varargin) 
%LINEAR_KERNEL Linear kernel for SVM functions

% Copyright 2004-2012 The MathWorks, Inc.

K = (u*v');
