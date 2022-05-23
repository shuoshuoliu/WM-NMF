function K = poly_kernel(u,v,polyOrder,varargin)
%POLY_KERNEL Polynomial kernel for SVM functions

% Copyright 2004-2012 The MathWorks, Inc.


if nargin < 3 || isempty(polyOrder)
    polyOrder = 3; %default order
else
     if ~isscalar(polyOrder) || ~isnumeric(polyOrder)
        error(message('stats:poly_kernel:BadPolyOrder'));
    end
    if polyOrder ~= floor(polyOrder) || polyOrder < 1
        error(message('stats:poly_kernel:PolyOrderNotInt'))
    end
end

dotproduct = (u*v');

K = dotproduct;

for i = 2:polyOrder
    K = K.*(1 + dotproduct);
end
