function m = geomean(x,dim,varargin)
%GEOMEAN Geometric mean.
%   M = GEOMEAN(X) returns the geometric mean of the values in X.  When X
%   is an n element vector, M is the n-th root of the product of the n
%   elements in X.  For a matrix input, M is a row vector containing the
%   geometric mean of each column of X.  For N-D arrays, GEOMEAN operates
%   along the first non-singleton dimension.
%
%   GEOMEAN(X,'all') is the geometric mean of all the elements of X.
%
%   GEOMEAN(X,DIM) takes the geometric mean along dimension DIM of X.
%
%   GEOMEAN(X,VECDIM) finds the geometric mean of the elements of X based 
%   on the dimensions specified in the vector VECDIM.
%
%   GEOMEAN(...,NANFLAG) specifies how NaN values are treated. 
%     'includenan' - Include NaN values when computing the geometric mean, 
%                    resulting in NaN (the default). 
%     'omitnan'    - Ignore all NaN values in the input.
%
%   See also MEAN, HARMMEAN, TRIMMEAN.

%   Copyright 1993-2018 The MathWorks, Inc.


if any(x(:) < 0) || ~isreal(x)
    error(message('stats:geomean:BadData'))
end

if nargin < 2 || isempty(dim)
    % Figure out which dimension sum will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Take the n-th root of the product of elements of X, along dimension DIM.
if nargin < 2
    m = exp(mean(log(x)));
else
    m = exp(mean(log(x),dim,varargin{:}));
end
