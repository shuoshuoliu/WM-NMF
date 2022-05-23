function y = iqr(x,dim)
%IQR Interquartile range. 
%   Y = IQR(X) returns the interquartile range of the values in X.  For
%   vector input, Y is the difference between the 75th and 25th percentiles
%   of X.  For matrix input, Y is a row vector containing the interquartile
%   range of each column of X.  For N-D arrays, IQR operates along the
%   first non-singleton dimension.
%
%   The IQR is a robust estimate of the spread of the data, since changes
%   in the upper and lower 25% of the data do not affect it.
%
%   IQR(X,'all') calculates the interquartile of all the elements in X.
%
%   IQR(X,DIM) calculates the interquartile range along the dimension DIM
%   of X.
%
%   IQR(X,VECDIM) calculates the interquartile of elements of X based on 
%   the dimensions specified in the vector VECDIM.
%
%   See also PRCTILE, STD, VAR.

%   Copyright 1993-2018 The MathWorks, Inc. 


if nargin == 1
    y = diff(prctile(x, [25; 75]));
else
    t = prctile(x, [25; 75],dim);
    if ((ischar(dim) && isrow(dim)) || ...
     (isstring(dim) && isscalar(dim) && (strlength(dim) > 0))) && ...
     strncmpi(dim,'all',max(strlength(dim), 1))
        dim = 1;
    end
    y = diff(t,[],min(dim));
end