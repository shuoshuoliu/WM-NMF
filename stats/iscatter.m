function h = iscatter(x,y,i,c,m,msize)
%ISCATTER Scatter plot grouped by index vector.
%   ISCATTER(X,Y,I,C,M,msize) displays a scatter plot of X vs. Y grouped
%   by the index vector I.  
%
%   No error checking.  Use GSCATTER instead.
%
%   See also GSCATTER, GPLOTMATRIX.

%   Copyright 1993-2014 The MathWorks, Inc.

% All functionality is provided by internal.stats.scatterplot
if nargin > 3
    c = convertStringsToChars(c);
end

if nargin > 4
    m = convertStringsToChars(m);
end

hh = internal.stats.scatterplot(gca, x, y, i ,c ,m, msize);

% Return the handles if desired.  They are arranged so that even if X
% or Y is a matrix, the first ni elements of hh(:) represent different
% groups, so they are suitable for use in drawing a legend.
if (nargout>0)
   h = hh(:);
end
