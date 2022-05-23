function m = harmmean(x,varargin)
%HARMMEAN Harmonic mean.
%   M = HARMMEAN(X) returns the harmonic mean of the values in X.  For
%   vector input, M is the inverse of the mean of the inverses of the
%   elements in X.  For matrix input, M is a row vector containing the
%   harmonic mean of each column of X.  For N-D arrays, HARMMEAN operates
%   along the first non-singleton dimension.
%
%   HARMMEAN(X,'all') is the harmonic mean of all the elements in X.
%
%   HARMMEAN(X,DIM) takes the harmonic mean along dimension DIM of X.
%
%   HARMMEAN(X,VECDIM) finds the harmonic mean of the elements of X based 
%   on the dimensions specified in the vector VECDIM.
%   HARMMEAN(...,NANFLAG) specifies how NaN values are treated. 
%     'includenan' - Include NaN values when computing the harmonic mean, 
%                    resulting in NaN (the default). 
%     'omitnan'    - Ignore all NaN values in the input.
%
%   See also MEAN, GEOMEAN, TRIMMEAN.

%   Copyright 1993-2018 The MathWorks, Inc. 

% Take the reciprocal of the mean of the reciprocals of X.
if nargin == 1
    m = 1 ./ mean(1./x);
else
    m = 1 ./ mean(1./x, varargin{:});
end
