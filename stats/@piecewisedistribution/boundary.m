function [p,q] = boundary(obj,j)
%BOUNDARY Boundary between segments of piecewise distribution.
%    [P,Q]=BOUNDARY(OBJ) returns the boundary points between the segments
%    of the piecewise distribution defined by OBJ.  P is a vector of the
%    probability values at each boundary.  Q is a vector of the quantile
%    values at each boundary.
%
%    [P,Q]=BOUNDARY(OBJ,J) returns P and Q for the Jth boundary.
%
%    See also PIECEWISEDISTRIBUTION, PIECEWISEDISTRIBUTION/NSEGMENTS.

%   Copyright 2006-2007 The MathWorks, Inc. 


if nargin==1
    p = obj.P;
    q = obj.Q;
else
    if ~isvector(j)
        j = j(:);
    end
    if any(~ismember(j,1:numel(obj.P)))
        error(message('stats:piecewisedistribution:boundary:BadIndex'));
    end
    p = obj.P(j);
    q = obj.Q(j);
end
