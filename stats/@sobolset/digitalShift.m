function X = digitalShift(obj, X)
%DIGITALSHIFT Apply a digital shift to points.
%   DIGITALSHIFT(P,X) applies a digital shift to the points in X if there
%   is one defined in the point set P.

%   Copyright 2007-2018 The MathWorks, Inc.


DS = obj.DigitalShifts;
if ~isempty(DS)
    X = internal.stats.sobolset.matrixBitXor(X, DS);
end
