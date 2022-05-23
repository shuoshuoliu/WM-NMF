function [D, C, F] = dfviewdata(dataset)
% DFVIEWDATA Helper function for Distribution Fitter viewdata panel
%
%    [X, Y, W] = DFVIEWDATA(DATASET)
%    returns the x, y and w values for the given dataset
%	 (in a manner that the Java GUI can use)

%   Copyright 2003-2016 The MathWorks, Inc. 


ds = handle(dataset);
D = ds.y;
C = ds.censored;
F = ds.frequency;
