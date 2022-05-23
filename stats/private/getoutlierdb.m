function theoutlierdb=getoutlierdb(varargin)


%   Copyright 2003-2004 The MathWorks, Inc.

theoutlierdb = dfgetset('theoutlierdb');

% Create a singleton class instance
if isempty(theoutlierdb)
   theoutlierdb = stats.outlierdb;
   dfgetset('theoutlierdb',theoutlierdb);
end
