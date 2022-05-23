function thedsdb=getdsdb(varargin)


%   Copyright 2003-2004 The MathWorks, Inc.

thedsdb = dfgetset('thedsdb');

% Create a singleton class instance
if isempty(thedsdb)
   thedsdb = stats.dsdb;
   dfgetset('thedsdb',thedsdb);
end


