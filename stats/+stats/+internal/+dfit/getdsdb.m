function thedsdb=getdsdb(varargin)
%GETDSDB returns the distributionFitter dataset database.


%   Copyright 2019-2020 The MathWorks, Inc.

thedsdb = dfgetset('thedsdb');

% Create a singleton class instance
if isempty(thedsdb)
   thedsdb = stats.dsdb;
   dfgetset('thedsdb',thedsdb);
end


