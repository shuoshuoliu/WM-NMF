function theoutlierdb=getoutlierdb(varargin)
% GETOUTLIERDB returns the distributionFitter exclusion rule database.

%   Copyright 2019-2020 The MathWorks, Inc.

theoutlierdb = dfgetset('theoutlierdb');

% Create a singleton class instance
if isempty(theoutlierdb)
   theoutlierdb = stats.outlierdb;
   dfgetset('theoutlierdb',theoutlierdb);
end
