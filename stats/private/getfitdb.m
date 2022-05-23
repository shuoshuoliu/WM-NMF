function fitdb=getfitdb(varargin)
% GETFITDB A helper function for DFITTOOL


% Copyright 2003-2004 The MathWorks, Inc.

thefitdb = dfgetset('thefitdb');

% Create a singleton class instance
if isempty(thefitdb)
   thefitdb = stats.fitdb;
end

dfgetset('thefitdb',thefitdb);
fitdb=thefitdb;
