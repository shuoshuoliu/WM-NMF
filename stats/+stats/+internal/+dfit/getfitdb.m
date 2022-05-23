function fitdb=getfitdb(varargin)
% %GETFITDB returns the distributionFitter fit database.


% Copyright 2019-2020 The MathWorks, Inc.

thefitdb = dfgetset('thefitdb');

% Create a singleton class instance
if isempty(thefitdb)
   thefitdb = stats.fitdb;
end

dfgetset('thefitdb',thefitdb);
fitdb=thefitdb;
