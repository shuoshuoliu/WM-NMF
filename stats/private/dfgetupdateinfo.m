function [name, distname, dataset, exrule, results] = dfgetupdateinfo(fit)
%DFGETUPDATEINFO GUI helper to delete an exclusion rule


%   Copyright 2003-2004 The MathWorks, Inc.

name = fit.name;
fittype = fit.fittype;
      
if strcmp(fittype, 'smooth')
    distname = 'nonparametric';
else %parametric
    distname = fit.distname;
end

dataset = fit.dataset;
exrule = fit.exclusionrulename;
results = fit.resultstext;
    







