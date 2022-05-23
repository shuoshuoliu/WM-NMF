function cbkclear
%CBKCLEAR Callback for Clear button


%   Copyright 2003-2004 The MathWorks, Inc.

% Clear all saved fits from the plot and notify fits manager
fitdb = stats.internal.dfit.getfitdb;
fit = down(fitdb);
while(~isempty(fit))
   fit.plot = 0;
   fit = right(fit);
end

% Clear all datasets from the plot and notify data sets manager
dsdb = stats.internal.dfit.getdsdb;
ds = down(dsdb);
while(~isempty(ds))
   ds.plot = 0;
   ds = right(ds);
end

stats.internal.dfit.updatexlim;
stats.internal.dfit.updateallplots;
stats.internal.dfit.updateylim;
