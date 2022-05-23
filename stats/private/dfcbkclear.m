function dfcbkclear
%DFCBKCLEAR Callback for Clear button


%   Copyright 2003-2004 The MathWorks, Inc.

% Clear all saved fits from the plot and notify fits manager
fitdb = getfitdb;
fit = down(fitdb);
while(~isempty(fit))
   fit.plot = 0;
   fit = right(fit);
end

% Clear all datasets from the plot and notify data sets manager
dsdb = getdsdb;
ds = down(dsdb);
while(~isempty(ds))
   ds.plot = 0;
   ds = right(ds);
end

dfupdatexlim;
dfupdateallplots;
dfupdateylim;
