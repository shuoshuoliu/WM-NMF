function setbinwidthrules(dataset, v0, v1, v2, v3, v4, v5, v6)
% DFSETBINWIDTHRULES Helper function for the dfittool set default bin width panel

%   Copyright 2003-2020 The MathWorks, Inc. 


binDlgInfo.rule = v0;
binDlgInfo.nbinsExpr = v1;
binDlgInfo.nbins = str2num(v1); %#ok<*ST2NM>
binDlgInfo.widthExpr = v2;
binDlgInfo.width = str2num(v2);
binDlgInfo.placementRule = v3;
binDlgInfo.anchorExpr = v4;
binDlgInfo.anchor = str2num(v4);
binDlgInfo.applyToAll = v5;
binDlgInfo.setDefault = v6;


if (v5 == true) % apply to all
    dsdb = stats.internal.dfit.getdsdb;
    ds = down(dsdb);
    while(~isempty(ds))
        ds.binDlgInfo = binDlgInfo;
        ds = right(ds);
    end
    stats.internal.dfit.updateallplots(true, false, true); 
elseif ~isempty(dataset)
    ds = handle(dataset);
    ds.binDlgInfo = binDlgInfo;
    clearplot(ds);
    updateplot(ds);
end

% Update the main fig axis limits to fit the new histograms
stats.internal.dfit.updatexlim;
stats.internal.dfit.updateylim;

if (v6 == true) % set default
    dfgetset('binDlgInfo', binDlgInfo);
end

dfgetset('dirty',true);   % session has changed since last save
