function [v0, v1, v2, v3, v4, v5, v6] = dfgetbinwidthdefaults()
% DFGETBINWIDTHDEFAULTS Helper function for the dfittool set default bin width panel

%   Copyright 2003-2004 The MathWorks, Inc. 


binDlgInfo = dfgetset('binDlgInfo');
v0 = binDlgInfo.rule;
v1 = binDlgInfo.nbinsExpr;
v2 = binDlgInfo.widthExpr;
v3 = binDlgInfo.placementRule;
v4 = binDlgInfo.anchorExpr;
v5 = binDlgInfo.applyToAll;
v6 = binDlgInfo.setDefault;
