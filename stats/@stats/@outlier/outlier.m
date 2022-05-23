function h = outlier(name, dataset, yl, yh, yle, yhe)
%

% Copyright 2003-2020 The MathWorks, Inc.

h = stats.outlier;

if nargin==0
   h.YLow='';
   h.YHigh='';
   h.YLowLessEqual=0;
   h.YHighGreaterEqual=0;
   h.name = '';
   h.dataset = '';
else
   h.YLow=yl;
   h.YHigh=yh;
   h.YLowLessEqual=yle;
   h.YHighGreaterEqual=yhe;
   h.dataset=dataset;
   
   % assumes name is unique
   h.name = name;
end

% add it to the list of outliers
s = settings;
if (s.stats.DistributionFitter.LegacyDistributionFitter.ActiveValue == false) 
    connect(h,stats.internal.dfit.getoutlierdb,'up');
else
    connect(h,dfswitchyard('getoutlierdb'),'up');
end

% Need to update session now, if name changes, and if object deleted
list(2) = handle.listener(h,findprop(h,'name'),'PropertyPostSet',@markdirty);
list(1) = handle.listener(h,'ObjectBeingDestroyed', @markdirty);
h.listeners=list;

markdirty;

%=============================================================================
function markdirty(varargin)
dfgetset('dirty',true);   % session has changed since last save
