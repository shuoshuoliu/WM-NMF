function ok = canplotdata(ds,dffig)
%CANPLOTDATA Determine if we can plot data in the current plot


%   Copyright 2003-2020 The MathWorks, Inc.

% No problem unless we have non-positive data
ok = true;
xlimits = ds.datalim;
xmin = xlimits(1);
if xmin>0
   return
end

% Get handle to control containing the distribution list
if nargin<2
   dffig = dfgetset('dffig');
end
hsel = getappdata(dffig,'selectioncontrols');

% No problem unless we have a probability plot
h = hsel(3);   % handle of display type control
choice = get(h,'Value');
ftypes = getappdata(h,'codenames');
ftype = ftypes{choice};
if ~isequal(ftype, 'probplot')
   return
end

% No problem unless distribution has support that excludes some data
ax = findobj(dffig, 'Tag', 'dfMainAxes');
distspec = getappdata(ax,'DistSpec');
if ~isempty(distspec)
    lobnd = distspec.support(1);
    strict = ~distspec.closedbound(1);
    if strict && lobnd>=xmin
        ok = false;
    elseif ~strict && lobnd>xmin
        ok = false;
    end
end
