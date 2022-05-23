function toggleaxlimctrl(dffig)
%DFTOGGLEAXLIMCTRL Toggle x and y axis limit controls on or off


%   Copyright 2003-2004 The MathWorks, Inc.

% Get handle to menu item, may be current object or may not
h = gcbo;
if ~isequal(get(h,'Tag'),'showaxlimctrl')
   h = findall(dffig,'Tag','showaxlimctrl');
end

% Get new state
onoff = stats.internal.dfit.on2off(get(h,'Checked'));
dfgetset('showaxlimctrl',onoff);

% Add or remove controls
stats.internal.dfit.axlimctrl(dffig,onoff)

% Change menu state
set(h,'Checked',onoff);

