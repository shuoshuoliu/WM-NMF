function togglelegend(dffig,onoff,leginfo)
%TOGGLELEGEND Toggle distribution fitter legend on or off


%   Copyright 2003-2007 The MathWorks, Inc.

% Get new state if not passed in -- note uimenu state reflects
% old state, and uitoggletool state reflects new state
if nargin<2 || isempty(onoff)
   h = gcbo;
   if isequal(get(h,'Type'),'uimenu')
      onoff = stats.internal.dfit.on2off(get(h,'Checked'));
   else
      onoff = get(h,'State');
   end
end
dfgetset('showlegend',onoff);

if nargin<3
    leginfo = {};
end

% Change menu state
h = findall(dffig,'Type','uimenu','Tag','showlegend');
if ~isempty(h), set(h,'Checked',onoff); end

% Change button state
h = findall(dffig,'Type','uitoggletool','Tag','showlegend');
if ~isempty(h), set(h,'State',onoff); end

% Change legend state
stats.internal.dfit.updatelegend(dffig, false, leginfo);
