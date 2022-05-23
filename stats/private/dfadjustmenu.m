function dfadjustmenu(dffig)
%DFADJUSTMENU Adjust contents of curve fit plot menus


%   Copyright 2003-2011 The MathWorks, Inc.

% Remove some menus entirely
h = findall(dffig, 'Type','uimenu', 'Parent',dffig);
removelist = {'figMenuEdit' 'figMenuInsert' 'figMenuDesktop'};
for j=1:length(removelist)
   h0 = findall(h,'flat', 'Tag',removelist{j});
   if (~isempty(h0))
      delete(h0);
      h(h==h0) = [];
   end
end

% Add or remove some items from other menus
% Fix FILE menu
h0 = findall(h,'flat', 'Tag','figMenuFile');
h1 = findall(h0, 'Type','uimenu', 'Parent',h0);
for j=length(h1):-1:1
   mtag = get(h1(j),'Tag');
   if isequal(mtag,'figMenuFileClose')
      m7 = h1(j);
      set(m7,'Label',getString(message('stats:dfstrings:label_CloseDistributionFitting')),'Tag','dfitMenuClose')
   elseif isequal(mtag,'printMenu')
      m5 = h1(j);
   else
      delete(h1(j));
      h1(j) = [];
   end
end
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_ImportData')), 'Position',1,...
      'Callback','dfittool(''import data'')', 'Tag','dfitMenuImportData');
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_ClearSession')),'Position',2,...
       'Callback','dfittool(''clear session'')','Separator','on', ...
       'Tag','dfitMenuImportClearSession');
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_LoadSession')), 'Position',3,...
      'Callback','dfittool(''load session'')', 'Tag','dfitMenuLoadSession');
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_SaveSession')), 'Position',4,...
           'Callback','dfittool(''save session'')', 'Tag','dfitMenuSaveSession');
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_GenerateCode')), 'Position',5,...
           'Callback','dfittool(''generate code'')', 'Tag','dfitMenuGenCode');

uimenu(h0, 'Label',getString(message('stats:dfstrings:label_DefineCustomDistributions')),'Position',6,...
           'Callback',{@dfcustomdist,'define'}','Separator','on', ...
           'Tag','dfitMenuDefineCustom');
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_ImportCustomDistributions')), 'Position',7, ...
           'Callback',{@dfcustomdist,'import'},'Tag','importcustom');

set(m5,'Position',8,'Separator','on');
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_PrintToFigure')), 'Position',9,...
           'Callback','dfittool(''duplicate'')', 'Tag','dfitMenuPrint2Fig');
set(m7,'Position',10,'Separator','on');

% Fix VIEW menu
h0 = findall(h,'flat', 'Tag','figMenuView');
h1 = findall(h0, 'Type','uimenu', 'Parent',h0);
delete(h1);
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_Legend')), 'Position',1,'Separator','off',...
           'Callback','dfittool(''togglelegend'')', 'Checked','on',...
           'Tag','showlegend');
dfgetset('showlegend','on');
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_Grid')), 'Position',2,...
           'Callback','dfittool(''togglegrid'')', 'Checked','off', ...
           'Tag','showgrid');
dfgetset('showgrid','off');
h1 = uimenu(h0, 'Label',getString(message('stats:dfstrings:label_ConfidenceLevel')),'Position',3,'Separator','on');
uimenu(h1, 'Label','9&0%', 'Position',1, ...
           'Callback','dfittool(''setconflev'',.90)','Tag','conflev90');
uimenu(h1, 'Label','9&5%', 'Position',2, 'Checked','on',...
           'Callback','dfittool(''setconflev'',.95)','Tag','conflev95');
uimenu(h1, 'Label','9&9%', 'Position',3, ...
           'Callback','dfittool(''setconflev'',.99)','Tag','conflev99');
uimenu(h1, 'Label',getString(message('stats:dfstrings:label_Other')), 'Position',4, ...
           'Callback','dfittool(''setconflev'',[])','Tag','conflevOther');
dfgetset('conflev',0.95);
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_ClearPlot')), 'Position',4,...
           'Callback','dfittool(''clear plot'')', 'Tag','dfitMenuClearPlot');

% Fix TOOLS menu
h0 = findall(h,'flat', 'Tag','figMenuTools');
h1 = findall(h0, 'Type','uimenu', 'Parent',h0);
for j=length(h1):-1:1
   mtag = get(h1(j),'Tag');
   if ~any(strcmp(mtag,{'figMenuPan' 'figMenuZoomOut' 'figMenuZoomIn'}))
     delete(h1(j));
     h1(j) = [];
   else
      set(h1(j),'Separator','off');
   end
end
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_AxesLimitControl')), 'Position',4, 'Separator','on', ...
           'Callback','dfittool(''toggleaxlimctrl'')', 'Checked','off', ...
           'Tag','showaxlimctrl');
dfgetset('showaxlimctrl','off');
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_DefaultAxesLimits')), 'Position',5, ...
           'Callback','dfittool(''defaultaxes'')', 'Tag','dfitMenuDefaultAxes');
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_SetDefaultBinRules')), 'Position',6, 'Separator','on', ...
           'Callback', @setDefaultBinWidthRules, 'Tag','setbinrules');
           

% Fix HELP menu
h0 = findall(h,'flat', 'Tag','figMenuHelp');
h1 = findall(h0, 'Type','uimenu', 'Parent',h0);
delete(h1);
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_StatisticsToolboxHelp')), 'Position',1,'Callback',...
       'doc stats', 'Tag','dfitMenuHelpTbx');
uimenu(h0, 'Label', getString(message('stats:dfstrings:label_DistributionFittingToolHelp')), 'Position',2,'Callback',...
        @(varargin) dfhelpviewer('distribution_fitting', 'dfittool'), ...
        'Tag','dfitMenuHelpDfit');
uimenu(h0, 'Label',getString(message('stats:dfstrings:label_Demos')), 'Position',3,'Separator','on','Callback',...
       'demo toolbox stat', 'Tag','dfitMenuDemos'); 

% ------------------------------------

function setDefaultBinWidthRules(varargin)
% SETDEFAULTBINWITHRULES Callback for Set Default Bin Rules

binWidth = awtinvoke('com.mathworks.toolbox.stats.BinWidth', 'getBinWidth');
awtinvoke(binWidth, 'displayDefaultBinWidth');
