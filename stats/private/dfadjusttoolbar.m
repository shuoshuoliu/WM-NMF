function dfadjusttoolbar(dffig)
%DFADJUSTTOOLBAR Adjust contents of distribution fitter plot toolbar


%   Copyright 2003-2012 The MathWorks, Inc.

h0 = findall(dffig,'Type','uitoolbar');
h1 = findall(h0,'Parent',h0);
czoom = [];
for j=length(h1):-1:1
   mtag = get(h1(j),'Tag');
   if any(strcmp(mtag,{'Exploration.Pan' 'Exploration.ZoomOut' 'Exploration.ZoomIn'}))
      czoom(end+1) = h1(j);
   elseif ~strcmp(mtag,'Standard.PrintFigure')
      delete(h1(j));
      h1(j) = [];
   else
     c1 = h1(j);
   end
end

% Add more icons especially for distribution fitter
state = dfgetset('showlegend');
if isempty(state), state = 'on'; end
c2 = uitoolfactory(h0,'Annotation.InsertLegend');
set(c2, 'State',state,...
        'TooltipString', getString(message('stats:dfstrings:tooltip_LegendOnOff')), ...
        'Separator','on',...
        'ClickedCallback','dfittool(''togglelegend'')',...
        'Tag','showlegend');
    
if exist('dficons.mat','file')==2
   icons = load('dficons.mat','icons');
   state = dfgetset('showgrid');
   if isempty(state), state = 'off'; end
   c3 = uitoggletool(h0, 'CData',icons.icons.grid,...
                    'State',state,...
                    'TooltipString', getString(message('stats:dfstrings:tooltip_GridOnOff')), ...
                    'Separator','off',...
                    'ClickedCallback','dfittool(''togglegrid'')',...
                    'Tag','showgrid');
   c4 = uipushtool(h0, 'CData',icons.icons.resetview,...
                    'TooltipString',getString(message('stats:dfstrings:tooltip_RestoreLimits')), ...
                    'Separator','off',...
                    'ClickedCallback','dfittool(''defaultaxes'')', ...
                    'Tag','defaultaxes');
   cnew = [c1 czoom c2 c3 c4]';
   
   set(h0,'Children',cnew(end:-1:1));
end
