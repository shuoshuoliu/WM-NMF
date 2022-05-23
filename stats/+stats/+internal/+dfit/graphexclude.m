function graphexclude(dsname,xlo,xhi)
%GRAPHEXCLUDE  Create graph for selecting (x,y) pairs to exclude
%   GRAPHEXCLUDE(EXCLUDEPANEL,DSNAME,LOBND,UPBND) creates a graph
%   tied to the exclusion panel EXCLUDEPANEL, for dataset DSNAME, with
%   current lower and upper bounds LOBND and UPBND.  It provides a graphical
%   way to modify those bounds.

% Copyright 2001-2020 The MathWorks, Inc.

% Use old figure if any, or create a new one with a plot of the data
t = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');
c = get(0,'Child');
f = findobj(c,'flat','Type','figure','Tag','dfexcludegraph');
set(0,'ShowHiddenHandles',t);
if ~isempty(f)
   subfig = f;
else
   subfig = setupfigure(dsname,xlo,xhi);
end
if isempty(subfig)
   return
end
figure(subfig)

% Adjust the patches to show the desired
ax = findobj(subfig, 'Tag', 'exclude_graphically_axes');
xlim = get(ax,'XLim');

% If bounds already exist, put them onto the graph
if nargin>=2 && ~isempty(xlo) && ~isinf(xlo) && xlo > xlim(1)
   addremovepatch(subfig,'lower',xlo,'add')
end
if nargin>=3 && ~isempty(xhi) && ~isinf(xhi) && xhi < xlim(2)
   addremovepatch(subfig,'upper',xhi,'add')
end

set(ax,'XLimMode','manual','YLimMode','manual');
dfgetset('dfsubfig',subfig);

return

% -------------- Create figure
function subfig = makefigure(x,y,dolegend,xlo,xhi)
figcolor = get(0,'defaultuicontrolbackgroundcolor');
subfig = uifigure('IntegerHandle','off',...
    'name',getString(message('stats:dfstrings:fig_DefineBoundary')), ...
    'numbertitle','off',...
    'color',figcolor,...
    'Tag','dfexcludegraph');

% Create grids
mainGrid = uigridlayout(subfig, [2, 1]);
mainGrid.ColumnWidth = {'1x'};
mainGrid.RowHeight = {'1x', 'fit', 1, 'fit'};

middleGrid = uigridlayout(mainGrid);
middleGrid.ColumnWidth = {'fit', '1x'};
middleGrid.RowHeight = {'fit'};
middleGrid.Layout.Row = 2;
middleGrid.Layout.Column = 1;
middleGrid.Padding = [0, 0, 0, 0];

horizontalLine = uipanel(mainGrid);
horizontalLine.BackgroundColor = [.80 .80 .80];
horizontalLine.Layout.Row = 3;
horizontalLine.Layout.Column = 1;

bottomGrid = uigridlayout(mainGrid);
bottomGrid.ColumnWidth = {'fit', '1x', 'fit'};
bottomGrid.RowHeight = {'fit'};
bottomGrid.Layout.Row = 4;
bottomGrid.Layout.Column = 1;
bottomGrid.Padding = [0,0,0,0];

middleLeftGrid = uigridlayout(middleGrid);
middleLeftGrid.ColumnWidth = {'fit'};
middleLeftGrid.RowHeight = {'fit', 'fit'};
middleLeftGrid.Layout.Row = 1;
middleLeftGrid.Layout.Column = 1;
middleLeftGrid.ColumnSpacing = 0;
middleLeftGrid.RowSpacing = 5;
middleLeftGrid.Padding = [0, 0, 0, 0];

middleRightGrid = uigridlayout(middleGrid);
middleRightGrid.ColumnWidth = {'1x','fit','fit'};
middleRightGrid.RowHeight = {'fit'};
middleRightGrid.Layout.Row = 1;
middleRightGrid.Layout.Column = 2;
middleRightGrid.Padding = [0, 0, 0, 10];

% Add panel for Axes
p = uipanel('Parent',mainGrid,...
    'Title', '',...
    'BorderType', 'none');
p.Layout.Row = 1;
p.Layout.Column = 1;

% Add axes, for now in default position
xlim = [min(x) max(x)];
if isfinite(xlo)
    xlim(1) = min(xlim(1), xlo);
end
if isfinite(xhi)
    xlim(2) = max(xlim(2), xhi);
end
xlim = xlim + .1 * [-1 1] * diff(xlim);
ax = axes('Parent',p, ...
    'Box','on', ...
    'HitTest','on',...
    'Tag', 'exclude_graphically_axes', ...
    'XLim',xlim, ...
    'YLim',[0, max(y)+1]);
axtoolbar(ax, {'export','pan','zoomin', 'zoomout','restoreview'});
disableDefaultInteractivity(ax);

% Add control buttons
iCreateDummyMinSizeButton(middleRightGrid, 1, 2);
iCreateDummyMinSizeButton(middleRightGrid, 1, 3);

lowerLimitButton = uibutton(middleRightGrid,...
    'Text',getString(message('stats:dfstrings:button_AddLowerLimit')), ...
    'Tag','lower', ...
    'ButtonPushedFcn', @buttoncallback);
lowerLimitButton.Layout.Row = 1;
lowerLimitButton.Layout.Column = 2;

upperLimitButton = uibutton(middleRightGrid,...
    'Text',getString(message('stats:dfstrings:button_AddUpperLimit')), ...
              'Tag','upper', ...
              'ButtonPushedFcn',@buttoncallback);
upperLimitButton.Layout.Row = 1;
upperLimitButton.Layout.Column = 3;

% Add help and close buttons
iCreateDummyMinSizeButton(bottomGrid, 1, 1);
iCreateDummyMinSizeButton(bottomGrid, 1, 3);

helpButton = uibutton(bottomGrid, ...
    'Text', getString(message('stats:dfittool:button_help')), ...
    'Tag', 'dfGraphExcludeHelpButton', ...
    'ButtonPushedFcn', @helpButton);
helpButton.Layout.Row = 1;
helpButton.Layout.Column = 1;

closeButton = uibutton(bottomGrid, ...
    'Text',getString(message('stats:dfstrings:button_Close')), ...
    'Tag','dfGraphExcludeCloseButton',...
    'ButtonPushedFcn',@done);
closeButton.Layout.Row = 1;
closeButton.Layout.Column = 3;

setappdata(subfig,'buttons',[closeButton upperLimitButton lowerLimitButton]);

% Place text as well
censoredOpenCirclesText = uilabel(middleLeftGrid, ...
    'Text', getString(message('stats:dfstrings:uitext_ObservedFilledCircles')), ...
    'Tag', 'dfGraphExcludeCensoredOpenCirclesLabel');
censoredOpenCirclesText.Layout.Row = 1;
censoredOpenCirclesText.Layout.Column = 1;

observedFilledCirclesText = uilabel(middleLeftGrid, ...
    'Text', getString(message('stats:dfstrings:uitext_CensoredOpenCircles')), ...
    'Tag', 'dfGraphExcludeObservedFilledCirclesLabel');
observedFilledCirclesText.Layout.Row = 2;
observedFilledCirclesText.Layout.Column = 1;

if ~dolegend
   censoredOpenCirclesText.Visible = 'off';
   observedFilledCirclesText.Visible = 'off';
end

set(subfig,'WindowButtonMotionFcn',@fixcursor);

% ------------------- helper function to set up figure
function subfig = setupfigure(dsname,xlo,xhi)
%SETUPFIGURE Set up figure to do graphical exclusion
% We're excluding based on data in one dataset
dsdb = stats.internal.dfit.getdsdb;
a = down(dsdb);
ds = [];
while(~isempty(a))
   if isequal(dsname,a.name)
      ds = a;
      break;
   end
   a = right(a);
end
if isempty(ds)
   subfig = [];
   return
end
   
[ydata,cens,freq] = getincludeddata(ds,[]); % get data w/o NaNs
if isempty(cens)
   cens = zeros(size(ydata));
end
if isempty(freq)
   freq = ones(size(ydata));
end

% Sort y and carry along the rest
[ydata,i] = sort(ydata);
cens = cens(i);
freq = freq(i);

% Create x and y vectors to plot
n = sum(freq);
x = zeros(n,1);
y = zeros(n,1);
g = zeros(n,1);
x(1:freq(1)) = ydata(1);
y(1:freq(1)) = (1:freq(1))';
g(1:freq(1)) = cens(1);
i = freq(1)+1;
for k=2:length(ydata)
   for j=1:freq(k)
      x(i) = ydata(k);
      g(i) = cens(k);
      if (i>1) && (x(i)==x(i-1))
         y(i) = y(i-1) + 1;
      else
         y(i) = 1;
      end
      i = i+1;
   end
end

% Make a figure to receive graph
dolegend = any(g==0) & any(g==1);
subfig = makefigure(x,y,dolegend, xlo, xhi);
ax = findobj(subfig, 'Tag', 'exclude_graphically_axes');

% Place data points into graph
t = (g==0);
if any(t)
   line('XData',x(t),'YData',y(t),'HitTest','off','PickableParts','none',...
          'Color','b','Marker','.','LineStyle','none',...
          'MarkerSize',24,'Parent',ax,'Tag','observed');
end
t = (g==1);
if any(t)
   line('XData',x(t),'YData',y(t),'HitTest','off','PickableParts','none',...
          'Color','b','Marker','o','LineStyle','none',...
          'Parent',ax,'Tag','censored');
end


% --------------------------------
function buttoncallback(varargin)

% Get some handles and dimensions
button = gcbo;
fig = gcbf;
buttontag = get(button,'Tag');
addremovepatch(fig,buttontag);
updateGUI;

% --------------------------------
function addremovepatch(fig,whichbound,xbnd,addremove)

% Get some handles and dimensions
ax = findobj(fig, 'Tag', 'exclude_graphically_axes');
xlim = get(ax,'XLim');
ylim = get(ax,'YLim');
dx = 0.05 * diff(xlim);
if nargin>=4 && isequal(addremove,'add')
   forceadd = true;
else
   forceadd = false;
end

if nargin<3
   if isequal(whichbound,'lower')
      xbnd = xlim(1) + dx;
   else
      xbnd = xlim(2) - dx;
   end
end

% Carry out requested action
if isequal(whichbound,'lower')
   hPatch = findall(fig,'Tag','lowerpatch');
   if isempty(hPatch) || forceadd
      otherpatch = findall(fig,'Tag','upperpatch');
      if ~isempty(otherpatch)
         % Never put new limit beyond the other limit
         otherx = get(otherpatch,'XData');
         otherx = otherx(2);
         xbnd = min(xbnd, xlim(1) + .9*(otherx-xlim(1)));
      end
      x = [xlim(1), xbnd,    xbnd,    xlim(1), xlim(1)];
      y = [ylim(1), ylim(1), ylim(2), ylim(2), ylim(1)];
      if isempty(hPatch)
         patch(x,y,[.8 .8 .8],'Parent',ax,'Tag','lowerpatch',...
            'FaceAlpha',0.6,'ButtonDownFcn',@startselect);
      else
         set(hPatch,'XData',x,'YData',y);
      end
         
      newtxt = getString(message('stats:dfstrings:button_RemoveLowerLimit'));
   else
      delete(hPatch);
      newtxt = getString(message('stats:dfstrings:button_AddLowerLimit'));
   end
else
   hPatch = findall(fig,'Tag','upperpatch');
   if isempty(hPatch) || forceadd
      otherpatch = findall(fig,'Tag','lowerpatch');
      if ~isempty(otherpatch)
         % Never put new limit beyond the other limit
         otherx = get(otherpatch,'XData');
         otherx = otherx(2);
         xbnd = max(xbnd, xlim(2) - .9*(xlim(2)-otherx));
      end
      x = [xlim(2), xbnd,    xbnd,    xlim(2), xlim(2)];
      y = [ylim(1), ylim(1), ylim(2), ylim(2), ylim(1)];
      if isempty(hPatch)
         patch(x,y,[.8 .8 .8],'Parent',ax,'Tag','upperpatch',...
            'FaceAlpha',0.6,'ButtonDownFcn',@startselect);
      else
         set(hPatch,'XData',x,'YData',y);
      end
      newtxt = getString(message('stats:dfstrings:button_RemoveUpperLimit'));
   else
      delete(hPatch);
      newtxt = getString(message('stats:dfstrings:button_AddUpperLimit'));
   end
end

% Update button text
button = findobj(fig,'Tag',whichbound);
set(button,'Text',newtxt);

% ------------- function to initiate graphical selection
function startselect(varargin)

% Get figure and axis handles, define functions to do and end selection
subfig = findall(0, 'Tag','dfexcludegraph');
ax = findobj(subfig, 'Tag', 'exclude_graphically_axes');

if(~strcmp(ax.InteractionContainer.CurrentMode,'none')) % exit callback if we are in pan/zoom mode
    return
end
    
% Get current exclusion limits, use axis limits if none
lims = get(ax,'XLim');
hPatch = findall(subfig,'Tag','lowerpatch');
if ~isempty(hPatch)
   x = get(hPatch,'XData');
   lims(1) = max(x(:));
end
hPatch = findall(subfig,'Tag','upperpatch');
if ~isempty(hPatch)
   x = get(hPatch,'XData');
   lims(2) = min(x(:));
end

% Save information for other functions
hPatch = gcbo;
set(subfig,'WindowButtonMotionFcn',{@movepatch hPatch},...
           'WindowButtonUpFcn',@endmove);
setappdata(ax,'limits',lims);
setappdata(ax,'objmoving',hPatch);


% ------------- function to update GUI
function updateGUI()
subfig = gcbf;

hPatch = findall(subfig,'Tag','lowerpatch');
if isempty(hPatch)
   xl = '';
else
   xl = get(hPatch,'XData');
   xl = num2str(xl(2));
end
hPatch = findall(subfig,'Tag','upperpatch');
if isempty(hPatch)
   xh = '';
else
   xh = get(hPatch,'XData');
   xh = num2str(xh(2));
end

ex = stats.internal.dfit.Exclude.getInstance();
ex.updateBoundsFields(xl, xh);

% ------------- function to complete graphical selection
function endmove(varargin)

% Turn off window functions to end selection
subfig = gcbf;
set(subfig,'WindowButtonMotionFcn',@fixcursor, 'WindowButtonUpFcn',[]);
updateGUI;

% ------------- move patch boundary
function movepatch(~,~,hPatch)

varargin
subfig = gcbf;
ax = findobj(subfig, 'Tag', 'exclude_graphically_axes');

% Get exclusion limits and axis limits
lims = getappdata(ax,'limits');
xlim = get(ax,'XLim');
delta = .01 * abs(xlim(2) - xlim(1));

% Extend patch to the current point, but within limits
cp = get(ax,'CurrentPoint');
x = cp(1);
if isequal(get(hPatch,'Tag'),'lowerpatch')
   lobnd = xlim(1) + min(delta, .9*(lims(2)-xlim(1)));
   upbnd = min(xlim(2),lims(2)) - delta;
   x = max(lobnd, min(x,upbnd));
   lims(1) = x;
else
   lobnd = max(xlim(1),lims(1)) + delta;
   upbnd = xlim(2) - min(delta, .9*(xlim(2)-lims(1)));
   x = min(upbnd, max(lobnd,x));
   lims(2) = x;
end

% Update saved limits, and x data for this patch
setappdata(ax,'limits',lims);
xdata = get(hPatch,'XData');
xdata(2:3) = x;
set(hPatch,'XData',xdata);

% --------------- set cursor if we're on something that can move
function fixcursor(~,eventData)
ptr = get(gcbf,'Pointer');
onpatch = isequal(get(eventData.HitObject,'Type'),'patch');
if isequal(ptr,'arrow')
    if onpatch
        set(gcbf,'Pointer','left');
    end
else
    if ~onpatch
        set(gcbf,'Pointer','arrow');
    end
end

% --------------- close figure
function done(varargin)

delete(gcbf);

% --------------- help button callback
function helpButton(varargin)

stats.internal.dfit.helpviewer('exclude_data', 'dfittool')

% --------------- helper to minimum-sized uibutton
function iCreateDummyMinSizeButton(grid, row, column)
dummyText = repmat('X',1,8);
minButton = uibutton(grid, 'Text', dummyText);
minButton.Layout.Row = row;
minButton.Layout.Column = column;
minButton.Visible = 'off';
