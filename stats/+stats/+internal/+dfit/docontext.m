function docontext(varargin)
%DOCONTEXT Perform context menu actions for Distribution Fitter

% Copyright 2019 The MathWorks, Inc.

% Special action to create context menus
if isequal(varargin{1},'create')
   makecontextmenu(varargin{2});
   return
end

% Get information about what invoked this function
obj = gcbo;
dffig = gcbf;
action = get(obj,'Tag');
h = gco(dffig);
if isempty(h), return; end

% Set up variables that define some menu items
[~, styles, markers] = getmenuitems;
styles{end+1} = 'none';

changed = true;   % did a line property change?

switch action

 % This case is triggered when we display the menu
 case {'fitcontext' 'datacontext'}
   % Store a handle to the object that triggered this menu
   set(obj,'UserData',h);

   hObject = get(h,'UserData');
   c = findall(obj,'Type','uimenu');
   hBounds = findall(c,'flat','Tag','confbounds');
   ftype = dfgetset('ftype');
   
   % Enable or disable as appropriate
   hMarker = findall(c,'flat','Tag','marker');
   if isequal(action,'datacontext')
      hLineStyle = findall(c,'flat','Tag','linestyle');
      hLineWidth = findall(c,'flat','Tag','linewidth');
      hBinRules = findall(c,'flat','Tag','binrules');
      if isequal(ftype,'probplot')
         set(hMarker,'Enable','on');
         set(hLineStyle,'Enable','off');
         set(hLineWidth,'Enable','off');
         set(hBinRules,'Enable','off');
      elseif isequal(ftype,'pdf')
         set(hMarker,'Enable','off');
         set(hLineStyle,'Enable','on');
         set(hLineWidth,'Enable','on');
         set(hBinRules,'Enable','on');
      else
         set(hMarker,'Enable','off');
         set(hLineStyle,'Enable','on');
         set(hLineWidth,'Enable','on');
         set(hBinRules,'Enable','off');
      end
   else
      if ~hObject.iscontinuous && isequal(ftype,'pdf')
         set(hMarker,'Enable','on');
      else
         set(hMarker,'Enable','off');
      end
   end
   try
      hasconfbounds = hObject.distspec.hasconfbounds;
   catch
      hasconfbounds = false;
   end
   if isequal(ftype,'pdf') || isequal(ftype,'probplot') || ...
      (isequal(action,'fitcontext') && ~hasconfbounds) || ...
      (isequal(action,'datacontext') && isequal(ftype,'icdf'))
      set(hBounds,'Enable','off');
   else
      set(hBounds,'Enable','on');
   end

   set(c,'Checked','off');

   % Fix check mark for confidence bounds
   if hObject.showbounds
      set(hBounds,'Checked','on');
   end

   % Fix check marks on line width and line style cascading menus
   w = get(h,'LineWidth');
   u = findall(c,'flat','Tag',num2str(w));
   if ~isempty(u)
      set(u,'Checked','on');
   end
   w = get(h,'LineStyle');
   u = findall(c,'flat','Tag',w);
   if ~isempty(u)
      set(u,'Checked','on');
   end
   w = get(h,'Marker');
   u = findall(c,'flat','Tag',w);
   if ~isempty(u)
      set(u,'Checked','on');
   end
   return
   
 % Remaining cases are triggered by selecting menu items
 case 'confbounds'
   hObject = get(h,'UserData');
   hObject.showbounds = ~hObject.showbounds;
   htag = get(h,'Tag');
   if isequal(htag,'dfdata')
       dataObj = stats.internal.dfit.Data.getInstance();
       dataObj.updateManageDataSetsTable();
   else % 'distfit'
       fitManagerObj = stats.internal.dfit.FitsManager.getInstance();
       fitManagerObj.updateFitsTable();
   end
   dfgetset('dirty',true);   % session has changed since last save

 case 'color'
   oldcolor = get(h,'Color');
   newcolor = uisetcolor(oldcolor);
   if ~isequal(oldcolor,newcolor)
      set(h,'Color',newcolor);
   end
   dfgetset('dirty',true);   % session has changed since last save

 case styles
   set(h,'LineStyle',action);
   dfgetset('dirty',true);   % session has changed since last save

 case markers
   if isequal(action,'point')
      msize = 12;
   else
      msize = 6;
   end
   set(h,'Marker',action,'MarkerSize',msize);
   dfgetset('dirty',true);   % session has changed since last save

% Hide a fit 
case 'hidefit'
    hndl = get(h,'UserData');
    hndl.plot = 0;
    fitManagerObj = stats.internal.dfit.FitsManager.getInstance();
    fitManagerObj.updateFitsTable();
    changed = false;
      
 % Delete a fit
 case 'deletefit'
    hndl = get(h,'UserData');
    fitManagerObj = stats.internal.dfit.FitsManager.getInstance();
    fitManagerObj.deleteFits(hndl);
    changed = false;

 % Edit a fit
 case 'editfit'
    htag = get(h,'Tag');
    if isequal(htag,'distfit')  % should always be true
        hndl = get(h,'UserData');
        fitEditor = hndl.fitframe;
        fitEditor.showPanel(); % bring FitEditor forward and make it visible
   end
   changed = false;

 % Hide data 
 case 'hidedata'
    hndl = get(h,'UserData');
    hndl.plot = 0;
    dataObj = stats.internal.dfit.Data.getInstance();
    dataObj.updateManageDataSetsTable();
    changed = false;

 % Bring up the "Set Bin Width Rules" dialog for this data set
 case 'binrules'
   htag = get(h,'Tag');
   if isequal(htag,'dfdata')
      hndl = get(h,'UserData');
      nm = get(hndl,'name');
      dataObj = stats.internal.dfit.Data.getInstance();
      dataObj.openBinWidthRulesPanel(nm);
   end
   changed = false;
 
 % If the menu item is a number, it is a line width
 otherwise
   j = str2num(action); %#ok<ST2NM>
   if ~isempty(j)
      set(h,'LineWidth',j);
   end
   dfgetset('dirty',true);   % session has changed since last save

end

if changed
   % Save plot info in the fit or data set object
   hObject = get(h,'UserData');
   savelineproperties(hObject);

   % Update legend
   stats.internal.dfit.updatelegend(dffig);
end


% ---------------------- helper to populate the shared context menu items
function makesharedcontextmenuitems(hMenu)
%MAKESHAREDCONTEXTMENUITEMS Populates a menu with the shared contextmenu items
uimenu(hMenu,'Label',getString(message('stats:dfstrings:label_Color')),'Tag','color','Callback',@stats.internal.dfit.docontext);

% Add menu items for line and marker control
umark = uimenu(hMenu,'Label',getString(message('stats:dfstrings:label_Marker')),'Tag','marker');
uwidth = uimenu(hMenu,'Label',getString(message('stats:dfstrings:label_LineWidth')),'Tag','linewidth');
ustyle = uimenu(hMenu,'Label',getString(message('stats:dfstrings:label_LineStyle')),'Tag','linestyle');

% Add menu items to control confidence bounds
uimenu(hMenu,'Label',getString(message('stats:dfstrings:label_ConfidenceBounds')),'Callback',@stats.internal.dfit.docontext,...
            'Tag','confbounds');

% Get menu item labels and tags
[sizes, styles, markers, slabels, mlabels] = getmenuitems;

for j=1:length(markers)
   uimenu(umark,'Label',mlabels{j},'Callback',@stats.internal.dfit.docontext,'Tag',markers{j});
end

% Sub-menus for line widths
for i = 1:length(sizes)
   val = num2str(sizes(i));
   uimenu(uwidth,'Label',val,'Callback',@stats.internal.dfit.docontext,'Tag',val);
end

% Sub-menus for line styles
for j=1:length(styles)
   uimenu(ustyle,'Label',slabels{j},'Callback',@stats.internal.dfit.docontext,'Tag',styles{j});
end

% ---------------------- helper to make context menu
function makecontextmenu(dffig)
%MAKECONTEXTMENU Creates context menu for Distribution Fitter figure

% Create context menus for fits, data curve, probability plot data curves
cFit = uicontextmenu('Parent',dffig,'Tag','fitcontext','Callback',@stats.internal.dfit.docontext);
makesharedcontextmenuitems(cFit);

% Create context menus for the data curve
cData = uicontextmenu('Parent',dffig,'Tag','datacontext','Callback',@stats.internal.dfit.docontext);
makesharedcontextmenuitems(cData);

% Add items for fit menus only
uimenu(cFit,'Label',getString(message('stats:dfstrings:label_HideFit')),'Tag','hidefit','Callback',@stats.internal.dfit.docontext,...
       'Separator','on');
uimenu(cFit,'Label',getString(message('stats:dfstrings:label_DeleteFit')),'Tag','deletefit','Callback',@stats.internal.dfit.docontext);
uimenu(cFit,'Label',getString(message('stats:dfstrings:label_EditFit')),'Tag','editfit','Callback',@stats.internal.dfit.docontext);

% Add items for data menus only
uimenu(cData,'Label',getString(message('stats:dfstrings:label_HideData')),'Tag','hidedata',...
       'Callback',@stats.internal.dfit.docontext,'Separator','on');
uimenu(cData,'Label',getString(message('stats:dfstrings:label_SetBinRules')),'Tag','binrules',...
       'Callback',@stats.internal.dfit.docontext,'Separator','on');

% -------------- helper to get menu item labels
function [sizes,styles,markers,slabels,mlabels] = getmenuitems
%GETMENUITEMS Get items for Distribution Fitter context menus
sizes = [0.5 1 2 3 4 5 6 7 8 9 10];
styles = {'-' '--' ':' '-.'};
markers = {'+' 'o' '*' '.' 'x' 'square' 'diamond' ...
        'v' '^' '<' '>' 'pentagram' 'hexagram'};
slabels = {'solid' 'dash' 'dot' 'dash-dot'};
mlabels = {'plus' 'circle' 'star' 'point' 'x-mark' 'square' 'diamond' ...
           'triangle (down)' 'triangle (up)' 'triangle (left)' ...
           'triangle (right)' 'pentagram' 'hexagram'};
