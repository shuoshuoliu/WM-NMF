function customdist(~,~,action)
%CUSTOM Callbacks for menu items related to custom distributions


%   Copyright 2003-2020 The MathWorks, Inc.

dft = stats.internal.dfit.DistributionFitting.getInstance();

switch(action)

  % --------------------------------
  case 'define'       % define a file of probability distribution specs
      % Start editing a new file with default contents
      fname = fullfile(matlabroot,'toolbox','stats','stats','private',...
          'dftoolinittemplate.m');
      txt = fileread(fname);
      matlab.desktop.editor.newDocument(txt);
      
      % Display a helpful message about what's going on
      msg = getString(message('stats:dfstrings:cellstr_DefineCustomEdit'));
      msgbox(msg,getString(message('stats:dfstrings:dlg_DefineCustom')),'none','modal');
      dfittool('adjustlayout'); % in case figure is disturbed above

  % --------------------------------
  case 'import'       % refresh the distribution list 
   % Remember current distributions
   olds = dfgetset('alldistributions');
   
   try
      makedist -reset
      dists = stats.internal.dfit.getdistributions([],true,true,true); % force refresh
      newrows = find(~ismember(lower({dists.name}),lower({olds.name})));
      errmsg = '';
   catch ME
      errmsg = ME.message;
   end
   
   % Revert to previous distribution list if anything bad happened
   if ~isempty(errmsg)
      dists = olds;
      errordlg(getString(message('stats:dfstrings:dlg_ErrorImportCustom',...
                       errmsg)),...
               getString(message('stats:dfstrings:dlg_ImportCustomDistributions')),'modal');
      newrows = [];
   end

   % Sort by name
   lowernames = lower(char(dists.name));
   [~, ind] = sortrows(lowernames);
   dists = dists(ind);
   newrows = find(ismember(ind,newrows));

   if isempty(errmsg)
      showresults({dists.name},newrows);
   end
   stats.internal.dfit.setdistributions(dft,dists);

   dfgetset('dirty',true);   % session has changed since last save
end

% ---------------------------------
function showresults(liststring,asterisk)
%SHOWRESULTS Stripped-down version of listdlg, just to show a list

promptstring = getString(message('stats:dfstrings:assignment_NewDistributionList'));

if nargin>=2
   for j=1:length(asterisk)
      liststring{asterisk(j)} = sprintf('%s *',liststring{asterisk(j)});
   end
   footnote = ~isempty(asterisk);
else
   footnote = false;
end

ex = get(0,'defaultuicontrolfontsize')*1.7;  % extent height per line
fp = get(0,'defaultfigureposition');
fus = 8;       % frame/uicontrol spacing
ffs = 8;       % frame/figure spacing
uh = 22;       % uicontrol button height
listsize = [160 300];
if footnote
   footnoteheight = 2*ex;
else
   footnoteheight = 0;
end

w = 2*(fus+ffs)+listsize(1);
h = 2*ffs+6*fus+ex+listsize(2)+uh + footnoteheight;
fp = [fp(1) fp(2)+fp(4)-h w h];  % keep upper left corner fixed

figcol = get(0,'defaultUicontrolBackgroundColor');
fig_props = { ...
    'Name'                   getString(message('stats:dfstrings:cellstr_ImportedDistributions'))  ...
    'Color'                  figcol ...
    'NumberTitle'            'off' ...
    'MenuBar'                'none' ...
    'WindowStyle'            'modal' ...
    'Visible'                'off' ...
    'IntegerHandle'          'off'    ...
    'HandleVisibility'       'callback' ...
    'Position'               fp   ...
    'Tag'                    'dfImportCustomDistributionFigure' ...
            };
fig = uifigure(fig_props{:});

posn = [ffs+fus     fp(4)-(ffs+fus+ex) ...
        listsize(1) ex];

uilabel(fig, 'Text', promptstring,...
          'HorizontalAlignment','left','Position',posn);

btn_wid = (fp(3)-2*(ffs+fus)-fus)/2;
liststring=cellstr(liststring);
uilistbox(fig, 'Position',[ffs+fus ffs+uh+4*fus+footnoteheight listsize],...
                    'Items',liststring,...
                    'BackgroundColor',figcol,...
                    'Tag','listbox',...
                    'Value',{});
if footnote
   uilabel(fig, 'Tag','footnote',...
             'Text',getString(message('stats:dfstrings:sprintf_ImportedOrChanged')),...
             'HorizontalAlignment','left',...
             'Position',[ffs+fus, ffs+fus+uh+footnoteheight/4, listsize(1), footnoteheight]);
end


uibutton(fig, 'Text','OK',...
                   'Tag','OK',...
                   'Position',[ffs+fus+listsize(1)/2-btn_wid/2 ffs+fus btn_wid uh],...
                   'ButtonPushedFcn','delete(gcbf)');

% make sure we are on screen
placetitlebar(fig)
set(fig, 'Visible','on');

function placetitlebar(fig)
%PLACETITLEBAR ensures that a figure's titlebar is on screen.

oldRootUnits = get(0, 'Units');
oldFigUnits = get(fig, 'Units');

set(0, 'Units', 'pixels');
set(fig, 'Units', 'pixels');
   
screenpos = get(0, 'Screensize');
outerpos = get(fig, 'Outerposition');
if outerpos(2) + outerpos(4) > screenpos(4)
    outerpos(2) = screenpos(4) - outerpos(4);
    set(fig, 'Outerposition', outerpos);
end
%restore units
set(0, 'Units', oldRootUnits);
set(fig, 'Units', oldFigUnits);

