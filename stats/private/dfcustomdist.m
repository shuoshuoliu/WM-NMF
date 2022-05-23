function dfcustomdist(~,~,action)
%DFCUSTOM Callbacks for menu items related to custom distributions


%   Copyright 2003-2011 The MathWorks, Inc.

dft = com.mathworks.toolbox.stats.DistributionFitting.getDistributionFitting;

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
      dists = dfgetdistributions([],true,true,true); % force refresh
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
   dfsetdistributions(dft,dists);

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
    'name'                   getString(message('stats:dfstrings:cellstr_ImportedDistributions'))  ...
    'color'                  figcol ...
    'resize'                 'off' ...
    'numbertitle'            'off' ...
    'menubar'                'none' ...
    'windowstyle'            'modal' ...
    'visible'                'off' ...
    'integerhandle'          'off'    ...
    'handlevisibility'       'callback' ...
    'position'               fp   ...
    'closerequestfcn'        'delete(gcbf)' ...
    'Dock'                   'off' ...
            };
fig = figure(fig_props{:});

posn = [ffs+fus     fp(4)-(ffs+fus+ex) ...
        listsize(1) ex];

uicontrol('style','text','string',promptstring,...
          'horizontalalignment','left','position',posn);

btn_wid = (fp(3)-2*(ffs+fus)-fus)/2;
liststring=cellstr(liststring);
uicontrol('style','listbox',...
                    'position',[ffs+fus ffs+uh+4*fus+footnoteheight listsize],...
                    'string',liststring,...
                    'backgroundcolor',figcol,...
                    'max',2,...
                    'tag','listbox',...
                    'value',[]);

if footnote
   uicontrol('style','text',...
             'tag','footnote',...
             'string',getString(message('stats:dfstrings:sprintf_ImportedOrChanged')),...
             'horizontalalignment','left',...
             'position',[ffs+fus, ffs+fus+uh+footnoteheight/4, listsize(1), footnoteheight]);
end


uicontrol('style','pushbutton',...
                   'string','OK',...
                   'tag','OK',...
                   'position',[ffs+fus+listsize(1)/2-btn_wid/2 ffs+fus btn_wid uh],...
                   'callback','delete(gcbf)');

% make sure we are on screen
placetitlebar(fig)
set(fig, 'visible','on');
