function outF=stepwise(X,y,in,penter,premove)
%STEPWISE Interactive tool for stepwise regression.
%   STEPWISE(X,Y) displays an interactive tool for creating a
%   regression model to predict the vector Y using a subset of the
%   predictors given by columns of the matrix X.  Initially no
%   predictors are included in the model, but you can click on
%   predictors to switch them into and out of the model.  STEPWISE
%   automatically includes a constant term in all models.
%
%   For each predictor in the model, its least squares coefficient
%   is plotted with a blue filled circle.  For each predictor not
%   in the model, a filled red circle indicates the coefficient it
%   would have if it were added to the model.  Horizontal bars indicate
%   90% (colored) and 95% (black) confidence intervals.
%
%   STEPWISE(X,Y,INMODEL,PENTER,PREMOVE) specifies the initial state
%   of the model and the confidence levels to use.  INMODEL is a logical
%   or index vector specifying the predictors that should be in the
%   initial model (default is none).  PENTER specifies the maximum
%   p-value for a predictor to be recommended for adding to the model
%   (default 0.05).  PREMOVE specifies the minimum p-value for a
%   predictor to be recommended for removal (default 0.10).
%
%   A NaN in either X or Y is treated as a missing value.  Rows
%   containing missing values are not used in the fit.
%
%   See also STEPWISEFIT.

%   Reference: Draper, Norman and Smith, Harry, Applied Regression Analysis,
%   Second Edition, John Wiley & Sons, Inc. 1981 pp. 307-312.

%   Copyright 1993-2020 The MathWorks, Inc.


if nargin==0
    % No data given, so launch stepwise using sample data.
    noargmsg = getString(message('stats:stepwise:dialogText_NoArgs'));
    noargtitle = getString(message('stats:stepwise:dialogTitle_NoArgs'));
    uiwait(msgbox(noargmsg,noargtitle,'modal'));

    s = load('hald');
    X = s.ingredients;
    y = s.heat;
    nargs = 2;
else
    nargs = nargin;
end

if nargs<2, error(message('stats:stepwise:InvalidSyntax')); end
if nargs<3 || isempty(in), in = false(1,size(X,2)); end
if nargs<4, penter = 0.05; end
if nargs<5, premove = max(penter,0.1); end

if penter>premove
   error(message('stats:stepwise:BadPValues'));
end
P = size(X,2);

if numel(y)~=length(y)
   error(message('stats:stepwise:InvalidData'));
end
if islogical(in)
   if length(in)~=P
      error(message('stats:stepwise:InvalidInmodelVector'));
   end
else
   if any(~ismember(in,1:P))
      error(message('stats:stepwise:InvalidInmodelVector'));
   end
   in = ismember((1:P),in);
end

% Remove NaN rows
[badin,wasnan,X,y] = statremovenan(X,y);
if (badin>0)
   error(message('stats:stepwise:LengthMismatch'))
end

% Create figure
F = makefig(P);

% Save info for callbacks
setappdata(F,'y',y);
setappdata(F,'allx',X);
setappdata(F,'in',in);
setappdata(F,'penter',penter);
setappdata(F,'premove',premove);
setappdata(F,'wasnan',wasnan);

% Perform regression and record it in the history
[in,stats] = stepmain(X,y,in,F);
addhistorystep(F,in,stats.rmse,stats.df0);
updateadvice(F,in);

if nargout>0
   outF = F;
end

% -----------------------------------------
function [in,stats] = stepmain(X,y,in,F)
% Perform steps until converged

P = length(in);

% Perform current fit
penter = getappdata(F,'penter');
premove = getappdata(F,'premove');
wasnan = getappdata(F,'wasnan');

[~,~,~,in,stats,nextstep] = ...
    stepwisefit(X,y,'inmodel',in,'penter',penter,'premove',premove,...
                'display','off','maxiter',0);

% NaN rows were already removed, so update the stats structure
stats.wasnan = wasnan;

% Update all parts of figure
setappdata(F,'in',in);
setappdata(F,'stats',stats);
setappdata(F,'advisedstep',nextstep);
updateplot(F,stats,in);
updatefooter(F,stats);
resizefig(F,[],P);

% -----------------------------------------
function addhistorystep(F,in,rmse,df0)
% Add computed step to the history record and plot

% Get current history state
H = getappdata(F,'history');
nhistorystats = 2;
if isempty(H)
   P = size(in,2);
   H = zeros(P+nhistorystats,0);
end

% Add one or more new states
nsteps = size(in,1);
col = size(H,2)+(1:nsteps);
H(1,col) = rmse;                    % root mean squared error
H(2,col) = df0;                     % number of terms "in" this model
H(nhistorystats+1:end,col) = in';   % boolean indicating "in" terms
setappdata(F,'history',H);
updatehistoryplot(F, col, rmse, in);

% -----------------------------------------
function stepcb(~,~,swap)

% See if this figure is ready to process a new action
F = gcbf;
ptr = get(F,'Pointer');
set(F,'Pointer','watch');

axes(getappdata(F,'plotaxes'));
dostep(F,swap);

set(F,'Pointer',ptr);


% -----------------------------------------
function dostep(F,swap)

% Swap this X variable, re-fit, and record in history
[X,y,in] = getdata(F);
in(swap) = ~in(swap);
[in,stats] = stepmain(X,y,in,F);
addhistorystep(F,in,stats.rmse,stats.df0);
updateadvice(F,in);


% -----------------------------------------
function F = makefig(P)

% Create a figure with a real set of axes and a fake set into which
% the coefficients table will be drawn

gray = get(0,'defaultuicontrolbackgroundcolor');
F = figure('Color',gray,'Visible','off','HandleVisibility','callback',...
           'IntegerHandle','off',...
           'Name',getString(message('stats:stepwise:dialogTitle_StepwiseRegression')),...
           'Units','points', 'NumberTitle','off','PaperPositionMode','auto',...
           'WindowStyle', 'Normal', ...
           'DockControls', 'off', ...
           'MenuBar','figure', ...
           'ToolBar','none', ...
           'Tag','stepwisefigure');
p = get(F,'Position');
p(3) = 550;
p(4) = 340;
set(F,'Position',p);

hMainAxes = axes('Units','points','Box','on',...
         'Ylim',[.5 P+.5],'YDir','rev','Parent',F,'Tag','mainAxes');

mainAxesTitle = getString(message('stats:stepwise:mainAxesTitle_Default'));
title(hMainAxes,mainAxesTitle,'FontSize',12);

hTableAxes = axes('Units','points','Box','on',...
          'YTick',[],'XTick',[],'Color',gray,'Parent',F,'Tag','tableAxes');
set(hTableAxes,'Ydir','rev','Ylim',[.5 P+.5],'XLim',0:1);
tableAxesTitle = getString(message('stats:stepwise:tableAxesTitle'));
title(hTableAxes,tableAxesTitle, ...
    'FontName', get(0,'FixedWidthFontName'), 'FontSize',12);

hHistoryAxes = axes('Units','points','Box','on',...
         'Xlim',[.5 1.5],'Visible','on','Parent',F,'tag','ha');
historyAxesTitle = getString(message('stats:stepwise:historyAxesTitle'));
title(hHistoryAxes,historyAxesTitle,'FontSize',12);
set(F,'WindowButtonMotionFcn',@histtips);

% Store axes handles into figure right away
setappdata(F,'plotaxes',hMainAxes);
setappdata(F,'tableaxes',hTableAxes);
setappdata(F,'historyaxes',hHistoryAxes);

% Figure out this font's text height by creating a temporary text object
temptext = getString(message('stats:stepwise:footerText_AdjRSQ',...
                                    sprintf('%g',0)));
h = text(0,0,temptext,'Parent',hMainAxes);
set(h,'Units','points');
e = get(h,'Extent');
delete(h);
textheight = 1.25*e(4);
setappdata(F,'textheight',textheight);

% Add a slider to the right of the top two axes
fpos = get(F,'Position');
apos = get(hMainAxes,'Position');

sliderwidth = 15;
spos = [fpos(3)-sliderwidth apos(2) sliderwidth apos(4)];
uicontrol(F,'Style','slider','Units','points','Position',spos,...
            'Value',1,'Callback',{@scrollfig P},'Tag','slider');

% Create a footer under the top two axes using a frame and more text objects
clr = get(F,'Color');
colwidth = (apos(3)-20)/3;
uipanel(F,'Units','points',...
          'BackgroundColor',clr,'Tag','frame',...
          'Position',[apos(1) 0 apos(3) 2.5*textheight]);

uicontrol(F,'Style','text','Units','points','Tag','intercept',...
          'BackgroundColor',clr,...
          'Position',[apos(1)+5 1.2*textheight colwidth textheight],...
          'HorizontalAlignment','left');
uicontrol(F,'Style','text','Units','points','Tag','rmse',...
          'BackgroundColor',clr,...
          'Position',[apos(1)+5 0.1*textheight colwidth textheight],...
          'HorizontalAlignment','left');

uicontrol(F,'Style','text','Units','points','Tag','rsq',...
          'BackgroundColor',clr,...
          'Position',[apos(1)+colwidth+5 1.2*textheight colwidth textheight],...
          'HorizontalAlignment','left');
uicontrol(F,'Style','text','Units','points','Tag','adjrsq',...
          'BackgroundColor',clr,...
          'Position',[apos(1)+colwidth+5 0.1*textheight colwidth textheight],...
          'HorizontalAlignment','left');

p = [apos(1)+2*colwidth+5 1.2*textheight colwidth textheight];
uicontrol(F,'Style','text','Units','points','Tag','fstat',...
          'BackgroundColor',clr,...
          'Position',p, 'HorizontalAlignment','left');
p = [apos(1)+2*colwidth+5 0.1*textheight colwidth textheight];
uicontrol(F,'Style','text','Units','points','Tag','pval',...
          'BackgroundColor',clr,...
          'Position',p, 'HorizontalAlignment','left');

% Create the buttons
buttonPos = [0 0 95 19.125];
hSaveButton = uicontrol(F,...
                        'String',getString(message('stats:stepwise:button_Export')),...
                        'Units','points',...
                        'Tag','savebutton','Position',buttonPos,...
                        'TooltipString',getString(message('stats:stepwise:tooltip_SaveResults')),...
                        'Callback',{@savebutton, F});
hStepButton = uicontrol(F,...
                        'String',getString(message('stats:stepwise:button_NextStep')),...
                        'Units','points',...
                        'Tag','stepbutton','Position',buttonPos,...
                        'TooltipString',getString(message('stats:stepwise:tooltip_NextStep')),...
                        'Callback',{@stepbutton 'one'});
hAutoButton = uicontrol(F,...
                        'String',getString(message('stats:stepwise:button_AllSteps')),...
                        'Units','points',...
                        'Tag','autobutton','Position',buttonPos,...
                        'TooltipString',getString(message('stats:stepwise:tooltip_AllSteps')),...
                        'Callback',{@stepbutton 'all'});
hAdviceLabel = uicontrol(F, 'Style', 'text',...
                        'String',getString(message('stats:stepwise:adviceLabel')),...
                        'Units','points','Position',[10 10 100 20],...
                        'HorizontalAlignment','left',...
                        'FontName', 'Helvetica','FontSize', 10);
extent = get(hAdviceLabel, 'Extent');
pos = get(hAdviceLabel, 'Position');
pos(3) = extent(3);
set(hAdviceLabel, 'Position', pos);

advicetext = getString(message('stats:stepwise:adviceText_NoMove'));
hAdviceText = uicontrol(F,'Style','text',...
                        'String',advicetext,...
                        'Units','points','Position',[10 10 100 20],...
                        'HorizontalAlignment','left', 'BackgroundColor', [1 1 .85],...
                        'FontName', 'Helvetica', 'FontSize', 10);
setappdata(F,'savebutton',hSaveButton);
setappdata(F,'stepbutton',hStepButton);
setappdata(F,'autobutton',hAutoButton);
setappdata(F,'advicetext',hAdviceText);
setappdata(F,'advicelabel',hAdviceLabel);

% Update menus to show only stepwise-related items
adjustmenu(F);

% Show figure now that it is created
set(F,'Visible','on');

% Set resize call back function
set(F,'ResizeFcn',{@resizefig P})

% -----------------------------------------
function scrollfig(hslider,~,P,varnum)
% inputs are slider handle, event, number of vars, desired var

F = get(hslider,'Parent');
hMainAxes = getappdata(F,'plotaxes');
hTableAxes = getappdata(F,'tableaxes');

% Get slider value, but force to an integer
v = get(hslider,'Value');
if v~=round(v)
   v = round(v);
end

% Adjust axis y limits to scroll lines on plot
smax = get(hslider,'Max');
nvars = max(1,P-smax+1);
ylim = [P-v-nvars+2-.5 P-v+1+.5];

% Make sure requested variable is visible
if nargin>=4 && ~isempty(varnum) && varnum>0
    if varnum<ylim(1) || varnum>ylim(2)
        ylo = max(1,varnum-floor(nvars/2));
        yhi = ylo+nvars-1;
        ylim = [ylo-0.5, yhi+0.5];
        v = P-ylo-nvars+2;
    end
end
set(hMainAxes,'YLim',ylim);
set(hslider,'Value',v);

% Hide parts of table not within view
h = getappdata(F,'tablehandles');
hy = 1:length(h);
set(h(hy<ylim(1) | hy>ylim(2)),'Visible','off');
set(h(hy>=ylim(1) & hy<=ylim(2)),'Visible','on');
in = getappdata(F,'in');
set(h(in),'Color',[0 0 1]);
set(h(~in),'Color','r');

set(hTableAxes,'YLim',ylim);


% -----------------------------------------
function resizefig(F,~,P)
hMainAxes = getappdata(F,'plotaxes');
fpos = get(F,'Position');
apos = get(hMainAxes,'Position');
textheight = getappdata(F,'textheight');

hSaveButton = getappdata(F,'savebutton');
hStepButton = getappdata(F,'stepbutton');
hAutoButton = getappdata(F,'autobutton');
hAdviceText = getappdata(F,'advicetext');
hAdviceLabel = getappdata(F,'advicelabel');
p = get(hSaveButton, 'Position');
buttonmargin = 15;
buttonwidth = p(3) + 2*buttonmargin;   % width of whole button area

% Get slider and its width; needed to lay out other things in figure
hslider = findobj(F,'Tag','slider');
spos = get(hslider,'Position');
sliderwidth = spos(3);
tablewidth = getappdata(F,'tablewidth');

% Lay out axes
footerheight = 2.5*textheight;
xlabelheight = 1.5*textheight;
titleheight = 2*textheight;
totalaxesheight = max(1,fpos(4)-footerheight-2*xlabelheight-2*titleheight);
topaxesheight = (2/3) * totalaxesheight;
bottomaxesheight = (1/3) * totalaxesheight;
axesbase = bottomaxesheight + footerheight + 2*xlabelheight + titleheight;
apos(2) = axesbase;
apos(3) = max(1,fpos(3) - apos(1) - sliderwidth - tablewidth - buttonwidth - 1);
apos(4) = max(1,topaxesheight);
set(hMainAxes,'Position',apos);

% Lay out table axes
hTableAxes = getappdata(F,'tableaxes');
set(hMainAxes,'Units','points');
a2pos = apos;
a2pos(1) = apos(1)+apos(3);
a2pos(3) = tablewidth;
set(hTableAxes,'Position',a2pos);
tableedge = a2pos(1) + a2pos(3);

% Adjust slider and axis limits if necessary
nvars = min(P, floor(apos(4)/(1.2*textheight)));
smin = 1;
smax = max(1,P-nvars+1);
if smin==smax
   sstep = [1 1];
else
   sstep = ([1 min(nvars,smax-smin)] ./ (smax-smin));
end
sval = max(1, min(P-nvars+1, get(hslider,'Value')));
if (smax == smin)
   set(hslider,'Min',smin,'Max',smax+eps,'SliderStep',sstep,'Value',sval,...
               'Visible','off');
else
   set(hslider,'Min',smin,'Max',smax,'SliderStep',sstep,'Value',sval,...
               'Visible','on');
end

scrollfig(hslider,[],P);
if ~isempty(hslider) && ishghandle(hslider)
   set(hslider,'Position',[tableedge apos(2) sliderwidth apos(4)]);
end

% Lay out button and text
labelp = get(hAdviceLabel,'Position');
labele = get(hAdviceLabel, 'Extent');
labelp(4) = labele(4);
labelp(1) = tableedge + sliderwidth + buttonmargin;
labelp(2) = apos(2) + apos(4) - labelp(4);
set(hAdviceLabel,'Position',labelp);
textp = get(hAdviceText,'Position');
textp(4) = labele(4);
textp(1) = labelp(1);
textp(2) = labelp(2) - textp(4);
set(hAdviceText,'Position',textp);
p = get(hStepButton,'Position');
buttonspacer = 5;
p(1) = textp(1);
p(2) = textp(2) - p(4) - buttonspacer;
set(hStepButton,'Position',p);
p(2) = p(2) - p(4) - (buttonspacer*2);
set(hAutoButton,'Position',p);
p(2) = min(apos(2), p(2) - p(4) - buttonspacer);
set(hSaveButton,'Position',p);

% Lay out footer
ftrwidth = fpos(3) - apos(1) - 20;
colwidth = (ftrwidth-20)/3;

h = findobj(F,'Tag','frame');
footerbase = bottomaxesheight + titleheight + xlabelheight;
framepos = [apos(1) footerbase ftrwidth 2.5*textheight];

if ~isempty(h) && ishghandle(h)
   set(h,'Position',framepos);
end

h = findobj(F,'Tag','intercept');
if ~isempty(h) && ishghandle(h)
    set(h,'Position',[apos(1)+5 footerbase+1.2*textheight colwidth textheight]);
end

h = findobj(F,'Tag','rmse');
if ~isempty(h) && ishghandle(h)
    set(h,'Position',[apos(1)+5 footerbase+0.1*textheight colwidth textheight]);
end

h = findobj(F,'Tag','rsq');
if ~isempty(h) && ishghandle(h)
   set(h,'Position',[apos(1)+colwidth+5 footerbase+1.2*textheight colwidth textheight]);
end

h = findobj(F,'Tag','adjrsq');
if ~isempty(h) && ishghandle(h)
   set(h,'Position',[apos(1)+colwidth+5 footerbase+0.1*textheight colwidth textheight]);
end

h = findobj(F,'Tag','fstat');
if ~isempty(h) && ishghandle(h)
   set(h,'Position',...
         [apos(1)+2*colwidth+5 footerbase+1.2*textheight colwidth textheight]);
end

h = findobj(F,'Tag','pval');
if ~isempty(h) && ishghandle(h)
   set(h,'Position',...
         [apos(1)+2*colwidth+5 footerbase+0.1*textheight colwidth textheight]);
end

% Lay out history axes
hHistoryAxes = getappdata(F,'historyaxes');
pos = get(hHistoryAxes,'Position');
pos(1) = apos(1);
pos(2) = xlabelheight;
pos(3) = framepos(1) + framepos(3) - apos(1);
pos(4) = max(1,bottomaxesheight);
set(hHistoryAxes,'Position',pos);
hLabel = get(hHistoryAxes,'YLabel');
set(hLabel,'String',getString(message('stats:stepwise:historyAxesYlabel')));
e = get(hLabel,'Extent');
xlim = get(hHistoryAxes,'XLim');
ratio = pos(3) ./ (xlim(2)-xlim(1));
leftside = pos(1) + (e(1)-xlim(1))*ratio;
if leftside<=0
   rightside = pos(1)+pos(3);
   pos(1) = 1 - (e(1)-xlim(1))*ratio;
   pos(3) = rightside - pos(1);
   set(hHistoryAxes,'Position',pos);
end


% -----------------------------------------
function updateplot(F,stats,in,penter,premove)

if nargin<4
    penter = getappdata(F,'penter');
end
if nargin<5
    premove = getappdata(F,'premove');
end
B = stats.B;
SE = stats.SE;

hMainAxes = getappdata(F,'plotaxes');
hTableAxes = getappdata(F,'tableaxes');
cla(hMainAxes);

% Create zero line to check for significance
P = length(in);
line(zeros(2,1),[-P;2*P],'Color','k','LineStyle',':','Marker','none',...
     'YLimInclude','off','Parent',hMainAxes);

% Collect names of coefficients and text sizes required
lbl = cell(P,1);
dfnow = stats.dfe;
ttl = text(0,0,'','Parent',hMainAxes);
set(ttl,'Visible','off','Units','points');
maxe = 0;
DF = dfnow - (in==0)';
minx = Inf;
maxx = -Inf;
for j=1:P
   % Generate coefficient name
   lbl{j} = getString(message('stats:stepwise:mainAxes_YTickLabel',sprintf('%d',j)));

   % Compute confidence interval
   if in(j)
      clr = [0 0 1];
      df = dfnow;
   else
      clr = [1 0 0];
      df = dfnow-1;
   end
   t1 = -tinv(premove/2,df);
   ci = B(j) + [-1;1] * t1 * SE(j);

   % Draw point and line
   cb = {@stepcb j};
   line(B(j),j,'Marker','.','LineStyle','none','Color',clr,...
           'MarkerSize',20,'ButtonDownFcn',cb,'Parent',hMainAxes);
   line(ci,[j;j],'Marker','none','LineStyle','-','Color',clr,...
           'LineWidth',2,'ButtonDownFcn',cb,'Parent',hMainAxes);
   minx = min(minx, min(ci));
   maxx = max(maxx, max(ci));
   if premove>penter
      t2 = -tinv(penter/2,df);
      ci2 = B(j) + [-1;1] * t2 * SE(j);
      line([ci2(1);ci(1);NaN;ci(2);ci2(2)],[j;j;j;j;j],'Marker','none',...
           'LineStyle','-','Color',[0 0 0],'LineWidth',1,...
           'ButtonDownFcn',cb,'Parent',hMainAxes);
      minx = min(minx, min(ci2));
      maxx = max(maxx, max(ci2));
   end

   % Measure space required for name
   set(ttl,'String',lbl{j});
   e = get(ttl,'Extent');
   maxe = max(maxe, e(3));
end
delete(ttl);

% Save numeric data for display
tstat = Inf(size(SE));
tstat(SE>0) = B(SE>0)./SE(SE>0);
pval = tcdf(tstat,DF);
pval = 2 * min(pval,1-pval);
tbl = [B tstat pval];

% Update table
h = getappdata(F,'tablehandles');
delete(h);
h = zeros(P,1);
for j=1:P
    cb = {@stepcb j};
    h(j) = text(0.02,j,sprintf('%12g %8.4f %8.4f',tbl(j,:)),...
        'Parent',hTableAxes,'Tag','coefftable','ButtonDownFcn',cb,...
        'HorizontalAlignment','left', ...
        'FontName',get(0,'FixedWidthFontName'));
end

% Measure the largest horizontal extent in points
e = get(h,'Extent');
e = reshape([e{:}],4,length(h));
[~,jmax] = max(e(3,:));

oldpos = get(h(jmax),'Position');
oldu = get(h(jmax),'Units');
set(h(jmax),'Units','points');
e = get(h(jmax),'Extent');
set(h(jmax),'Units',oldu,'Position',oldpos);
tablewidth = e(3) / 0.95;
setappdata(F,'tablewidth',tablewidth);
setappdata(F,'tablehandles',h);

% Set y tick labels equal to names
set(hMainAxes,'YTick',1:P,'YTickLabel',lbl,'Box','on');

% Set x axis limits now so they won't change as we scroll
if isinf(maxx) || isinf(minx)
    minx = min(B);
    maxx = max(B);
end
dx = (maxx - minx)/10;
set(hMainAxes, 'XLim',[minx-dx, maxx+dx]);

% Remember space required
apos = get(hMainAxes,'Position');
apos(1) = maxe + 10;
set(hMainAxes,'Position',apos);


% -----------------------------------------
function updateadvice(F,in)

buttons = [getappdata(F,'stepbutton') getappdata(F,'autobutton')];

% Figure out which term to move, enable buttons if any
termnum = getappdata(F,'advisedstep');
if termnum==0
   txt = getString(message('stats:stepwise:adviceText_NoMove'));
   set(buttons,'Enable','off');
else
   if in(termnum)==1
      txt = getString(message('stats:stepwise:adviceText_MoveOut',termnum));
   else
      txt = getString(message('stats:stepwise:adviceText_MoveIn',termnum));
   end
   set(buttons,'Enable','on');
end

% Update advice text
hAdviceText = getappdata(F,'advicetext');
set(hAdviceText,'String', txt);
extent = get(hAdviceText, 'Extent');
pos = get(hAdviceText, 'Position');
pos(3) = extent(3);
set(hAdviceText, 'Position', pos);

% Update highlighting -- yellow background for recommended term
hText = getappdata(F,'tablehandles');
set(hText,'BackgroundColor',get(F,'Color'));
if termnum>0
   set(hText(termnum),'BackgroundColor',[1 1 .85]);
   
   % Make sure this term is visible
   if isequal(get(hText(termnum),'Visible'),'off')
       hslider = findobj(F,'Tag','slider');
       P = length(in);
       scrollfig(hslider,[],P,termnum);
   end
end


% -----------------------------------------
function updatefooter(F,stats)
h = findobj(F,'tag','intercept');
if ~isempty(h) && ishghandle(h)
   set(h,'String',getString(message('stats:stepwise:footerText_Intercept',...
                                    sprintf('%g',stats.intercept))));
end
h = findobj(F,'tag','rmse');
if ~isempty(h) && ishghandle(h)
   set(h,'String',getString(message('stats:stepwise:footerText_RMSE',...
                                    sprintf('%g',stats.rmse))));
end
h = findobj(F,'tag','rsq');
if ~isempty(h) && ishghandle(h)
   rsq = 1-(stats.SSresid/stats.SStotal);
   set(h,'String',getString(message('stats:stepwise:footerText_RSQ',...
                                    sprintf('%g',rsq))));
end
h = findobj(F,'tag','adjrsq');
if ~isempty(h) && ishghandle(h)
   rsq = 1-stats.SSresid/stats.SStotal;
   adjrsq = 1 - (1-rsq)*((stats.df0+stats.dfe)/stats.dfe);
   set(h,'String',getString(message('stats:stepwise:footerText_AdjRSQ',...
                                    sprintf('%g',adjrsq))));
end
h = findobj(F,'tag','fstat');
if ~isempty(h) && ishghandle(h)
   set(h,'String',getString(message('stats:stepwise:footerText_F',...
                                    sprintf('%g',stats.fstat))));
end
h = findobj(F,'tag','pval');
if ~isempty(h) && ishghandle(h)
   set(h,'String',getString(message('stats:stepwise:footerText_PVAL',...
                                    sprintf('%g',stats.pval))));
end

% -----------------------------------------
function updatehistoryplot(F,modnums,rmse,in)
%UPDATEHISTORYPLOT Update plot from history variable

% Add a dot for each new model
hHistoryAxes = getappdata(F,'historyaxes');
for j=1:length(modnums)
   line(modnums(j),rmse(j),'Marker','.','MarkerSize',20,...
        'Parent',hHistoryAxes,'Tag','history','UserData',in(j,:),...
        'ButtonDownFcn',{@historycb in(j,:)});
end

% Update x axis limits to include all models
set(hHistoryAxes,'XLim',[.5 modnums(end)+.5]);

% Avoid non-integer x tick labels
set(hHistoryAxes,'XTickMode','auto','XTickLabelMode','auto');
xtick = get(hHistoryAxes,'XTick');
xlabels = get(hHistoryAxes,'XTickLabel');
okticks = (xtick == floor(xtick));
set(hHistoryAxes, 'XTick',xtick(okticks), ...
                  'XTickLabel',deblank(xlabels(okticks,:)));


% -----------------------------------------
function historycb(~,~,in)

% Launch a new stepwise gui starting from this point in the history
F = gcbf;
ptr = get(F,'Pointer');
set(F,'Pointer','watch');

% Get old position (source GUI) and new position (just-created GUI)
y = getappdata(F,'y');
X = getappdata(F,'allx');

newF = stepwise(X,y,in);

% Turn on scaling in new figure if scaling is on in old figure.
if getappdata(F,'doscale')
    doscale([], [], newF);
end

newPos = get(newF,'Position');
oldPos = get(F,'Position');

% Offset new position from old
outerPos = get(F,'OuterPosition');
frameHeight = oldPos(2) - outerPos(2);
titleHeight = outerPos(4) - oldPos(4) - frameHeight;
newPos(1) = oldPos(1) + titleHeight;
newPos(2) = oldPos(2) - titleHeight;

% Make sure it is on-screen
oldu = get(0,'Units');
set(0,'Units','points');
screenSize = get(0,'ScreenSize');
set(0,'Units',oldu);
if newPos(1)+newPos(3) > screenSize(3)-titleHeight
   newPos(1) = titleHeight;
end
if (newPos(2) < titleHeight)
   newPos(2) = screenSize(4) - newPos(4) - titleHeight;
end

set(newF,'Position',newPos);

set(F,'Pointer',ptr);

% -----------------------------------------
function addvarmenu(~,~,F)
%ADDVARMENU Callback for menu item to produce added variable plot

% Prompt for variable to plot
[~,~,in] = getdata(F);
nvars = length(in);
namelist = cell(1,nvars);
for i = 1:nvars
    namelist{i} = getString(message('stats:stepwise:addvar_ListString_Default', ...
        strjust(int2str(i),'left')));
end

[vnum,ok] = listdlg('Name',getString(message('stats:stepwise:addvar_Dlg_Title')),...
                    'ListString',namelist,'SelectionMode','single', ...
                    'ListSize', [400, 300]);

% Create the plot
if ok
   stepaddvarplot(vnum,F);
end

% -----------------------------------------
function stepaddvarplot(vnum,F)
%STEPADDVARPLOT Produce added variable plot

if nargin<2, F = gcf; end

[X,y,in] = getdata(F);
stats = getappdata(F,'stats');

fig = figure('NumberTitle', 'off', ...
             'Name', getString(message('stats:stepwise:addvar_Plot_Title')), ...
             'IntegerHandle', 'off', ...
             'Tag', 'addvarFigure');

% Restore NaN rows if any
wasnan = getappdata(F,'wasnan');
stats.wasnan = wasnan;

addedvarplot(X,y,vnum,in,stats,fig);

% -----------------------------------------
function doscale(hMenu,~,F)
%DOSCALE Toggle "scale inputs" state

% Update saved boolean
oldstate = getappdata(F,'doscale');
newstate = ~oldstate;
setappdata(F,'doscale',newstate);

% Check or un-check menu
if newstate
   newval = 'on';
else
   newval = 'off';
end

% When stepwise is called from the history and scaling is on, doscale
% is called with an empty hMenu.
if isempty(hMenu)
    hMenu = findall(F, 'type', 'uimenu', 'Tag', 'stepwiseMenuScale');    
end

set(hMenu,'Checked',newval);


% Perform new fit without adding it to the history
[X,y,in] = getdata(F);
in = stepmain(X,y,in,F);

% Update plot title to indicate scaling status
hMainAxes = getappdata(F,'plotaxes');
hTitle = get(hMainAxes,'Title');
if newstate
   ttl = getString(message('stats:stepwise:mainAxesTitle_Scaled'));
else
   ttl = getString(message('stats:stepwise:mainAxesTitle_Default'));
end
set(hTitle,'String',ttl);

% Update advice to refresh highlighting
updateadvice(F,in)


% -----------------------------------------
function [X,y,in] = getdata(F)
%GETDATA Get saved data arrays from appdata

% Get response data and boolean in/out flags
y = getappdata(F,'y');
in = getappdata(F,'in');

% Get X data
if getappdata(F,'doscale')
   % Retrieve or calculate scaled version of X
   X = getappdata(F,'scaledx');
   if isempty(X)
      X = getappdata(F,'allx');
      sx = std(X);
      X = X./sx(ones(size(X,1),1),:);
      setappdata(F,'scaledx',X);
   end
else
   X = getappdata(F,'allx');
end

% -----------------------------------------
function stepbutton(~,~,action)
%STEPBUTTON Call-back for the two step buttons

F = gcbf;
switch(action)
 case 'one'
   termnum = getappdata(F,'advisedstep');
   if termnum>0
      dostep(F,termnum);
   end

 case 'all'
   doallsteps(F)
end


% -----------------------------------------
function doallsteps(F)
%DOALLSTEPS Start automatic stepwise regression from the current state

[X,y,in] = getdata(F);

% Perform series of steps
penter = getappdata(F,'penter');
premove = getappdata(F,'premove');
[~,~,~,in,stats,nextstep,history] = ...
    stepwisefit(X,y,'inmodel',in,'penter',penter,'premove',premove,...
                'display','off');

% If any steps were taken, update stored information and figure
inmat = history.in;
if size(inmat,1)>0
   % Remember current state
   setappdata(F,'in',in);
   setappdata(F,'stats',stats);
   setappdata(F,'advisedstep',nextstep);

   % Update all parts of figure
   updateplot(F,stats,in);
   updatefooter(F,stats);
   resizefig(F,[],length(in));
   addhistorystep(F,history.in,history.rmse,history.df0);
   updateadvice(F,in);
end

%-----------------------------------------------------
function savebutton(~,ignore, F)
%SAVEBUTTON Save to workspace

% labels =  {'Coefficients', 'Confidence Intervals', 'Terms In' 'Terms Out' ...
%            'Statistics',   'Coefficient Table',    'History'};
labels = { getString(message('stats:stepwise:ExportLabel_Coef')), ...
           getString(message('stats:stepwise:ExportLabel_CI')), ...
           getString(message('stats:stepwise:ExportLabel_TermsIn')), ...
           getString(message('stats:stepwise:ExportLabel_TermsOut')), ...
           getString(message('stats:stepwise:ExportLabel_Stats')), ...
           getString(message('stats:stepwise:ExportLabel_CoefTable')), ...
           getString(message('stats:stepwise:ExportLabel_History')), ...
         };
     
varnames = {'beta', 'betaci', 'in' ,'out', 'stats', 'coeftab', 'history'};

in = getappdata(F,'in');

p = length(in);
notin = 1:p;
inmodel = find(in);
notin(inmodel) = [];

stats = getappdata(F,'stats');
b = stats.B;
seb = stats.SE;
df = stats.dfe;
penter = getappdata(F, 'penter');
sigprob = 1 - penter/2;
crit = tinv(sigprob,df);

beta = zeros(p,1);
delta = zeros(p,1);
if ~isempty(inmodel)
   delta(inmodel) = seb(inmodel)*crit;
   beta(inmodel) = b(inmodel);
end

% Create structure holding statistics under coefficient display
outst.intercept = stats.intercept;
outst.rmse      = stats.rmse;
outst.rsq       = 1-stats.SSresid/stats.SStotal;
outst.adjrsq    = 1 - (1-outst.rsq)*((stats.df0+stats.dfe)/stats.dfe);
outst.fstat     = stats.fstat;
outst.pval      = stats.pval;

% Create structure holding coefficient table
coeftab.beta  = stats.B;
coeftab.se    = stats.SE;
coeftab.tstat = stats.TSTAT;
coeftab.pval  = stats.PVAL;

% Create structure of history information
H = getappdata(F,'history');
hist.rmse  = H(1,:)';
hist.nvars = H(2,:)';
hist.in    = H(2:end,:)';

items = {beta, [beta-delta beta+delta], inmodel, notin, outst, coeftab, hist};
export2wsdlg(labels, varnames, items, ...
    getString(message('stats:stepwise:ExportTitle')));

resizefig(F,ignore,p);  % needed to fix resize issues cased by above call

%-----------------------------------------------------
function histtipsonoff(varargin)
%HISTIPSONOFF Turn History Tips on or off

h = gcbo;
F = gcbf;
ischecked = get(h,'Checked');
if isequal(ischecked,'on')
   set(h,'Checked','off')
   set(F,'WindowButtonMotionFcn','');
else
   set(h,'Checked','on')
   set(F,'WindowButtonMotionFcn',@histtips);
end

%-----------------------------------------------------
function histtips(varargin)
%HISTTIPS Display fit tips for models on history graph

msg = '';
eventData = varargin{2};
h = eventData.HitObject;
F = gcbf;
if ~isequal(F,ancestor(h,'figure'))
    % For some reason the figure under the cursor is not the one
    % that launched us, so don't try to change the history tips.
    return
end
hTip = findobj(F,'tag','histtip');
if ~isempty(h) && isequal(get(h,'Tag'),'history') ...
               && isequal(get(h,'Type'),'line')
   % Figure out if the cursor is on something we know how to label
   ax = get(h,'Parent');
   x = get(h,'XData');
   y = get(h,'YData');

   % See if the previous tip (if any) is still valid
   if ~isempty(hTip) && ishghandle(hTip)
      oldx = get(hTip,'UserData');
      if (x == oldx)
         return
      end
   end

   % No, create a new tip
   xlim = get(ax,'XLim');
   dx = diff(xlim) * 0.02;
   in = get(h,'UserData');
   if ~any(in)
      msg = 'None';
   else
      msg = sprintf('X%d, X%d, X%d, X%d, X%d,\n',find(in));
      if msg(end)==sprintf('\n') %#ok<SPRINTFN> 
         msg(end-1:end) = [];   % remove ',\n'
      else
         msg(end-2:end) = [];   % remove ', X'
      end
   end
   msg = getString(message('stats:stepwise:HistoryDataTip',msg));
end

% If we can't label this thing, delete the label components
if isempty(msg)
   if ~isempty(hTip)
      delete(hTip);
   end
   return
end

% Otherwise we need to create the proper label
% Create the text object if missing
if isempty(hTip) || ~ishghandle(hTip)
   yellow = [1 1 .85];
   hTip = text(x,y,'','Color','k','VerticalAlignment','bottom',...
                       'Parent',ax, 'Tag','histtip','HitTest','off','PickableParts','none',...
                       'Backgroundcolor',yellow,'Margin',3,...
                       'Edgecolor','k');
end

% Position the text so it is not clipped, then write the label
if (x<sum(xlim)/2)
   ha = 'left';
else
   ha = 'right';
   dx = - dx;
end

set(hTip,'Position',[x+dx y 0],'String',msg,'UserData',x,...
          'HorizontalAlignment',ha,'VerticalAlignment','bottom');

%-----------------------------------------------------
function helpCallback(~, ~)
% display help

fname = [docroot '/toolbox/stats/stats.map'];
if exist(fname,'file')
    helpview(fname,'stepwise_demo')
else
    disp('Unable to display help for Stepwise');
end

%-----------------------------------------------------
function adjustmenu(F)
%ADJUSTMENU Adjust contents of stepwise regression figure menus

% Create a new menu with stepwise-related controls
hStepwise = uimenu(F, ...
    'Label',getString(message('stats:stepwise:FigMenu_Stepwise')), ...
    'Tag','stepwiseMenuStepwise','Position',6);
uimenu(hStepwise, ...
    'Label',getString(message('stats:stepwise:StepwiseMenu_Scale')), ...
    'Checked','off',...
    'Callback',{@doscale F}, 'Tag','stepwiseMenuScale');
uimenu(hStepwise, ...
    'Label',getString(message('stats:stepwise:StepwiseMenu_AddVar')), ...
    'Callback',{@addvarmenu F}, 'Tag','stepwiseMenuAddedVar');
uimenu(hStepwise, ...
    'Label',getString(message('stats:stepwise:StepwiseMenu_HistoryTips')), ...
    'Checked','on',...
    'Callback',{@histtipsonoff}, 'Tag','stepwiseMenuHistTips');
setappdata(F,'doscale',false);

% Remove some menus entirely
h = findall(F, 'Type','uimenu', 'Parent',F);
h0 = findall(h,'flat', 'Tag','figMenuInsert');
if (~isempty(h0))
   delete(h0);
   h(h==h0) = [];
end
h0 = findall(h,'flat', 'Tag','figMenuView');
if (~isempty(h0))
   delete(h0);
   h(h==h0) = [];
end

% Add or remove some items from other menus
% Fix Tools menu
h0 = findall(h,'flat', 'Tag','figMenuTools');
h1 = findall(h0, 'Type','uimenu', 'Parent',h0);
for j=length(h1):-1:1
   mtag = get(h1(j),'Tag');
   if ~ismember(mtag, {'figMenuZoomOut' 'figMenuZoomIn'})
     delete(h1(j));
     h1(j) = [];
   end
end

% Fix Edit menu
% On PC, leave only Copy Figure and Copy Option
% On other platforms, delete the entire Edit menu
h0 = findall(h,'flat', 'Tag','figMenuEdit');
if ispc
    h1 = findall(h0, 'Type','uimenu', 'Parent',h0);
    for j=length(h1):-1:1
       mtag = get(h1(j),'Tag');
       if ~ismember(mtag, {'figMenuEditCopyOptions' 'figMenuEditCopyFigure'})
          delete(h1(j));
          h1(j) = [];
       else
          set(h1(j),'Separator','off');
       end
    end
else
   delete(h0);
   h(h==h0) = [];
end

% Fix File menu
h0 = findall(h,'flat', 'Tag','figMenuFile');
delete(findall(h0, 'Tag', 'figMenuMFileSaveAs', 'Parent',h0));

% Add Stepwise help to help menu (add a separator to old "first" item)
helpMenu = findall(h,'flat', 'Tag','figMenuHelp');
firstItem = findall(helpMenu, 'Parent', helpMenu, 'Position', 1);
set(firstItem, 'Separator', 'on');
uimenu(helpMenu, ...
    'Label',getString(message('stats:stepwise:FigMenu_Help')), ...
    'Position', 1, ...
    'Callback', {@helpCallback}, 'Tag','stepwiseMenuHelp');

