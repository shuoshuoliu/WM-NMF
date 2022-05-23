function [h,ax,BigAx] = gplotmatrix(varargin)
%GPLOTMATRIX  Scatter plot matrix with grouping variable.
%   GPLOTMATRIX(X,Y,G) creates a matrix of scatter plots of the columns of
%   X against the columns of Y, grouped by G.  If X is P-by-M and Y is
%   P-by-N, GPLOTMATRIX will produce a N-by-M matrix of axes.  If you omit
%   Y or specify it as [], the function graphs X vs. X.  G is a grouping
%   variable that determines the marker and color assigned to each point in
%   each matrix, and it can be a categorical variable, vector, string
%   array, string matrix, or cell array of strings.  Alternatively G can be
%   a cell array of grouping variables (such as {G1 G2 G3}) to group the
%   values in X by each unique combination of grouping variable values.
%
%   Use the data cursor to read precise values from the plot, as well as
%   the observation number and the values of related variables.
%
%   GPLOTMATRIX(X,Y,G,CLR,SYM,SIZ) specifies the colors, markers, and size
%   to use.  CLR is a string of color specifications, and SYM is a string
%   of marker specifications.  Type "help plot" for more information.  For
%   example, if SYM='o+x', the first group will be plotted with a circle,
%   the second with plus, and the third with x. SIZ is a marker size to use
%   for all plots.  By default, the colors are 'bgrcmyk', the marker is
%   '.', and the marker size depends on the number of plots and the size of
%   the figure window.
%
%   GPLOTMATRIX(X,Y,G,CLR,SYM,SIZ,DOLEG) lets you control whether legends
%   are created.  Set DOLEG to 'on' (default) or 'off'.
%
%   GPLOTMATRIX(X,Y,G,CLR,SYM,SIZ,DOLEG,DISPOPT) lets you control how to
%   fill the diagonals in a plot of X vs. X.  Set DISPOPT to 'none' to
%   leave them blank, 'hist' to plot histograms of all data points,
%   'stairs'(default if there is more than one group) to display the
%   outlines of grouped histograms, 'grpbars' to plot grouped histogram
%   bars, or 'variable' to write the variable names.
%
%   GPLOTMATRIX(X,Y,G,CLR,SYM,SIZ,DOLEG,DISPOPT,XNAM,YNAM) specifies XNAM
%   and YNAM as the names of the X and Y variables.  Each must be a
%   character array, string array or cell array of strings of the
%   appropriate dimension.
%
%   GPLOTMATRIX(PARENT,...) plots into the parent specified as either a
%   figure or a uipanel.
%
%   [H,AX,BigAx] = GPLOTMATRIX(...) returns an array of handles H to the
%   plotted points; a matrix AX of handles to the individual subaxes; and a
%   handle BIGAX to big (invisible) axes framing the subaxes.  The third
%   dimension of H corresponds to groups in G.  If DISPOPT is 'hist',
%   'stairs' or 'grpbars', AX contains one extra row of handles to
%   invisible axes in which the histograms are plotted. BigAx is left as
%   the CurrentAxes so that a subsequent TITLE, XLABEL, or YLABEL will be
%   centered with respect to the matrix of axes.
%
%   Example:
%      load carsmall; X =
%      [MPG,Acceleration,Displacement,Weight,Horsepower]; varNames = {'MPG'
%      'Acceleration' 'Displacement' 'Weight' 'Horsepower'};
%      gplotmatrix(X,[],Cylinders,'bgrcm',[],[],'on','hist',varNames);
%
%   See also GRPSTATS, GSCATTER, PLOTMATRIX.

%   Copyright 1993-2019 The MathWorks, Inc.



narginchk(1,11);
nin = nargin;
if nin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end 
% Allow uipanel/figure input as the first argument (parent)
if(isa(varargin{1},'matlab.ui.Figure') ||...
        isa(varargin{1},'matlab.ui.container.Panel'))
    parent = varargin{1};
    varargin = varargin(2:end);
    nin = nin-1;
else
    % Check if a figure exists to be plot into, or create one
    parent = ancestor(gca,'figure');
end

% Use inputParser to parse inputs
iP = inputParser;
addOptional(iP,'x',[],@(x) isnumeric(x) || isdatetime(x) || isduration(x));
addOptional(iP,'y',[],@(x) isnumeric(x) || isdatetime(x) || isduration(x));
addOptional(iP,'g',[],@(x) iscategorical(x)||isnumeric(x)||iscell(x)||...
    ischar(x)||islogical(x));
addOptional(iP,'clr',[],@(x) ischar(x)||isnumeric(x)||iscell(x));
addOptional(iP,'symb',[],@(x) ischar(x)||isnumeric(x)||iscell(x));
addOptional(iP,'siz',[],@isnumeric);
addOptional(iP,'doleg',[],@(x) ischar(x)||isnumeric(x) ||islogical(x));
addOptional(iP,'dispopt',[],@(x) ischar(x)||isnumeric(x));
addOptional(iP,'xnam',[],@(x) ischar(x)||isnumeric(x)||iscell(x));
addOptional(iP,'ynam',[],@(x) ischar(x)||isnumeric(x)||iscell(x));
parse(iP,varargin{:});

% Assign parsed inputs
x = iP.Results.x;
y = iP.Results.y;
g = iP.Results.g;
clr = iP.Results.clr;
sym = iP.Results.symb;
siz = iP.Results.siz;
doleg = iP.Results.doleg;
dispopt = iP.Results.dispopt;
xnam = iP.Results.xnam;
ynam = iP.Results.ynam;

if (nin < 2), y = []; end
if isempty(y) % gplotmatrix(x)  
  rows = size(x,2); cols = rows;
  y = x;
  XvsX = true;
else % gplotmatrix(x,y)
  rows = size(y,2); cols = size(x,2);
  XvsX = false;
end
if (nin > 2) && ~isempty(g)
    if iscategorical(g) 
        g = removecats(g);
    end
    [g,gn] = mgrp2idx(g,size(x,1),',');
    ng = max(g);  
else
   % Set g to be a vector of ones instead of empty so that datatip callback
   % can index into it
   g = ones(size(x,1),1);
   gn = [];
   ng = 1;
end

% Default colors, markers, etc.
if (nin < 4) || isempty(clr), clr = 'bgrcmyk'; end
if isnumeric(clr) && isrow(clr)
    if length(clr) == 3
        % if it is a row vector, replicate it to a 2-row matrix and pass to
        % internal.stats.cycleLineProperties
        clr = [clr;clr];
    else
        error(message('stats:internal:colorStringToRGB:ValueMustBe3ElementVector'));
    end
end
clr = internal.stats.cycleLineProperties(ng,clr);
clr = internal.stats.colorStringToRGB(clr);

if (nin < 5) || isempty(sym), sym = '.'; end
if (nin < 6), siz = []; end
if (nin < 7) || isempty(doleg), doleg = 'on'; end
if (nin < 8) || isempty(dispopt)
    if ng > 1
        dispopt = 's';
    else
        dispopt = 'g';       
    end
end
if (nin < 9) || isempty(xnam)
   xnam = {};
else
   if ischar(xnam) && (size(xnam,1)==cols)
       xnam = cellstr(xnam);
   elseif iscellstr(xnam) && (numel(xnam)==cols)
      % ok
   else
      error(message('stats:gplotmatrix:XnamSizeMismatch'));
   end
end
if (XvsX)
   ynam = xnam;
elseif (nin < 10) || isempty(ynam)
   ynam = {};
else
   if ischar(ynam) && (size(ynam,1)==rows)
       ynam = cellstr(ynam);
   elseif iscellstr(ynam) && (numel(ynam)==rows)
      % ok
   else
      error(message('stats:gplotmatrix:YnamSizeMismatch'));
   end
end

dispopt = internal.stats.getParamVal(dispopt,{'hist','stairs','grpbars','variable','none'},'DISPOPT');
% What should go into the plot matrix?
doleg = internal.stats.parseOnOff(doleg,'''Legend''');
doleg = (doleg==1) && (~XvsX || (rows>1)) && ~isempty(gn);
dohist = XvsX && (dispopt(1)=='h');
doghistBar = (XvsX && (dispopt(1)=='g'));
doghistStair = (XvsX && (dispopt(1)=='s'));
donames = (XvsX && (dispopt(1)=='v'));

% Don't plot anything if either x or y is empty
if isempty(rows) || isempty(cols)
   if nargout>0, h = []; ax = []; BigAx = []; end
   return
end

if ~ismatrix(x) || ~ismatrix(y)
   error(message('stats:gplotmatrix:MatrixRequired'));
end
if size(x,1)~=size(y,1)
  error(message('stats:gplotmatrix:XYSizeMismatch'));
end
if (~isempty(g)) && (length(g) ~= size(x,1))
  error(message('stats:gplotmatrix:XGSizeMismatch'));
end

% Check that x is numeric before passing to histogram
histNorm = 'pdf';
if (isdatetime(x) || isduration(x)) && (dohist || doghistBar || doghistStair)
    histNorm = 'count';
end

% Create/find BigAx and make it invisible
if(isa(parent,'matlab.ui.Figure'))
    clf(parent);
    BigAx = newplot;
else % Parent is a uipanel
    % To maintain the behavior for figures: Delete all children of panel
    % and create a new axes
    delete(parent.Children);
    BigAx = axes;
end
hold_state = ishold;
set(BigAx,'Visible','off','color','none','Parent',parent)

if (isempty(siz))
   siz = repmat(get(0,'defaultlinemarkersize'), size(sym));
   if any(sym=='.')
      units = get(BigAx,'units');
      set(BigAx,'units','pixels');
      pos = get(BigAx,'Position');
      set(BigAx,'units',units);
      siz(sym=='.') = max(1,min(15, ...
                       round(15*min(pos(3:4))/size(x,1)/max(rows,cols))));
   end
end

% Store global data for datatips into BixAx
ginds = cell(1,ng);
for i=1:ng
    ginds{i} = find(g==i);
end

setappdata(BigAx,'ginds',ginds);
setappdata(BigAx,'xnam',xnam);
setappdata(BigAx,'ynam',ynam);
setappdata(BigAx,'x',x);
setappdata(BigAx,'y',y);
setappdata(BigAx,'XvsX',XvsX);
setappdata(BigAx,'gn',gn);
setTickChanged(BigAx,false);

% Make datatips show up in front of axes
dcm_obj = datacursormode(ancestor(BigAx,'figure'));

dataCursorBehaviorObj = hgbehaviorfactory('DataCursor');
set(dataCursorBehaviorObj,'UpdateFcn',@gplotmatrixDatatipCallback);

% Create and plot into axes
ax2filled = false(max(rows,cols),1);
ax2 = gobjects(size(ax2filled));
pos = get(BigAx,'Position');
width = pos(3)/cols;
height = pos(4)/rows;
space = .02; % 2 percent space between axes
pos(1:2) = pos(1:2) + space*[width height];
[m,n,k] = size(y); %#ok<ASGLU>
xlim = repmat(cat(3,min(x,[],1),max(x,[],1)),[rows 1 1]);
ylim = repmat(cat(3,min(y,[],1)',max(y,[],1)'),[1 cols 1]);

for i=rows:-1:1
   for j=1:1:cols
      axPos = [pos(1)+(j-1)*width pos(2)+(rows-i)*height ...
               width*(1-space) height*(1-space)];
      ax(i,j) = axes('Position',axPos, 'visible', 'on', 'Box','on','Parent',parent);

      if ((i==j) && XvsX)
          if (dohist||doghistBar||doghistStair)
              histax = axes('Position',axPos,'Parent',parent);
              ax2(j) = histax;
              ax2filled(j) = true;
              if dohist
                  hhdiag(i,1) = histogram(reshape(y(:,i,:),[m k]),'BinMethod','scott','DisplayStyle','bar','Norm',histNorm);
                  if ng == 1
                      set(hhdiag(i),'FaceColor',clr(1,:));
                  end
              elseif doghistBar
                  hhdiag(i,:) = internal.stats.plotGroupedHist(reshape(y(:,i,:),[m k]),g,'Color',clr,'DisplayStyle','bar','Norm',histNorm);
              elseif doghistStair
                  hhdiag(i,:) = internal.stats.plotGroupedHist(reshape(y(:,i,:),[m k]),g,'Color',clr,'DisplayStyle','stairs','Norm',histNorm);
              end
              % Additionally scatter plot the data in the axes that is used
              % to display the tick-labels. % Set Marker to 'n' for none
              internal.stats.scatterplot(ax(i,j),reshape(x(:,j,:),[m k]), ...
                  reshape(y(:,i,:),[m k]), g, clr, 'n', siz);
              set(histax, 'YAxisLocation', 'right', ...
                  'Visible','off', 'XTick',[], 'YTick',[], ...
                  'XGrid','off', 'YGrid','off', ...
                  'XTickLabel','', 'YTickLabel','');
              axis(histax,'tight');
              xlim(i,j,:) = get(histax,'xlim');
              set(histax, 'ylim', get(histax,'ylim').*[1 1.05]);
          end
      else
         hhij = internal.stats.scatterplot(ax(i,j),reshape(x(:,j,:),[m k]), ...
                              reshape(y(:,i,:),[m k]), ...
                              g, clr, sym, siz);
         hh(i,j,1:length(hhij)) = hhij;
         axis(ax(i,j),'tight');
         ylim(i,j,:) = get(ax(i,j),'ylim');
         xlim(i,j,:) = get(ax(i,j),'xlim');

         % Store information for gname
         set(ax(i,j), 'UserData', {'gscatter' x(:,j,:) y(:,i,:) g});
         % Attach data cursor
         for q=1:ng
             hgaddbehavior(hh(i,j,q),dataCursorBehaviorObj);
             setappdata(hh(i,j,q),'dtcallbackdata',{BigAx,q,i,j});
         end
      end
      set(ax(i,j),'xlimmode','auto', 'ylimmode','auto', ...
                  'xgrid','off', 'ygrid','off')
   end
end

% Fill in histogram handles
if XvsX && (dohist || doghistBar || doghistStair)
    for i=1:rows
        hh(i,i,:) = hhdiag(i,:);
    end
end

xlimmin = min(xlim(:,:,1),[],1); xlimmax = max(xlim(:,:,2),[],1);
ylimmin = min(ylim(:,:,1),[],2); ylimmax = max(ylim(:,:,2),[],2);

% Set all the limits of a row or column to be the same and leave just a 5%
% gap between data and axes.
inset = .05;
if isnumeric(ylimmin)
    for i=1:rows
        set(ax(i,1),'ylim',[ylimmin(i,1) ylimmax(i,1)])
        dy = diff(get(ax(i,1),'ylim'))*inset;
        set(ax(i,:),'ylim',[ylimmin(i,1)-dy ylimmax(i,1)+dy])
    end
end

if isnumeric(xlimmin)
    for j=1:cols
        set(ax(1,j),'xlim',[xlimmin(1,j) xlimmax(1,j)])
        dx = diff(get(ax(1,j),'xlim'))*inset;
        set(ax(:,j),'xlim',[xlimmin(1,j)-dx xlimmax(1,j)+dx])
        if ax2filled(j)
            set(ax2(j),'xlim',[xlimmin(1,j)-dx xlimmax(1,j)+dx])
        end
    end
end

% Label plots one way or the other
figparent = ancestor(parent,'figure');
if (donames && ~isempty(xnam))
   for j=1:cols
      set(figparent,'CurrentAxes',ax(j,j));
      h = text(0.5,0.5, -.1,...
          xnam{j},'HorizontalAlignment','center',...
          'VerticalAlignment','middle','Units','normalized');
   end
else
   if ~isempty(xnam)
      for j=1:cols, xlabel(ax(rows,j),xnam{j}); end
   end
   if ~isempty(ynam)
      for i=1:rows, ylabel(ax(i,1),ynam{i}); end
   end
end

% Ticks and labels on outer plots only
set(ax(1:rows-1,:),'xticklabel','')
set(ax(:,2:cols),'yticklabel','')
set(BigAx,'userdata',ax,'tag','PlotMatrixBigAx');

setTickChanged(BigAx,false);

addlistener(ax(1:rows,1),   'YTickLabel','PostSet',@(varargin)setTickChanged(BigAx,true));
addlistener(ax(rows,1:cols),'XTickLabel','PostSet',@(varargin)setTickChanged(BigAx,true));

% Set XTick and YTick only if data is numeric, this is automatically set
% for datetime arrays afterwards
if(isnumeric(get(ax(rows,1),'xtick')))
    set(BigAx,'XTick',get(ax(rows,1),'xtick'));
end
if(isnumeric(get(ax(rows,1),'ytick')))
    set(BigAx,'YTick',get(ax(rows,1),'ytick'));
end          

% Create legend if requested; base it on the top right plot
if (doleg)
   gn = gn(ismember(1:size(gn,1),g),:);
   legend(ax(1,cols),gn);
end

% Make BigAx the CurrentAxes
set(figparent,'CurrentAx',BigAx)
if ~hold_state
   set(figparent,'NextPlot','replace')
end

% Also set Title and X/YLabel visibility to on and strings to empty
set([get(BigAx,'Title'); get(BigAx,'XLabel'); get(BigAx,'YLabel')], ...
 'String','','Visible','on')

for i=1:cols
    axToBeLinked = ax(:,i);
    hz = zoom;
    if any(ax2filled)
        hlinkx(i) = linkprop([axToBeLinked;ax2(i)],{'XLim','XScale'});
        setAxesZoomMotion(hz,ax2(i),'horizontal');        
    else
        hlinkx(i) = linkprop(axToBeLinked,{'XLim','XScale'});        
    end
end
for j=1:rows
    axToBeLinked = ax(j,:);
    hlinky(j) = linkprop(axToBeLinked,{'YLim','YScale'});
end
setappdata(BigAx,'LinkPropXLim',hlinkx);
setappdata(BigAx,'LinkPropYLim',hlinky);

Fcallback = @(varargin)localLabelCallback(BigAx,ax,ax2,ax2filled,rows,cols,varargin{:});
localLabelCallback(BigAx,ax,ax2,ax2filled,rows,cols)
set(figparent,'SizeChanged',Fcallback);
hp = pan;
set(hp,'ActionPreCallback',@(varargin)localSizeChangedCallback(BigAx,ax,rows,cols,varargin{:}));
set(hz,'ActionPostCallback',Fcallback);
set(hp,'ActionPostCallback',Fcallback);

setTickChanged(BigAx,false);

if nargout~=0
  h = hh;
  if any(ax2filled)
     ax = [ax; ax2(:)'];
  end
end

% The following nested functions exist only to support loading of figures
% that included these callbacks
    function gplotmatrixSizeChangedCallback(varargin) %#ok<DEFNU> 
        localSizeChangedCallback(BigAx,ax,rows,cols,varargin{:})
    end
    function gplotmatrixLabelCallback(varargin) %#ok<DEFNU> 
        localLabelCallback(BigAx,ax,ax2,ax2filled,rows,cols,varargin{:})
    end
end

% --- remember if ticks change externally
function setTickChanged(BigAx,tf)
    setappdata(BigAx,'TickChanged',tf);
end
function tf = getTickChanged(BigAx)
    tf = getappdata(BigAx,'TickChanged');
    if isempty(tf)
        tf = false;
    end
end

% -----------------------------
function datatipTxt = gplotmatrixDatatipCallback(obj,evt)

target = get(evt,'Target');
ind = get(evt,'DataIndex');
pos = get(evt,'Position');

dtcallbackdata = getappdata(target,'dtcallbackdata');
[BigAx,gnum,row,col] = dtcallbackdata{:};

ginds = getappdata(BigAx,'ginds');
xnam = getappdata(BigAx,'xnam');
ynam = getappdata(BigAx,'ynam');
xdat = getappdata(BigAx,'x');
ydat = getappdata(BigAx,'y');
XvsX = getappdata(BigAx,'XvsX');
gn = getappdata(BigAx,'gn');
% If NumGroups(gn) is empty, set them as {1}, to be indexed into later
if isempty(gn)
    gn = {1};
end
gind = ginds{gnum};
obsind = gind(ind);

xvals = xdat(obsind,:);
yvals = ydat(obsind,:);

x = xvals(col);
y = yvals(row);
[xNum,yNum] = matlab.graphics.internal.makeNumeric(obj,x,y);

if xNum~=pos(1) || yNum~=pos(2)
    % Something is inconsistent, display default datatip.
    datatipTxt = {sprintf('X: %s',convertToString(pos(1))),sprintf('Y: %s',convertToString(pos(2)))};
else
    if isempty(xnam)
        xnam = cell(size(xdat,2),1);
        for i = 1:size(xdat,2)
            xnam{i} = getString(message('stats:gplotmatrix:VarLabelXvar',convertToString(i)));
        end
    end
    if isempty(ynam)
        ynam = cell(size(ydat,2),1);
        for i = 1:size(ydat,2)
            ynam{i} = getString(message('stats:gplotmatrix:VarLabelYvar',convertToString(i)));
        end
    end

    % Generate datatip text.
    datatipTxt = {
        [xnam{col},': ',convertToString(x)],...
        [ynam{row},': ',convertToString(y)],...
        '',...
        getString(message('stats:gplotmatrix:LabelObservation',convertToString(obsind))),...
        };

    if ~isempty(gn)
        datatipTxt{end+1} = getString(message('stats:gplotmatrix:LabelGroup', gn{gnum}));
    end
    datatipTxt{end+1} = '';

    xnamTxt = cell(length(xvals),1);
    for i=1:length(xvals)
        xnamTxt{i} = [xnam{i} ': ' convertToString(xvals(i))];
    end
    datatipTxt = {datatipTxt{:}, xnamTxt{:}};
    
    if ~XvsX
        ynamTxt = cell(length(yvals),1);
        for i=1:length(yvals)
            ynamTxt{i} = [ynam{i} ': ' convertToString(yvals(i))];
        end
        datatipTxt = {datatipTxt{:}, ynamTxt{:}};
    end

end
end

% -----------------------------
function localLabelCallback(BigAx,ax,ax2,ax2filled,rows,cols,~,~)
if ~all(isgraphics(ax(:)))
    return;
end

% Check for whether the figure may have been loaded from a file without
% callbacks. If so, restore them here.
if ~isempty(gcbf)
    % Make sure hz1 and ax2 have same figure parent. If axes are filled
    % then fig is ancestor of axes, else fig is either gcbf or gcf.
    if any(ax2filled)
        ax2fill = ax2(ax2filled);
        fig = ancestor(ax2fill(1),'figure');
    else
        fig = [];
    end
    if isempty(fig)
        if isempty(get(groot,'CurrentFigure')) || get(groot,'CurrentFigure') == gcbf
            fig = gcbf;
        else
            fig = get(groot,'CurrentFigure');
        end
    end
    hp1 = pan(fig);
    if isempty(get(hp1,'ActionPreCallback'))
        hz1 = zoom(fig);
        set(hp1,'ActionPreCallback',@(varargin)localSizeChangedCallback(BigAx,ax,rows,cols,varargin{:}));
        set(hz1,'ActionPostCallback',@(varargin)localLabelCallback(BigAx,ax,ax2,ax2filled,rows,cols,varargin{:}));
        set(hp1,'ActionPostCallback',@(varargin)localLabelCallback(BigAx,ax,ax2,ax2filled,rows,cols,varargin{:}));
        if any(ax2filled)
            setAxesZoomMotion(hz1,ax2,'horizontal');        
        end
    end
end

TickChanged = getTickChanged(BigAx);
if ~TickChanged
    set(ax(1:rows,1),'YTickLabelMode','auto');
    set(ax(rows,1:cols),'XTickLabelMode','auto');
    
    if rows>=2
        for ii = 1:rows
            rangeY = diff(ax(ii,1).YLim);
            if ax(ii,1).YTick(1)- ax(ii,1).YLim(1) < rangeY*0.05 && ii~=rows
                ax(ii,1).YTickLabel{1} = '';
            end
            if ax(ii,1).YLim(2)- ax(ii,1).YTick(end)< rangeY*0.05 && ii~=1
                ax(ii,1).YTickLabel{end} = '';
            end
        end
        for jj = 1:cols
            htx1 = text(0,0,ax(rows,jj).XTickLabel(1),'Units','characters','Parent',ax(rows,jj));
            dx1 = htx1.Extent(3);
            delete(htx1);
            htx2 = text(0,0,ax(rows,jj).XTickLabel(end),'Units','characters','Parent',ax(rows,jj));
            dx2 = htx2.Extent(3);
            delete(htx2);
            rangeX = diff(ax(rows,jj).XLim);
            if ax(rows,jj).XTick(1) - ax(rows,jj).XLim(1)< rangeX*dx1*0.02 && jj~=1
                ax(rows,jj).XTickLabel{1} = '';
            end
            if ax(rows,jj).XLim(2)- ax(rows,jj).XTick(end)< rangeX*dx2*0.02 && jj~=cols
                ax(rows,jj).XTickLabel{end} = '';
            end
        end
        
        setTickChanged(BigAx,false); % change was internal
    end
end
end

% -----------------------------
function localSizeChangedCallback(BigAx,ax,rows,cols,~,~) 
    if ~all(isgraphics(ax(:)))
        return;
    end
    set(ax(1:rows,1),'YTickLabelMode','auto');
    set(ax(rows,1:cols),'XTickLabelMode','auto');
    setTickChanged(BigAx,false); % change was internal
end


function str = convertToString(val)
if(isnumeric(val))
    str = num2str(val);
else
    str = char(string(val));
end
end