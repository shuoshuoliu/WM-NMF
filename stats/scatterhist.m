function h = scatterhist(x,y,varargin)
%SCATTERHIST 2D scatter plot with marginal histograms.
%   SCATTERHIST(X,Y) creates a 2D scatterplot of the data in the vectors X
%   and Y, and puts a univariate histogram on the horizontal and vertical
%   axes of the plot.  X and Y must be the same length.
%
%   SCATTERHIST(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
%   parameters and their values to control how the plot is made. Valid
%   parameters are the following:
%
%       'NBins'    - A positive integer value, a two-element vector, or a
%                    cell array of size 1-by-1, 1-by-2, or 2-by-1,
%                    specifying the number of bins for each group in the X
%                    and Y histograms.  All numbers should be positive
%                    integers greater than or equal to 2.
%
%      'Location'  - A string controlling the location of the marginal
%                    histograms within the figure.  'SouthWest' (the
%                    default) plots the histograms below and to the left of
%                    the scatterplot, 'SouthEast' plots them below and to
%                    the right, 'NorthEast' above and to the right, and
%                    'NorthWest' above and to the left.
%
%      'Direction' - A string controlling the direction of the marginal
%                    histograms in the figure. 'in' (the default) plots the
%                    histograms with bars directed in towards the
%                    scatterplot, 'out' plots the histograms with bars
%                    directed out away from the scatterplot.
%
%      'Group'     - A grouping variable. SCATTERHIST(X,Y,'Group',G)
%                    creates a 2D GSCATTER plot instead of a SCATTER plot,
%                    and the marginal histograms are replaced by grouped
%                    histograms.
%
%      'PlotGroup' - Indicates if grouped histograms or grouped kernel
%                    density plots are created when a grouping variable is
%                    given. Choices:
%       'on' or true     - This is the default when a grouping variable is
%                          specified. 
%       'off' or false   - The histograms or kernel density plots of the
%                          whole data set will be created if it is set to
%                          false.
% 
%      'Style'     - Grouped histogram display style, Choices:
%       'stairs'    - display a stairstep plot, which displays the outline
%                     of the histogram without filling the interior. This
%                     is default if there is more than one group. 
%       'bar'       - displays a histogram bar plot. This is default if
%                     there is only one group or when PlotGroup is set to
%                     false.
%
%      'Kernel'    - Indicates if histograms are replaced by grouped kernel
%                    density plots. Choices are:
%       'off'       - This is the default. 
%       'on'        - Plot kernel density instead of histograms or grouped
%                      histograms.
%       'overlay'   - A kernel density plot of the overall marginal
%                      distribution is overlaid on the histogram or grouped
%                      histograms.
%
%      'Color'     - A string or a matrix of RGB values that specifies the
%                    marker color of each group. For
%                    example, 'r' or {'red'} if not grouped, or 'rbk' if
%                    grouped. See the ColorSpec documentation for a table
%                    that lists the predefined colors and their RGB
%                    equivalents.
%
%      'Bandwidth' - The bandwidth of the kernel smoothing window. The
%                    default is empty and SCATTERHIST finds the optimal
%                    bandwidths for estimating normal densities. However,
%                    you may want to choose smaller values to reveal small
%                    features, and it must be a matrix of size 2 by K,
%                    where K is the number of unique groups. The first row
%                    gives the bandwidth of each group in X, and the second
%                    row gives the bandwidth of each group in Y.
%
%      'Legend'    - Indicates if a legend is created when a grouping
%                    variable is given. Choices:
%       'on' or true     - Sets legend visible. This is the default when a
%                          grouping variable is specified.
%       'off' or false   - Sets legend invisible. This is the default when
%                          no grouping variable is specified
%
%      'Parent'    - A uipanel or a figure specifying the parent container 
%                    of the plot.
%
%   The following parameters have effects on kernel density plots only.
%      'LineStyle' - A cell vector of valid line style strings that specify
%                    the 'LineStyle' property of each kernel density line.
%                    For example, {'-',':','-.'}. See documentation of PLOT
%                    for a list of valid line style strings.
%                   
%      'LineWidth' - A scalar or a vector of numbers that specify the
%                    'LineWidth' property of each kernel density line. For
%                    example, [1 2 3].
% 
%      'Marker'    - A string that specifies the markers used for each
%                    group in the gscatter plot. See the documentation of
%                    PLOT for a list of valid markers. The default is 'o',
%                    a circle.
%
%      'MarkerSize'- A numerical vector that specifies the size of the each
%                    marker. The default value is 6.
%
%   NOTE:SCATTERHIST cycles through the available values of the parameters:
%   'LineWidth', 'LineStyle', 'Color', 'Marker' and 'MarkerSize' when the
%   total number of groups is greater than the number of specified values.
%
%
%   H = SCATTERHIST(...) returns a vector of three axes handles for the
%   scatterplot, the plot along the horizontal axis, and the plot along the
%   vertical axis, respectively.
%
%   SCATTERHIST(X,Y,NBINS) is supported for backwards compatibility.
%
%   Example:
%      Independent normal and lognormal random samples
%         x = randn(1000,1);
%         y = exp(.5*randn(1000,1));
%         scatterhist(x,y)
%      Marginal uniform samples that are not independent
%         u = copularnd('Gaussian',.8,1000);
%         scatterhist(u(:,1),u(:,2))
%      Mixed discrete and continuous data
%         load('carsmall');
%         scatterhist(Weight,Cylinders,'NBins',[10 3])
%      Scatter plot with kernel density plots  
%         load fisheriris.mat;
%         x = meas(:,1);
%         y = meas(:,2);
%         scatterhist(x,y,'Group',species,'Kernel','on');
%      Draw the scatterhist in a uipanel
%         load carsmall
%         hp = uipanel('position',[0 .5 .5 .5]);
%         scatterhist(Weight,MPG,'Parent',hp);
%      Scatter plot with grouped histogram plots  
%         load fisheriris.mat;
%         x = meas(:,1);
%         y = meas(:,2);
%         scatterhist(x,y,'Group',species);
%
%   See also gscatter, plot, hist, histfit, grpstats, grp2idx, ksdensity

%   Copyright 2006-2018 The MathWorks, Inc.


narginchk(2, Inf);
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if ~isvector(x) || ~isnumeric(x) || ~isvector(y) || ~isnumeric(y)
    error(message('stats:scatterhist:BadXY'));
end
if numel(x)~=numel(y)
    error(message('stats:scatterhist:BadXY'));
end
x = x(:);
y = y(:);

if nargin == 3
    if nargout == 0
        scatterhist(x,y,'nbins',varargin{1});      
    else
        h = scatterhist(x,y,'nbins',varargin{1});
    end
    return;
else
    % Parse name/value pairs
    pnames = {'nbins','location','direction','group','plotgroup',...
        'style','kernel','color', 'marker', 'markersize',...
        'linestyle','linewidth','legend', 'bandwidth', 'Xfunc', 'Yfunc', 'parent'};
    
    dflts =  {  [],  'SouthWest',  'in',   [], true, ...
        'stair', 'off',  [],   'o'         6        {'-',':','-.','--'}, 2, ...
        true,        [],           []        []         []} ;
    [nbins,location,direction,group,plotGroup,style,kernel,clr,marker, ...
        siz, lStyle, lw,doLegend, ww, Xfunc, Yfunc, parent, setFlag]...
        = internal.stats.parseArgs(pnames, dflts, varargin{:});
    
    % Validate values
    if setFlag.nbins &&...
            (numel(nbins)>2||(~iscell(nbins)&&~internal.stats.isIntegerVals(nbins))...  % test positive integer value or two-element vector
            || (iscell(nbins) && (~isempty(nbins) && ~internal.stats.isIntegerVals(nbins{1})||...  % test 1x1 cell 
                                    (numel(nbins)==2 && ~internal.stats.isIntegerVals(nbins{2}))...   % test 1x2 or 2x1 cell 
                                 )...
               )... 
            )     
        error(message('stats:scatterhist:BadBins'));
    end  
    if isempty(group) && iscell(nbins)
        error(message('stats:scatterhist:BadBinsWithoutGroup'));
    end

    % Validate values         
    styleKWDs = {'stairs','bar'};
    style = internal.stats.getParamVal(style,styleKWDs,'''Style''');
    
    kernelKWDs = {'on','off','overlay'};
    kernel = internal.stats.getParamVal(kernel,kernelKWDs,'''Kernel''');
    
    if ~(isnumeric(siz)&& all(siz>0))
        error(message('stats:scatterhist:BadNumericLineProperty',...
            '''MarkerSize'''));
    end
    
    if ischar(lStyle) && size(lStyle,1) == 1
        lStyle = cellstr(lStyle);
    end
    
    if ~iscellstr(lStyle)
        error(message('stats:scatterhist:BadLineStyle'));
    end
    if ~(isnumeric(lw)&& all(lw>0))
        error(message('stats:scatterhist:BadNumericLineProperty',...
            '''LineWidth'''));
    end
    
    plotGroup = internal.stats.parseOnOff(plotGroup,'''Plotgroup''');
    doLegend = internal.stats.parseOnOff(doLegend,'''Legend''');
end

if ~isempty(group)    
    [gInt,gn,~] = grp2idx(group);
    if numel(gInt) ~= numel(x)
        error(message('stats:scatterhist:BadGroup'));
    end
    k = numel(gn); % number of unique groups   
    if strcmpi(kernel,'overlay') && ~setFlag.linestyle
        lStyle = {'-.'};
    end
    [marker,siz,lStyle,lw] = ...
        internal.stats.cycleLineProperties(k,marker,siz,lStyle,lw);
else
    k = 1;
    doLegend = false;
end


if isempty(clr)
    clr = lines(k);
else
    if isnumeric(clr) && isrow(clr)
        if length(clr) == 3
        % if it is a row vector, replicate it to a 2-row matrix and pass to
        % internal.stats.cycleLineProperties
            clr = [clr;clr];
        else
            error(message('stats:internal:colorStringToRGB:ValueMustBe3ElementVector'));
        end
    end
    clr = internal.stats.cycleLineProperties(k,clr);
    clr = internal.stats.colorStringToRGB(clr);
end

if isempty(nbins)
    xbin = [];
    ybin =[];  
elseif isscalar(nbins)
    xbin = nbins;  % specified number of bins, same for x and y
    ybin = nbins;
else
    xbin = nbins(1); % specified number of bins, x and y different
    ybin = nbins(2);
end

% Set the bandwidth of kernel smoothing window
if isempty(ww)
    ww1 = [];
    ww2 = [];
elseif isnumeric(ww) && isequal(size(ww),[2,k]) && (all(~isnan(ww(:))))
    ww1 = ww(1,:);
    ww2 = ww(2,:);
else
    error(message('stats:scatterhist:BadBandwidth',k));
end

% Setup marginal distribution plotting function
if isempty(Xfunc) || isempty(Yfunc)
    if ~plotGroup || isempty(group)||k==1
        grp = ones(1,numel(x)); % grp for histogram and kernel
        if ~setFlag.style      
            style = 'bar';
        end
    else
        grp = group;
    end

    switch kernel
        case 'on'
            Xfunc = @()internal.stats.plotGroupedKSDensity(x,grp,...
                'Color',clr,'LineWidth',lw,'LineStyle',lStyle,'Width',ww1);
            Yfunc = @()internal.stats.plotGroupedKSDensity(y,grp,...
                'Color',clr,'LineWidth',lw,'LineStyle',lStyle,'Width',ww2);   
            
        case 'overlay'
            if plotGroup && k>1
                color = clr;
            else
                color = 'r';
            end
            Xfunc = @()xOverlay();
            Yfunc = @()yOverlay();
                
        case 'off'
                Xfunc = @()internal.stats.plotGroupedHist(x,grp,'NBins',xbin,...
                    'Color',clr,'DisplayStyle',style,'Norm','pdf','AxisOn',false,'PlotGroup',plotGroup);
                Yfunc = @()internal.stats.plotGroupedHist(y,grp,'NBins',ybin,...
                    'Color',clr,'DisplayStyle',style,'Norm','pdf','AxisOn',false,'PlotGroup',plotGroup);
    end
end
function xOverlay()
    internal.stats.plotGroupedHist(x,grp,'Color',clr,'NBins',xbin,...
        'DisplayStyle',style,'Norm','pdf','AxisOn',false,'PlotGroup',plotGroup);   
    hold on 
    internal.stats.plotGroupedKSDensity(x,grp,'Color',color,'LineStyle',lStyle,'LineWidth',lw,'Width',ww1);
    hold off
end
function yOverlay()
    internal.stats.plotGroupedHist(y,grp,'Color',clr,'NBins',ybin,...
        'DisplayStyle',style,'Norm','pdf','AxisOn',false,'PlotGroup',plotGroup);  
    hold on
    internal.stats.plotGroupedKSDensity(y,grp,'Color',color,'LineStyle',lStyle,'LineWidth',lw,'Width',ww2);
    hold off
end

% Set up positions for the plots
switch lower(direction)
    case 'in'
        inoutSign = 1;
    case 'out'
        inoutSign = -1;
    otherwise
        error(message('stats:scatterhist:BadDirection'));
end

switch lower(location)
    case {'ne' 'northeast'}
        scatterPosn = [.1 .1 .55 .55];
        scatterXAxisLoc = 'top'; scatterYAxisLoc = 'right';
        histXSign = -inoutSign;
        histYSign = -inoutSign;
    case {'se' 'southeast'}
        scatterPosn = [.1 .35 .55 .55];
        scatterXAxisLoc = 'bottom'; scatterYAxisLoc = 'right';
        histXSign = inoutSign;
        histYSign = -inoutSign;
    case {'sw' 'southwest'}
        scatterPosn = [.35 .35 .55 .55];
        scatterXAxisLoc = 'bottom'; scatterYAxisLoc = 'left';
        histXSign = inoutSign;
        histYSign = inoutSign;
    case {'nw' 'northwest'}
        scatterPosn = [.35 .1 .55 .55];
        scatterXAxisLoc = 'top'; scatterYAxisLoc = 'left';
        histXSign = -inoutSign;
        histYSign = inoutSign;
    otherwise
        error(message('stats:scatterhist:BadLocation'));
end

if isempty(parent)
    parent = ancestor(gca,{'figure','uipanel'});
end

% Clear the figure/uipanel for the new plot
if isequal(get(parent,'Type'),'figure')
    clf;
else
    hChildren = get(parent,'Children');
    if ~isempty(hChildren)
        beingDeleted = strcmpi('on', get(hChildren,'BeingDeleted')); 
        delete(hChildren(~beingDeleted));
    end
end
origVisible = get(parent,'Visible');
origVisibleMode = get(parent,'VisibleMode');
set(parent,'Visible','off');

hScatter = axes('parent',parent);
hYHist = axes('parent',parent);
hXHist = axes('parent',parent);

% Use the function handles Xfunc and Yfunc to plot the marginal
% distribution

% 'xfunc' and 'yfunc' are NOT external parameters. The function handles can
% reference functions of two kinds of signatures. The second signature may
% change in the future.
set(gcf,'NextPlot','Add');
set(gcf,'CurrentAxes',hXHist);
if nargin(Xfunc) ==0
    Xfunc();
else 
    Xfunc(x,group,'Color',clr,'LineStyle',lStyle);
end

set(gcf,'CurrentAxes',hYHist);
if nargin(Yfunc)==0
    Yfunc();
else
    Yfunc(y,group,'Color',clr,'LineStyle',lStyle);
end

view(hYHist,270,90); % Rotate the yhist plot
if histXSign<0
    set(hXHist,'YDir','reverse');    
end
if histYSign>0
    set(hYHist,'YDir','reverse');
end

% Scatter plot by gscatter 
set(gcf,'CurrentAxes',hScatter);
if isempty(inputname(1))
    xname = 'x';
else
    xname = inputname(1);
end
if isempty(inputname(2))
    yname = 'y';
else
    yname = inputname(2);
end
hlines = gscatter(x,y,group,clr,marker,siz,'on',xname,yname);

% Create Legend and set visibility accordingly
hLegend = findobj(parent,'Tag','legend');
hLegendToggle = findall(gcf,'Tag','Annotation.InsertLegend');
if ~isempty(hLegend)&& ~doLegend
    set(hLegend,'Visible','off');
    set(hLegendToggle,'State','off');
end
localSetLegend();

% Change legend toggle button call back functions
set(hLegendToggle,'ClickedCallback',@legendtogglecallback);

% adjust axes settings
set(hScatter,'Box','on');
axis(hScatter,'tight');
axis(hXHist,'off','tight');
axis(hYHist,'off','tight');

% Make sure that the axes do not overlay or overlap on the plots so that
% the edges of the xhist and yhist are not being cropped by the HG renderer
bestLimFunc = @(x,y)[min(x(1),y(1)), max(x(2),y(2))];
XhistXlim = bestLimFunc(get(hXHist,'Xlim'),get(hScatter,'Xlim'));
YhistXlim = bestLimFunc(get(hYHist,'Xlim'),get(hScatter,'Ylim'));
set(hScatter,'Xlim', XhistXlim + [-0.01, 0.01]*range(XhistXlim));
set(hScatter,'Ylim', YhistXlim + [-0.01, 0.01]*range(YhistXlim));

XhistYlim = get(hXHist,'YLim');
YhistYlim = get(hYHist,'YLim');
set(hXHist,'YLim',XhistYlim + [-0.01, 0.01]*range(XhistYlim)); 
set(hYHist,'YLim',YhistYlim + [-0.01, 0.01]*range(YhistYlim));

% Invisible texts to help calculate positions
txt1 = text(0,0,'','Visible','off','HandleVisibility','off');
txt2 = text(1,1,'','Visible','off','HandleVisibility','off');

set(hScatter,'Position',scatterPosn, 'XAxisLocation',scatterXAxisLoc, ...
             'YAxisLocation',scatterYAxisLoc, 'Tag','scatter');
set(hXHist,'Tag','xhist');
set(hYHist,'Tag','yhist');

% Sync axis
scatterhistXYLimCallback();

% Add listeners to resize or relimit histograms when the scatterplot changes
l1 = addlistener(hScatter,'MarkedClean',@scatterhistPositionCallback);
l2 = addlistener(hScatter,'MarkedClean',@scatterhistXYLimCallback);

addlistener(hlines(1),'ObjectBeingDestroyed',@gscatterDeleteCallback);

% Disable unrelated toolbar ui buttons
internal.stats.removeFigToolbarButton(gcf,{'Exploration.Brushing',...
                                           'Annotation.InsertColorbar',...
                                           'Exploration.Rotate',...
                                           'DataManager.Linking'});
% Disable pan on xhist and yhist plots
hpan = pan;
hpan.setAllowAxesPan([hXHist,hYHist],false);
% Disable zoom on xhist and yhist plots
hzoom = zoom;
hzoom.setAllowAxesZoom([hXHist,hYHist],false);
% Disable 3D rotation on xhist and yhist plots
hrotate = rotate3d;
hrotate.setAllowAxesRotate([hXHist,hYHist],false);

% Set current axis to the scatter axis
scatterhistPositionCallback();
set(gcf,'CurrentAxes',hScatter);
set(parent,'Visible',origVisible,'VisibleMode',origVisibleMode);

% Elect out of basic Fitting
set(hggetbehavior(hScatter,'datadescriptor'),'Enable',false)
set(hggetbehavior(hXHist,'datadescriptor'),'Enable',false)
set(hggetbehavior(hYHist,'datadescriptor'),'Enable',false)

% Make Axes Toolbar invisible for histogram axes
hXHist.Toolbar.Visible = 'off';
hYHist.Toolbar.Visible = 'off';

if nargout>0
    h  = [hScatter hXHist hYHist];
end

%--------------------------------------------------------------------------
% Nested subfunctions
    function gscatterDeleteCallback(~,~)
        % No cleanup is needed if the parent is deleted
        if strcmpi(get(parent,'BeingDeleted'),'on') || ...
            strcmpi(get(ancestor(parent,'figure'),'BeingDeleted'),'on')
                return;
        end
          
        % Reset the main scatter axes to default position
        if any(strcmpi('off',get(findobj(parent,'Type','axes'), 'BeingDeleted')))
            set(gca,'Position',get(0,'defaultaxesposition'));
        end
        
        % Delete the listeners
        delete(l1);
        delete(l2);
        % Delete the two hist axis
        if ishghandle(hXHist) && strcmpi(get(hXHist,'BeingDeleted'),'off')
            delete(hXHist);
        end
        if ishghandle(hYHist) && strcmpi(get(hYHist,'BeingDeleted'),'off')
            delete(hYHist);
        end
    end

    function scatterhistPositionCallback(~,~)
    % Position Callback Function
    
    posn = getrealposition(hScatter,txt1,txt2);

    % drawnow ensures that OuterPosition is updated
    drawnow('expose');

    oposn = get(hScatter,'OuterPosition');
    if ~ishandle(hLegend)
        posnLegend = [];
    else
        posnLegend = get(hLegend,'Position');
    end
    switch lower(location)
        case {'sw' 'southwest'}
            % vertically: margin, histogram, margin/4, scatterplot, margin
            vmargin = min(max(1 - oposn(2) - oposn(4), 0), oposn(2));
            posnHistX = [posn(1) vmargin posn(3) oposn(2)-1.25*vmargin];
            % horizontally: margin, histogram, margin/4, scatterplot, margin
            hmargin = min(max(1 - oposn(1) - oposn(3), 0), oposn(1));
            posnHistY = [hmargin posn(2) oposn(1)-1.25*hmargin posn(4)];
            % legend position
            if ~isempty(posnLegend)
                posnLegend = [max(oposn(1)-posnLegend(3),0.1*hmargin),...
                    max(oposn(2)-posnLegend(4),0.1*vmargin),...
                    posnLegend(3),posnLegend(4)];
            end
        case {'ne' 'northeast'}
            % vertically: margin, scatterplot, margin/4, histogram, margin
            vmargin = max(oposn(2), 0);
            posnHistX = [posn(1) oposn(2)+oposn(4)+.25*vmargin ...
                posn(3) 1-oposn(2)-oposn(4)-1.25*vmargin];
            % horizontally: margin, scatterplot, margin/4, histogram, margin
            hmargin = max(oposn(1), 0);
            posnHistY = [oposn(1)+oposn(3)+.25*hmargin    posn(2)...
                1-oposn(1)-oposn(3)-1.25*hmargin  posn(4)];
            % legend position
            if ~isempty(posnLegend)
                posnLegend = [...
                    min(oposn(1)+oposn(3),1-posnLegend(3)-0.1*hmargin),...
                    min(oposn(2)+oposn(4),1-posnLegend(4)-0.1*hmargin),...
                    posnLegend(3),...
                    posnLegend(4)];
            end
        case {'se' 'southeast'}
            % vertically: margin, histogram, margin/4, scatterplot, margin
            vmargin = max(1 - oposn(2) - oposn(4), 0);
            posnHistX = [posn(1) vmargin posn(3) oposn(2)-1.25*vmargin];
            % horizontally: margin, scatterplot, margin/4, histogram, margin
            hmargin = max(oposn(1), 0);
            posnHistY = [oposn(1)+oposn(3)+.25*hmargin     posn(2)...
                1-oposn(1)-oposn(3)-1.25*hmargin  posn(4)];
            % legend position
            if ~isempty(posnLegend)
                posnLegend = [...
                    min(oposn(1)+oposn(3), 1-posnLegend(3)-0.1*hmargin),...
                    max(oposn(2)-posnLegend(4),0.1*vmargin),...
                    posnLegend(3),...
                    posnLegend(4)];
            end
        case {'nw' 'northwest'}
            % vertically: margin, scatterplot, margin/4, histogram, margin
            vmargin = max(oposn(2), 0);
            posnHistX = [posn(1) oposn(2)+oposn(4)+.25*vmargin...
                posn(3) 1-oposn(2)-oposn(4)-1.25*vmargin];
            % horizontally: margin, histogram, margin/4, scatterplot, margin
            hmargin = max(1 - oposn(1) - oposn(3), 0);
            posnHistY = [hmargin posn(2) oposn(1)-1.25*hmargin posn(4)];
            % legend position
            if ~isempty(posnLegend)
                posnLegend = [max(oposn(1)-posnLegend(3),0.1*hmargin),...
                    min(oposn(2)+oposn(4),1-posnLegend(4)-0.1*hmargin),...
                    posnLegend(3),...
                    posnLegend(4)];
            end
    end
    posnHistX = max(posnHistX,[0 0 .05 .05]);
    posnHistY = max(posnHistY,[0 0 .05 .05]);

    if newPosition(posnHistX, get(hXHist, 'Position'))
        set(hXHist,'Position',posnHistX);
    end
    if newPosition(posnHistY, get(hYHist, 'Position'))
        set(hYHist,'Position',posnHistY);
    end
    if ~isempty(posnLegend) && newPosition(posnLegend, get(hLegend, 'Position'))...
            && strcmp(hLegend.PositionMode,'auto')
        % Only update the location if the user has not manually changed it
        set(hLegend, 'Position',posnLegend);
        set(hLegend, 'PositionMode', 'auto');
    end

    scatterhistXYLimCallback();

    end % scatterhistPositionCallback function

% -----------------------------
    function scatterhistXYLimCallback(~,~)
    %callback function to sync the x axes of xhist and yhist with the axes
    %of scatter plot.
    
        limx = get(hScatter,'XLim');
        limy = get(hScatter,'YLim');
        if ~isequal(limx,get(hXHist,'XLim'))
            set(hXHist,'XLim',limx);
        end
        if ~isequal(limy,get(hYHist,'XLim'))
            set(hYHist,'XLim',limy);
        end
        end


% ------------------------------
    function localSetLegend(~,~)
    if isempty(lStyle) % If no line object, no line added to legends
        return;
    elseif numel(findobj(hXHist,'Type','line','-depth',1)) <3
        % there should be more than 2 lines
        return;
    end

    hLegLines = findobj(hLegend,'Type','line');
    for i = 0:2:(numel(hLegLines)-2)
        set(hLegLines(end-i),'LineStyle',lStyle{i/2+1},...
            'LineWidth',lw(i/2+1));
    end
    end

% ------------------------------
    function legendtogglecallback(~,~)
    if ~ishandle(hLegend)         
        hLegend = findobj(gcf,'Type','legend');
    end
    if strcmp(get(hLegend,'Visible'),'off')
        set(hLegend,'Visible','on')
    else
        set(hLegend,'Visible','off')
    end
    end

end % SCATTERHIST main function


% SUBFUNCTIONS:
% -----------------------------
function p = getrealposition(a,txt1,txt2)
p = get(a,'Position');

% For non-warped axes (as in "axis square"), recalculate another way
if isequal(get(a,'WarpToFill'),'off')
    pctr = p([1 2]) + 0.5 * p([3 4]);
    xl = get(a,'XLim');
    yl = get(a,'YLim');
    
    % Use text to get coordinate (in points) of southwest corner
    set(txt1,'Units','data');
    set(txt1,'Position',[xl(1) yl(1)]);
    set(txt1,'Units','pixels');
    pSW = get(txt1,'Position');
    
    % Same for northeast corner
    set(txt2,'Units','data');
    set(txt2,'Position',[xl(2) yl(2)]);
    set(txt2,'Units','pixels');
    pNE = get(txt2,'Position');
    
    % Re-create position
    % Use min/max/abs in case one or more directions are reversed
    p = [min(pSW(1),pNE(1)), ...
         max(pSW(2),pNE(2)), ...
         abs(pNE(1)-pSW(1)), ...
         abs(pNE(2)-pSW(2))];
    p = hgconvertunits(ancestor(a,'figure'),p, ...
             'pixels','normalized',ancestor(a,'figure'));
    
    % Position to center
    p = [pctr(1)-p(3)/2, pctr(2)-p(4)/2, p(3), p(4)];
end
end % getrealposition function
    
function diff = newPosition(r, current)
% return true iff positions r and current differ by more than some small tolerance        
tol = 1e-5;
diff = max(abs(r - current)) > tol;
end
