function h = polytool(x,y,n,alpha,xname,yname)
%POLYTOOL Fits a polynomial to (x,y) data and displays an interactive graph.
%   POLYTOOL(X,Y,N,ALPHA) is a prediction plot that provides a Nth degree
%   polynomial curve fit to (x,y) data. It plots a 100(1 - ALPHA) percent
%   global confidence interval for predictions as two red curves.  The
%   default value for N is 1 (for a linear model) and the default value
%   for ALPHA is 0.05.  For any value X, the default curves define a
%   nonsimultaneous 95% confidence interval for a new observation of Y
%   taken at that X value.  (This differs from NLPREDCI, which computes
%   confidence intervals for the regression function.)
%
%   POLYTOOL(X,Y,N,ALPHA,XNAME,YNAME) labels the x and y values on the 
%   graphical interface using the strings XNAME and YNAME. Specify N and 
%   ALPHA as [] to use their default values. 
%
%   H = POLYTOOL(X,Y,N,ALPHA,XNAME,YNAME) outputs a vector of handles, H,
%   to the line objects (data, fit, lower bounds and upper bounds) in the plot.
%
%   You can drag the dashed black reference line and watch the predicted
%   values update simultaneously. Alternatively, you can get a specific
%   prediction by typing the "X" value into an editable text field.
%   A text field at the top allows you to change the degree of the 
%   polynomial fit. Use the push button labeled Export to move 
%   specified variables to the base workspace. Use the Bounds menu
%   to change the type of confidence bounds. Use the Method menu to
%   change the fitting method from least squares to robust.
   
%   Copyright 1993-2016 The MathWorks, Inc.


if nargin > 0
    x = convertStringsToChars(x);
end

if nargin > 4
    xname = convertStringsToChars(xname);
end

if nargin > 5
    yname = convertStringsToChars(yname);
end

action = 'start';
if (nargin > 0)
   if (ischar(x))
      action = x;
   end
end

%On initial call weed out NaNs
if (strcmp(action, 'start'))
   if (nargin == 0)
      noargmsg = getString(message('stats:polytool:dialogText_NoArgs'));
      x = (1:10)';
      y = 50 + 4*x - 0.75*x.^2 + randn(10,1);
      polytool(x,y);
      uiwait(warndlg(noargmsg, getString(message('stats:polytool:dialogTitle_PolynomialFitting')),'modal'));
      return
   end

   if (nargin < 2)
      error(message('stats:polytool:TooFewInputs'));
   end
   [anybad, wasnan, y, x] = statremovenan(y,x);
   if (anybad==2)
     error(message('stats:polytool:InputSizeMismatch'));
   end
   if (any(wasnan(:)))
      warning(message('stats:polytool:RemovedMissingValues'))
   end
end

%On recursive calls get all necessary handles and data.
if ~strcmp(action,'start')

   if nargin == 2,
      flag = y;
   end
   
   poly_fig = gcf;
   if ~strcmp(get(poly_fig,'Tag'),'polyfig')
      return
   end
   
   poly_axes = get(poly_fig,'CurrentAxes');
   ud = get(poly_fig,'UserData');
   if isempty(ud) || ~isstruct(ud)
       return
   end
   
   uiHandles    = ud.uicontrolHandles;
   x_field      = uiHandles.x_field;
   y_field      = uiHandles.y_field;
   degreetext   = uiHandles.degreetext;

   stats        = ud.stats;
   x            = ud.data(:,1);
   y            = ud.data(:,2);
   wasnan       = ud.wasnan;
   
   fitline        = ud.plotHandles.fitline;
   reference_line = ud.plotHandles.reference_line;
 
   degree         = str2double(get(degreetext,'String'));
   
   xrange = get(poly_axes,'XLim');
   newx = str2double(get(x_field,'String'));  
end

switch action

case 'start',

if nargin < 3, 
   degree = 1;
else
   if isempty(n)
      n = 1;
   end
   degree = n;
end

model = DegreeToModel(degree);

if nargin < 4 || isempty(alpha), 
    alpha = 0.05;
end
if isempty(alpha)
   alpha = 0.05;
end

if nargin >= 5
   xstr = xname;
   if isempty(xname)
      xstr = getString(message('stats:polytool:xlabel_Default'));
   end  
else
   xstr = getString(message('stats:polytool:xlabel_Default'));
end

if nargin == 6
   ystr = yname;
else
   ystr = getString(message('stats:polytool:ylabel_Default'));
end


% Create plot figure and axes.
poly_fig = figure('Visible','off','Tag','polyfig', ...
                  'NumberTitle','off', 'IntegerHandle','off', ...
                  'MenuBar','figure', ...
                  'Units','Normalized','ToolBar','none');
set(poly_fig,'Name',getString(message('stats:polytool:figureTitle',model)));
poly_axes = axes;

% This figure is not compatible with brushing
delete(uigettool(poly_fig,'Exploration.Brushing'))
delete(uigettool(poly_fig,'DataManager.Linking'))
delete(findall(poly_fig,'Tag','figDataManagerBrushTools'))
delete(findall(poly_fig,'Tag','figBrush'))
delete(findall(poly_fig,'Tag','figLinked'))

% Fit data.
nobs = max(size(y));
ncoeffs = 2:(nobs-1);
dflist = nobs - ncoeffs;
stats.crit = tinv(1 - alpha/2,dflist);
stats.fcrit = sqrt(ncoeffs .* finv(1-alpha, ncoeffs, dflist));
stats = dofit(x,y,degree,stats,0);

% Define the first of the UserData components.
ud.stats = stats;
ud.data  = [x(:) y(:)];
ud.wasnan = wasnan;
ud.texthandle = [];
ud.simflag = 0;        % nonsimultaneous confidence intervals
ud.obsflag = 1;        % for new observation, not for mean
ud.confflag = 1;       % show confidence bounds
ud.robflag = 0;        % not a robust fit
% Call local function to create plots.
[ud.plotHandles, newx] = makeplots(x,y,stats,degree,poly_axes,ud);
fitline = ud.plotHandles.fitline;
% Call local function to create uicontrols.
ud.uicontrolHandles = makeuicontrols(newx,xstr,ystr,stats,degree,ud);
setconf(ud);           % initialize settings on Bounds menu
setfit(ud);            % initialize settings on Method menu

set(poly_fig,'WindowButtonMotionFcn',@(varargin) motionFun(0,varargin{:}),...
             'WindowButtonDownFcn',@downFun,...
             'UserData',ud,'Visible','on','HandleVisibility','callback');
% Finished with polytool initialization.

case 'edittext',
    if isnan(newx)
      entry = get(x_field,'String');
      newx = get(x_field,'UserData');
      set(x_field,'String',num2str(double(newx)));
      msg = getString(message('stats:polytool:errorDialog_InvalidX',entry));
      localerrordlg(msg);
      return;
    end

    [minx, maxx] = plotlimits(xrange);
    if newx > maxx
        newx = maxx;
        set(x_field,'String',num2str(double(newx)));
    end
    if newx < minx
        newx = minx;
        set(x_field,'String',num2str(double(newx)));
    end

    [newy, deltay] = getyvalues(newx,stats,degree,ud);
	setyfields(newy,deltay,y_field,ud.confflag)
	setreferencelines(newx,newy,reference_line);

case 'change_degree',
    maxdf = length(stats.crit);
    if degree >= maxdf+1
       degree = get(degreetext,'UserData');
       set(degreetext,'String',num2str(double(degree)));
       localerrordlg(getString(message('stats:polytool:errorDialog_CannotFit')));
       return;
    end   
    if isempty(degree) || isnan(degree)
       degree = get(degreetext,'UserData');
       set(degreetext,'String',num2str(double(degree)));
       localerrordlg(getString(message('stats:polytool:errorDialog_NonPosDeg')));
       return;
    end
    if floor(degree) ~= degree || degree <= 0
       degree = get(degreetext,'UserData');
       set(degreetext,'String',num2str(double(degree)));
       localerrordlg(getString(message('stats:polytool:errorDialog_NonPosDeg')));
       return;
    end
    ud.stats = dofit(x,y,degree,ud.stats,ud.robflag);
    updateplot(poly_axes,x,y,ud,degree,fitline,reference_line,newx,y_field);

    set(degreetext,'UserData',degree);
    model = DegreeToModel(degree);
    set(poly_fig,'UserData',ud);
    set(poly_fig,'Name',getString(message('stats:polytool:figureTitle',model)));

case 'output',
    bmf = get(poly_fig,'WindowButtonMotionFcn');
    bdf = get(poly_fig,'WindowButtonDownFcn');
    set(poly_fig,'WindowButtonMotionFcn','','WindowButtonDownFcn','');

    % Strings used in export dialog
    labels = {getString(message('stats:polytool:dialogText_Parameters')), ...
        getString(message('stats:polytool:dialogText_ParameterCI')), ...
        getString(message('stats:polytool:dialogText_Prediction')), ...
        getString(message('stats:polytool:dialogText_PredictionCI')), ...
        getString(message('stats:polytool:dialogText_Residuals'))};
    varnames = {'beta', 'betaci', 'yhat', 'yci', 'residuals'};
    dialogTitle = getString(message('stats:polytool:dialogTitle_Export'));

	structure = stats.structure;
	beta = stats.beta;
    R = structure.R;
    df = structure.df;
    rmse = structure.normr / sqrt(df);
    Rinv = eye(size(R))/R;
    dxtxi = sqrt(sum(Rinv'.*Rinv'));
    dbeta = dxtxi*stats.crit(degree)*rmse;
    yhat = polyval(beta,x);
	[newy, deltay] = getyvalues(newx,stats,degree,ud);

    fullresid = y-yhat;
    if (any(wasnan))
       fullresid = statinsertnan(wasnan,fullresid);
    end
    
    items = {beta, [beta-dbeta; beta+dbeta], newy, [newy-deltay newy+deltay], fullresid};
    export2wsdlg(labels, varnames, items, dialogTitle);
    set(poly_fig,'WindowButtonMotionFcn', bmf);
    set(poly_fig,'WindowButtonDownFcn', bdf);

case 'conf',
  if (nargin > 1), conf = flag; end
  switch(conf)
   case 1, ud.simflag = 1;             % simultaneous
   case 2, ud.simflag = 0;             % nonsimultaneous
   case 3, ud.obsflag = 0;             % for the mean (fitted line)
   case 4, ud.obsflag = 1;             % for a new observation
   case 5, ud.confflag = ~ud.confflag; % no confidence intervals
  end
  setconf(ud);
  set(poly_fig, 'UserData', ud);

  % Update bounds
    updateplot(poly_axes,x,y,ud,degree,fitline,reference_line,newx,y_field);

case 'method',
   if (nargin > 1), fittype = flag; end
   switch(fittype)
    case 1, ud.robflag = 0;             % least squares (not robust)
    case 2, ud.robflag = 1;             % robust fit
   end
   setfit(ud);

   % Update
   ud.stats = dofit(x,y,degree,ud.stats,ud.robflag);
   set(poly_fig, 'UserData', ud);
   updateplot(poly_axes,x,y,ud,degree,fitline,reference_line,newx,y_field);
end

if nargout == 1
   h = fitline;
end

function stats = dofit(x,y,degree,stats,robflag)
%DOFIT    Do fitting, either least squares or robust
if (~robflag)
   [stats.beta, stats.structure] = polyfit(x,y,degree);
else
   % Construct Vandermonde matrix
   V = ones(length(x),degree+1);
   if (degree>0), V(:,end-1) = x; end
   for j = degree-1:-1:1
      V(:,j) = x.*V(:,j+1);
   end
   [beta, other] = robustfit(V,y,'bisquare',[],0);
   
   % Compute robust version of stuff polyfit places into structure
   S.R = other.R;
   S.df = other.dfe;
   S.normr = sqrt(other.dfe) * other.s;
   stats.beta = beta';
   stats.structure = S;
end


function setyfields(newy,deltay,y_field,confflag)
% SETYFIELDS Local function for changing the Y field values.
set(y_field(1),'String',num2str(double(newy)));
if (isnan(deltay) || ~confflag)
   set(y_field(2),'String','');
   set(y_field(3),'String','');
else
   set(y_field(2),'String',num2str(double(deltay)));
   set(y_field(3),'String',' +/-');
end

function setreferencelines(newx,newy,reference_line)
% SETREFERENCELINES Local function for moving reference lines on polytool figure.
xrange = get(gca,'Xlim');
yrange = get(gca,'Ylim');
set(reference_line(1),'YData',yrange,'Xdata',[newx newx]);   % vertical reference line.
set(reference_line(2),'XData',xrange,'Ydata',[newy newy]);   % horizontal reference line.
 

function [yhat, deltay] = getyvalues(x,stats,degree,ud)
% GETYVALUES Local function for generating Y predictions and confidence intervals.
S = stats.structure;
[yhat, deltay]=polyval(stats.beta,x,S);
df = S.df;
normr = S.normr;
if (ud.obsflag || S.normr==0)
   dy = deltay;
else
   e = (deltay * sqrt(df) / normr).^2;
   e = sqrt(max(0, e-1));
   dy = normr/sqrt(df)*e;
end
if (ud.simflag)
   crit = stats.fcrit(degree);
else
   crit = stats.crit(degree);
end
deltay = dy * crit;


function [minx, maxx] = plotlimits(x)
% PLOTLIMITS   Local function to control the X-axis limit values.
maxx = max(x(:));
minx = min(x(:));
xrange = maxx - minx;
maxx = maxx + 0.025 * xrange;
minx = minx - 0.025 * xrange;

function [plotHandles, newx] = makeplots(x,y,stats,degree,poly_axes,ud)
% MAKEPLOTS   Local function to create the plots in polytool.
[minx, maxx] = plotlimits(x);
xfit = linspace(minx,maxx,41)';
[yfit, deltay] = getyvalues(xfit,stats,degree,ud);

set(poly_axes,'UserData','poly_axes','XLim',[minx maxx],'Box','on');
set(poly_axes,'NextPlot','add','Position',[.21 .18 .75 .72]);

% Plot prediction function with uncertainty bounds.
plotHandles.fitline = plot(x,y,'b+',xfit,yfit,'g-',xfit,yfit-deltay,'r--',...
                           xfit,yfit+deltay,'r--');
% set tags
fitline = plotHandles.fitline;
set(fitline(1),'Tag','data');
set(fitline(2),'Tag','fit');
set(fitline(3),'Tag','lowerBound');
set(fitline(4),'Tag','upperBound');

yrange = get(poly_axes,'YLim');
xrange = get(poly_axes,'XLim');

% Plot Reference Lines
newx = xfit(21);
newy = yfit(21);
reference_line = zeros(1,2);
reference_line(1) = plot([newx newx],yrange,'k--','Tag','VerticalReferenceLine');
reference_line(2) = plot(xrange,[newy newy],'k:','Tag','HorizontalReferenceLine');
set(reference_line(1),'ButtonDownFcn','polytool(''down'')');
plotHandles.reference_line = reference_line;

function uicontrolHandles = makeuicontrols(newx,xstr,ystr,stats,degree,ud)
% MAKEUICONTROLS   Local function to create uicontrols for polytool. 

if strcmp(computer,'MAC2')
  offset = -0.01;
else
  offset = 0;
end

fcolor = get(gcf,'Color');

xfieldp = [.50 .07 .15 .05];
yfieldp = [.01 .45 .13 .04]; 

[newy, deltay] = getyvalues(newx,stats,degree,ud);

uicontrolHandles.x_field = uicontrol('Style','edit','Units','normalized','Position',xfieldp,...
    'String',num2str(double(newx)),...
    'BackgroundColor','white','UserData',newx,...
    'CallBack','polytool(''edittext'')','String',num2str(double(newx)),...
    'Tag','XValues');

uicontrol('Style','text','Units','normalized',...
    'Position',yfieldp + [0 0.20 0 0],...
    'String',ystr,...
    'BackgroundColor',fcolor,'ForegroundColor','k','Tag','YValuesLabel');

uicontrolHandles.y_field(1)= uicontrol('Style','text','Units','normalized',...
    'Position',yfieldp + [0 .14 0 0],...
    'String',num2str(double(newy)),...
    'BackgroundColor',fcolor,'ForegroundColor','k','Tag','YValues');

uicontrolHandles.y_field(3)= uicontrol('Style','text','Units','normalized',...
    'Position',yfieldp + [0 .07 -0.01 0], ...
    'String',getString(message('stats:polytool:uicontrol_PlusMinus')),...
    'BackgroundColor',fcolor,'ForegroundColor','k');

uicontrolHandles.y_field(2)= uicontrol('Style','text','Units','normalized',...
    'Position',yfieldp,...
    'String',num2str(double(deltay)),...
    'BackgroundColor',fcolor,'ForegroundColor','k','Tag','YValueBounds');

uicontrol('Style','text','Units','normalized',...
    'Position',xfieldp - [0 0.07 0 0],'ForegroundColor','k',...
    'BackgroundColor',fcolor,...
    'String',xstr,...
    'Tag','XValuesLabel');

uicontrol('Style','Text',...
    'String',getString(message('stats:polytool:uicontrol_Degree')),...
    'Units','normalized',...
    'Position',xfieldp + [-0.03 0.87+offset -0.05 0],'ForegroundColor','k',...
    'BackgroundColor',fcolor);

uicontrolHandles.degreetext=uicontrol('Style','Edit',...
    'String',num2str(double(degree)),...
    'Units','normalized','UserData',degree,'BackgroundColor','white',...
    'Position',xfieldp + [0.07 0.87 -0.10 0],...
    'CallBack','polytool(''change_degree'')',...
    'Tag','Degree');

uicontrolHandles.outputpopup=uicontrol(...
    'String', getString(message('stats:polytool:uicontrol_Export')),...
    'Units','normalized','Position',[0.01 0.07 0.17 0.05],...
    'CallBack','polytool(''output'')','Tag','ExportButton');


uicontrol('Style','Pushbutton','Units','normalized',...
    'Position',[0.01 0.01 0.17 0.05],'Callback','close',...
    'String',getString(message('stats:polytool:uicontrol_Close')),...
    'Tag','CloseButton');

% Create menu for controlling confidence bounds
f = uimenu(...
    'Label',getString(message('stats:polytool:menuTitle_Bounds')),...
    'Position', 4, 'UserData','conf', 'Tag','Bounds');
uimenu(f,...
    'Label',getString(message('stats:polytool:menuLabel_Simultaneous')),...
    'Tag','Simultaneous','Callback','polytool(''conf'',1)', 'UserData',1);

uimenu(f,...
    'Label',getString(message('stats:polytool:menuLabel_Nonsimultaneous')),...
    'Tag','Nonsimultaneous','Callback','polytool(''conf'',2)', 'UserData',2);
uimenu(f,...
    'Label',getString(message('stats:polytool:menuLabel_Curve')),...
    'Tag','Curve', 'Separator','on',...
    'Callback','polytool(''conf'',3)', 'UserData',3);
uimenu(f,...
    'Label',getString(message('stats:polytool:menuLabel_Observation')),...
    'Tag','Observation', ...
    'Callback','polytool(''conf'',4)', 'UserData',4);
uimenu(f,...
    'Label',getString(message('stats:polytool:menuLabel_None')),...
    'Tag','None', 'Separator','on', ...
    'Callback','polytool(''conf'',5)', 'UserData',5);

% Create menu for controlling fit type
f = uimenu(...
    'Label',getString(message('stats:polytool:menuTitle_Method')),...
    'Position', 5, 'UserData','method','Tag','Method');
uimenu(f,...
    'Label',getString(message('stats:polytool:menuLabel_LeastSquares')),...
    'Tag','LeastSquares', ...
    'Callback','polytool(''method'',1)', 'UserData',1);
uimenu(f,...
    'Label',getString(message('stats:polytool:menuLabel_Robust')),...
    'Tag','Robust', ...
    'Callback','polytool(''method'',2)', 'UserData',2);

% -----------------------
function setconf(ud)
%SETCONF Update menus to reflect confidence bound settings

ma = findobj(gcf, 'Type','uimenu', 'UserData','conf');
hh = findobj(ma, 'Type','uimenu');

set(hh, 'Checked', 'off');          % uncheck all menu items

hh = findobj(ma, 'Type', 'uimenu', 'UserData', 2-ud.simflag);
if (length(hh) == 1), set(hh, 'Checked', 'on'); end

hh = findobj(ma, 'Type', 'uimenu', 'UserData', 3+ud.obsflag);
if (length(hh) == 1), set(hh, 'Checked', 'on'); end

hh = findobj(ma, 'Type', 'uimenu', 'UserData', 5);
if (~ud.confflag), set(hh, 'Checked', 'on'); end

% -----------------------
function setfit(ud)
%SETFIT Update menus to reflect fit type

ma = findobj(gcf, 'Type','uimenu', 'UserData','method');
hh = findobj(ma, 'Type','uimenu');

set(hh, 'Checked', 'off');          % uncheck all menu items

hh = findobj(ma, 'Type', 'uimenu', 'UserData', 1+ud.robflag);
if (length(hh) == 1), set(hh, 'Checked', 'on'); end

% -----------------------------------
function updateplot(poly_axes,x,y,ud,degree,fitline,reference_line,newx,y_field)
%UPDATEPLOT Update plot after degree or confidence setting change

[minx, maxx] = plotlimits(x);
xfit = linspace(minx,maxx,41);
[yfit, deltay] = getyvalues(xfit,ud.stats,degree,ud);
if (ud.confflag)
   ylo = yfit - deltay;
   yhi = yfit + deltay;
   [miny, maxy] = plotlimits([yhi;ylo]);
else
   ylo = NaN(size(yfit));
   yhi = ylo;
   [miny, maxy] = plotlimits(yfit);
end
ymax = max(y);
ymin = min(y);
delta = .05 * (ymax - ymin);
miny = min(miny, ymin - delta);
maxy = max(maxy, ymax + delta);

set(fitline(1),'XData',x,'YData',y);
set(fitline(2),'XData',xfit,'YData',yfit);
set(fitline(3),'XData',xfit,'YData',ylo);
set(fitline(4),'XData',xfit,'YData',yhi);

set(poly_axes,'YLim',[miny maxy]);
[newy, deltay] = getyvalues(newx,ud.stats,degree,ud);
setreferencelines(newx,newy,reference_line);
if (~ud.confflag), deltay = NaN; end	
setyfields(newy,deltay,y_field,ud.confflag)

% -----------------------------------
function localerrordlg(msg)
uiwait(warndlg(msg,getString(message('stats:polytool:dialogTitle_PolynomialFitting')),'modal'))

% -----------------------------------
function downFun(varargin)
poly_fig = gcbf;
p = get(poly_fig,'CurrentPoint');
if p(1) < .21 || p(1) > .96 || p(2) > 0.86 || p(2) < 0.18
   return
end

% Find mouse click location, update text fields and ref line
set(poly_fig,'Pointer','crosshair');
poly_axes = get(poly_fig, 'CurrentAxes');
cp = get(poly_axes,'CurrentPoint');
[minx, maxx] = plotlimits(get(poly_axes,'Xlim'));
newx = max(minx, min(maxx, cp(1,1)));

ud = get(poly_fig,'UserData');

uiHandles    = ud.uicontrolHandles;
x_field      = uiHandles.x_field;
y_field      = uiHandles.y_field;
degreetext   = uiHandles.degreetext;
degree       = str2double(get(degreetext,'String'));

[newy, deltay] = getyvalues(newx,ud.stats,degree,ud);
set(x_field,'String',num2str(double(newx)));
setyfields(newy,deltay,y_field,ud.confflag)

setreferencelines(newx,newy,ud.plotHandles.reference_line);

% Update motion function to indicate button was down
set(gcf,'WindowButtonMotionFcn',@(varargin) motionFun(1,varargin{:}),...
        'WindowButtonUpFcn',@upFun);

% ----
function motionFun(flag,varargin)
poly_fig = gcbf;
ud = get(poly_fig,'UserData');
poly_axes = get(poly_fig,'CurrentAxes');
xrange = get(poly_axes,'XLim');
yrange = get(poly_axes,'YLim');
uiHandles    = ud.uicontrolHandles;
x_field      = uiHandles.x_field;
y_field      = uiHandles.y_field;
newx = str2double(get(x_field,'String'));  
degreetext   = uiHandles.degreetext;
degree       = str2double(get(degreetext,'String'));

[minx, maxx] = plotlimits(xrange);
 
if flag==0
    % Button was up, set cursor if over ref line
    cursorstate = get(poly_fig,'Pointer');
    cp = get(poly_axes,'CurrentPoint');
    cx = cp(1,1);
    cy = cp(1,2);
    fuzz = 0.01 * (maxx - minx);
    online = cy > yrange(1) & cy < yrange(2) & cx > newx - fuzz & cx < newx + fuzz;
    if online && strcmp(cursorstate,'arrow'),
        set(poly_fig,'Pointer','crosshair');
    elseif ~online && strcmp(cursorstate,'crosshair'),
        set(poly_fig,'Pointer','arrow');
    end
else
    % Button was down, drag ref line to new location
    cp = get(poly_axes,'CurrentPoint');
    if ~isinaxes(cp, poly_axes)
        if flag == 1
            set(poly_fig,'Pointer','arrow');
            set(poly_fig,'WindowButtonMotionFcn',@(varargin) motionFun(2,varargin{:}));
        end
        return;
    elseif flag == 2
        set(poly_fig,'Pointer','crosshair');
        set(poly_fig,'WindowButtonMotionFcn',@(varargin) motionFun(1,varargin{:}));
    end
    newx = max(minx, min(maxx, cp(1,1)));

    [newy, deltay] = getyvalues(newx,ud.stats,degree,ud);
    set(x_field,'String',num2str(double(newx)));
    setyfields(newy,deltay,y_field,ud.confflag);

    reference_line = ud.plotHandles.reference_line;
    setreferencelines(newx,newy,reference_line);
end

% ----
function upFun(varargin)
% Reset motion function to indicate button is now up
set(gcbf,'WindowButtonMotionFcn',@(varargin) motionFun(0,varargin{:}),...
             'WindowButtonUpFcn','');

% ----
function model = DegreeToModel(degree)
% convert numeric degree to model string

switch degree
    case 1
        model = getString(message('stats:polytool:text_Linear'));
    case 2
        model = getString(message('stats:polytool:text_Quadratic'));
    case 3
        model = getString(message('stats:polytool:text_Cubic'));
    case 4
        model = getString(message('stats:polytool:text_Quartic'));
    case  5
        model = getString(message('stats:polytool:text_Quintic'));
    otherwise
        model = getString(message('stats:polytool:text_NthOrder', ...
            num2str(double(degree))));
end
