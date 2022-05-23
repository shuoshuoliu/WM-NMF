function nlintool(x,y,model,beta0,varargin)
%NLINTOOL Interactive graphical tool for nonlinear fitting and prediction.
%   NLINTOOL(X,Y,FUN,BETA0) is a prediction plot that provides a 
%   nonlinear curve fit to (x,y) data. It plots 95% global confidence 
%   bands for predictions as two red curves. 
%
%   Y is a vector.  X is a vector or matrix with the same number of rows
%   as Y.  FUN is a function that accepts two arguments, a coefficient
%   vector and an array of X values, and returns a vector of fitted Y
%   values.  BETA0 is a vector containing initial guesses for the
%   coefficients.
%
%   NLINTOOL(...,'PARAM1',val1,'PARAM2'val2,...) specifies one or more
%   of the following parameter name/value pairs:
%
%      'alpha'    An ALPHA value for 100(1-ALPHA) percent confidence bands
%      'xname'    A string array, cell array of strings, or a character
%                 array, containing X variable names
%      'yname'    Y variable name
%      'plotdata' Either 'on' to show data on plot, or 'off' (default) for
%                 fit only.  'on' is not allowed with multiple X variables.
%
%   You can drag the dotted white reference line and watch the predicted
%   values update simultaneously.  Alternatively, you can get a specific
%   prediction by typing the "X" value into an editable text field.
%   Use the Export button to move specified variables to the base workspace.
%   Use the Bounds menu to change the type of confidence bounds.
%
%   Examples
%   --------
%   FUN can be specified using @:
%      nlintool(x, y, @myfun, b0)
%   where MYFUN is a MATLAB function such as:
%      function yhat = myfun(beta, x)
%      b1 = beta(1);
%      b2 = beta(2);
%      yhat = 1 ./ (1 + exp(b1 + b2*x));
%
%   FUN can also be an anonymous function:
%      fun = @(b,x) 1 ./ (1 + exp(b(1) + b(2)*x));
%      nlintool(x, y, fun, b0)
%
%   See also NLINFIT, NLPARCI, NLPREDCI, NLMEFIT.
   
%   Copyright 1993-2016 The MathWorks, Inc.


if nargin > 2
    model = convertStringsToChars(model);
end

if nargin > 4
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if (nargin == 0)
    action = 'start';
elseif ~ischar(x) 
    action = 'start';
else
    action = x;
end

%On recursive calls get all necessary handles and data.
if ~strcmp(action,'start')
   if nargin == 2,
      flag = y;
   end
   
   nlin_fig = gcbf;
   if isempty(nlin_fig) || ~isequal(get(nlin_fig,'Tag'),'nlinfig')
       return
   end
   ud = get(nlin_fig,'Userdata');
   
   R             = ud.R;
   beta          = ud.beta;
   model         = ud.model;
   x_field       = ud.x_field;
   nlin_axes     = ud.nlin_axes;
   xsettings     = ud.xsettings;
   y_field1      = ud.y_field(1);
   y_field2      = ud.y_field(2);
   last_axes     = ud.last_axes;  
   n             = size(x_field,2);

   xrange         = zeros(n,2);

   for k = 1:n         
       xrange(k,1:2) = get(nlin_axes(k),'XLim');
   end
   newy = str2double(get(y_field1,'String'));  
   newdeltay = str2double(get(y_field2,'String'));  
end

switch action

case 'start'

if (nargin == 0)
   msg = sprintf([getString(message('stats:nlintool:UseWithYourOwnDataType')) ...
                  '   nlintool(x,y,fun,beta0)\n\n' ...
                  getString(message('stats:nlintool:TypehelpClickOK'))]);
   s = load('reaction');
   x = s.reactants;
   y = s.rate;
   b = s.beta;
   nlintool(x, y, 'hougen', b);
   warndlg(msg, getString(message('stats:nlintool:NonlinearFitting')),'modal');
   return
end

if (nargin<4)
    error(message('stats:nlintool:TooFewInputs'));
end

if (size(x,1)==1), x = x(:); end
if (size(y,1)==1), y = y(:); end
wasnan = (isnan(y) | any(isnan(x),2));
if (any(wasnan))
   x(wasnan,:) = [];
   y(wasnan)=[];
end

% Defaults for optional inputs
alpha = 0.05;
xname = '';
yname = '';
if size(x,2)>1
   plotdata = 'off';
else
   plotdata = 'on';
end

% Deal with grandfathered syntax:  nlintool(x,y,fun,b0,alpha,xname,yname)
if ~isempty(varargin) && isnumeric(varargin{1})
   if length(varargin)>3
       error(message('stats:nlintool:BadSyntax'));
   end
   if ~isempty(varargin{1})
       alpha = varargin{1};
   end
   if length(varargin)>=2
       xname = varargin{2};
   end
   if length(varargin)>=3
       yname = varargin{3};
   end
else
   okargs =   {'alpha' 'xname' 'yname' 'plotdata'};
   defaults = {alpha   xname   yname   plotdata};
   [alpha xname yname plotdata] = ...
                   internal.stats.parseArgs(okargs,defaults,varargin{:});
end

if ~isequal(plotdata,'off')
   if ~isequal(plotdata,'on')
       error(message('stats:nlintool:BadPlotDataNotOnOff'));
   elseif size(x,2)>1
       warning(message('stats:nlintool:BadPlotDataMultipleX'));
       plotdata = 'off';
   end
end

if isempty(xname)
    xname = strcat({'X'}, strtrim(cellstr(int2str((1:size(x,2))'))));
else
    xname = cellstr(xname); % char matrix -> cellstr
end
if isempty(yname)
    yname = getString(message('stats:nlintool:PredictedY'));
end

% Fit nonlinear model.  Catch any error so we can restore the warning
% state before re-throwing the error.
[lastwmsg,lastwid] = lastwarn;
lastwstate = warning;
lastwarn('');
warning('off'); %#ok<WNOFF>

try
   [ud.beta, residuals, J] = nlinfit(x,y,model,beta0);
catch ME
   warning(lastwstate);
   rethrow(ME);
end

wmsg = lastwarn;
lastwarn(lastwmsg,lastwid);
warning(lastwstate);

p = numel(ud.beta);
df = max(size(y)) - p;

ud.residuals = residuals;
ud.rmse = sqrt(sum(residuals.*residuals)./df);
ud.crit_f = sqrt(p*finv(1 - alpha,p,df))*ud.rmse;

% Get crit_fobs for predictions at new observations (not fitted curve).
try
    %----------------------------------------------------------------------
    %|  When using a non-constant error model use the 'WeightedJacobian'  |
    %|  option and pass in J returned by nlinfit. Also,                   |                                    |
    %|  (1) For fixed weights, use 'ErrorVarianceFit' = 1./weights and    |
    %|      'ErrorVariancePred' as ones(length(weights),1).               |
    %|  (2) Otherwise, use defaults for these parameters.                 |
    %----------------------------------------------------------------------
    % Use J returned by nlinfit.
    sch = internal.stats.getscheffeparam('UnweightedJacobian',J);
catch ignore %#ok<NASGU>
    % Use a conservative setting.
    sch = (numel(ud.beta)+1);
end
ud.crit_fobs = sqrt(sch*finv(1 - alpha,sch,df))*ud.rmse;

ud.crit_t = tinv(1 - alpha/2,df)*ud.rmse;
ud.crit = ud.crit_f;
ud.model = model;   
ud.simflag = 1;        % simultaneous confidence intervals
ud.obsflag = 0;        % not for new observation (for mean)
ud.confflag = 1;       % show confidence bounds

[~,R] = qr(J,0);
ud.R = R;

% Set positions of graphic objects
maxx = max(x);
minx = min(x);
xrange = maxx - minx;
xlims = [minx - 0.025 .* xrange; maxx + 0.025 .* xrange]';

n = size(x,2);

nlin_axes      = zeros(n,1);
fitline        = zeros(3,n);
reference_line = zeros(n,2);

xfit      = double(xrange(ones(41,1),:)./40);
xfit(1,:) = minx;
xfit      = cumsum(xfit);

avgx      = mean(x);
xsettings = avgx(ones(41,1),:);
if (isa(model, 'inline'))
   modname = getString(message('stats:nlintool:NonlinearFitOf',formula(model)));
elseif (isa(model, 'function_handle'))
   modname = getString(message('stats:nlintool:NonlinearFit'));
else
   modname = getString(message('stats:nlintool:NonlinearFitModel',model));
end
nlin_fig = figure('Units','Normalized','Interruptible','on','Position',[0.05 0.35 0.90 0.5],...
             'NumberTitle','off', 'IntegerHandle','off', ...
             'Name',modname,'Tag','nlinfig',...
             'MenuBar','figure', ...
             'ToolBar','none');
yextremes = zeros(n, 2);
addplotdata = isequal(plotdata,'on');

% This figure is not compatible with brushing
delete(uigettool(nlin_fig,'Exploration.Brushing'))
delete(uigettool(nlin_fig,'DataManager.Linking'))
delete(findall(nlin_fig,'Tag','figDataManagerBrushTools'))
delete(findall(nlin_fig,'Tag','figBrush'))
delete(findall(nlin_fig,'Tag','figLinked'))

for k = 1:n
   % Create an axis for each input (x) variable
    axisp   = [.18 + (k-1)*.80/n .22 .80/n .68];

    nlin_axes(k) = axes;
    set(nlin_axes(k),'XLim',xlims(k,:),'Box','on','NextPlot','add',...
        'Position',axisp,'GridLineStyle','none');
    if k>1
       set(nlin_axes(k),'Yticklabel',[]);
    end

    % Show data if requested, which implies just one column in x
    if addplotdata
       line(x(:,k),y, 'Linestyle','none', 'Marker','o', 'Parent',nlin_axes(k));
    end

   % Calculate y values for fitted line plot.
   ith_x      = xsettings;
   ith_x(:,k) = xfit(:,k);
   [yfit, deltay] = PredictionsPlusError(ith_x,ud);
 
   % Plot prediction line with confidence intervals
   set(nlin_fig,'CurrentAxes',nlin_axes(k));
   fitline(1:3,k) = plot(ith_x(:,k),yfit,'g-',ith_x(:,k),yfit-deltay,'r--',ith_x(:,k),yfit+deltay,'r--');

   % Calculate data for vertical reference lines
   yextremes(k,:) = get(gca,'YLim');
end  % End of the plotting loop over all the axes.
statbrushlink(fitline,false);

ymin = min(yextremes(:,1));
ymax = max(yextremes(:,2));

[newy, newdeltay] = PredictionsPlusError(avgx,ud);

for k = 1:n
   set(nlin_axes(k),'YLim',[ymin ymax]);
   xlimits = get(nlin_axes(k),'XLim');
   %   Create Reference Lines
   xvert = repmat(avgx(k), size(yextremes(k,:)));
   set(nlin_fig,'CurrentAxes',nlin_axes(k));
   reference_line(k,1) = plot(xvert,[ymin ymax],'--');
   reference_line(k,2) = plot(xlimits,[newy newy],':');
end
statbrushlink(reference_line,false);

% Reference lines should not cause the axis limits to extend
set(reference_line(:),'XLimInclude','off','YLimInclude','off');

uihandles = MakeUIcontrols(xname,nlin_fig,newy,newdeltay,yname,avgx);
ud.x_field = uihandles.x_field; 
ud.y_field  = uihandles.y_field;

ud.texthandle = [];

ud.xfit = xfit;
ud.nlin_axes = nlin_axes;
ud.xsettings = xsettings;
ud.fitline = fitline;
ud.reference_line = reference_line;
ud.last_axes = 0;
ud.wasnan = wasnan;
setconf(ud);           % initialize settings on Bounds menu
set(nlin_fig,'UserData',ud,'HandleVisibility','callback', ...
             'BusyAction','queue', ...
             'WindowButtonDownFcn',@downFun,...
             'WindowButtonUpFcn',@upFun, ...
             'WindowButtonMotionFcn',@(varargin) motionFun(0,varargin{:}));

if ~isempty(wmsg)
   warndlg(wmsg, getString(message('stats:nlintool:NonlinearFitting')),'modal');
end

% Finished with plot startup function.

case 'edittext',
   cx    = str2double(get(x_field(flag),'String'));
   
   if isnan(cx)
      set(x_field(flag),'String',num2str(double(xsettings(1,flag))));
      % Create Bad Settings Warning Dialog.
      warndlg(getString(message('stats:nlintool:TypeOnlyNumbers')), getString(message('stats:nlintool:NonlinearFitting')),'modal');
      return
   end  
   xl = get(nlin_axes(flag),'Xlim');
   if cx < xl(1) || cx > xl(2)
      % Create Bad Settings Warning Dialog.
      warndlg(getString(message('stats:nlintool:NumberOutsideRange')), getString(message('stats:nlintool:NonlinearFitting')),'modal');
      set(x_field(flag),'String',num2str(double(xsettings(1,flag))));
      return
   end
   
   last_axes = flag;
   ud.last_axes = last_axes;
   ud = updateplot(nlin_fig, ud);
   last_axes = 0;
   ud.last_axes = last_axes;
   set(nlin_fig,'Userdata',ud);       

case 'output',
    bmf = get(nlin_fig,'WindowButtonMotionFcn');
    bdf = get(nlin_fig,'WindowButtonDownFcn');
    set(nlin_fig,'WindowButtonMotionFcn','');
    set(nlin_fig,'WindowButtonDownFcn','');

    labels = { getString(message('stats:nlintool:Parameters')),...
               getString(message('stats:nlintool:ParameterCI')),...
               getString(message('stats:nlintool:Prediction')),...
               getString(message('stats:nlintool:PredictionCI')),...
               getString(message('stats:nlintool:RMSE')),...
               getString(message('stats:nlintool:Residuals'))};
              
    varnames = {'beta', 'betaci', 'ypred', 'ypredci', 'rmse', 'residuals'};

    fullresids = NaN(size(ud.wasnan));
    fullresids(~ud.wasnan) = ud.residuals;

    % code copied from nlparci
    Rinv = R\eye(size(R));
    diag_info = sum(Rinv.*Rinv,2);
    deltab = sqrt(diag_info) .* ud.crit_t;
    parci = [(beta(:) - deltab) (beta(:) + deltab)];
    yci = [newy-newdeltay newy+newdeltay];
        
    items = {beta, parci, newy, yci, ud.rmse, fullresids};  
    export2wsdlg(labels, varnames, items, ...
         getString(message('stats:nlintool:ExportToWorkspace')));
    set(nlin_fig,'WindowButtonMotionFcn',bmf);
    set(nlin_fig,'WindowButtonDownFcn',bdf);

case 'conf',
   if (nargin > 1), conf = flag; end
   switch(conf)
    case 1, ud.simflag = 1;             % simultaneous
    case 2, ud.simflag = 0;             % non-simultaneous
    case 3, ud.obsflag = 0;             % for the mean (fitted line)
    case 4, ud.obsflag = 1;             % for a new observation
    case 5, ud.confflag = ~ud.confflag; % no confidence intervals
   end

   if (ud.simflag)
       if (ud.obsflag)
           % Simultaneous confidence interval for new observations.
           ud.crit = ud.crit_fobs;
       else
           % Simultaneous confidence interval for fitted curve.
           ud.crit = ud.crit_f;
       end
   else
      ud.crit = ud.crit_t;
   end

   % Update menu and bounds
   setconf(ud);
   updateplot(nlin_fig, ud);
end


function [y, deltay] = PredictionsPlusError(x,ud)
% Local function for Predicting a response with error bounds.
crit  = ud.crit;
R     = ud.R;
beta  = ud.beta;
model = ud.model;

% Predict Response
y = feval(model,beta,x);
   
% Calculate Error Bounds.
fdiffstep = eps(class(beta)).^(1/3);
J=zeros(length(y),numel(beta)); 
for k = 1:numel(beta)
    delta = zeros(size(beta));
    if (beta(k) == 0)
       nb = sqrt(norm(beta));
       delta(k) = fdiffstep * (nb + (nb==0));
    else
       delta(k) = fdiffstep * beta(k);
    end
    yplus = feval(model,beta+delta,x);
    J(:,k) = (yplus - y)/delta(k);
end

E = J / R;
if (ud.obsflag)
   deltay = sqrt(sum(E.*E,2) + 1)*crit;  % predict new observation
else
   deltay = sqrt(sum(E.*E,2))*crit;      % estimate regression function
end

function uihandles = MakeUIcontrols(xname,nlin_fig,newy,newdeltay,yname,avgx)
% Local function for Creating uicontrols for nlintool.
fcolor = get(nlin_fig,'Color');
yfieldp = [.01 .45 .10 .04];
uihandles.y_field(1) = uicontrol(nlin_fig,'Style','text','Units','normalized',...
         'Position',yfieldp + [0 .14 0 0],'String',num2str(double(newy)),...
         'ForegroundColor','k','BackgroundColor',fcolor);

uihandles.y_field(2) = uicontrol(nlin_fig,'Style','text','Units','normalized',...
         'Position',yfieldp,'String',num2str(double(newdeltay)),...
         'ForegroundColor','k','BackgroundColor',fcolor);

uicontrol(nlin_fig,'Style','text','Units','normalized',...
         'Position',yfieldp + [0 .07 -0.01 0],'String',' +/-',...
       'ForegroundColor','k','BackgroundColor',fcolor);

uicontrol(nlin_fig,'String', getString(message('stats:nlintool:uicontrol_Export')),...
              'Units','pixels','Position',[20 45 100 25],...
              'CallBack','nlintool(''output'')');
uicontrol(nlin_fig,'Style','Pushbutton','Units','pixels',...
               'Position',[20 10 100 25],'Callback','close','String',getString(message('stats:nlintool:uicontrol_Close')));

uicontrol(nlin_fig,'Style','text','Units','normalized',...
        'Position',yfieldp + [0 0.21 0 0],'BackgroundColor',fcolor,...
        'ForegroundColor','k','String',yname);

n = length(xname);
for k = 1:n
   xfieldp = [.18 + (k-0.5)*.80/n - 0.5*min(.5/n,.15) .09 min(.5/n,.15) .07];
   xtextp  = [.18 + (k-0.5)*.80/n - 0.5*min(.5/n,.18) .02 min(.5/n,.18) .05];
   uicontrol(nlin_fig,'Style','text','Units','normalized',...
        'Position',xtextp,'BackgroundColor',fcolor,...
        'ForegroundColor','k','String',xname{k});

   uihandles.x_field(k)  = uicontrol(nlin_fig,'Style','edit','Units','normalized',...
         'Position',xfieldp,'String',num2str(double(avgx(k))),...
         'BackgroundColor','white','CallBack',['nlintool(''edittext'',',int2str(k),')']);
end

% Create menu for controlling confidence bounds
f = uimenu('Label',getString(message('stats:nlintool:label_Bounds')), 'Position', 4, 'UserData','conf');
uimenu(f,'Label',getString(message('stats:nlintool:label_Simultaneous')), ...
       'Callback','nlintool(''conf'',1)', 'UserData',1);
uimenu(f,'Label',getString(message('stats:nlintool:label_NonSimultaneous')), ...
       'Callback','nlintool(''conf'',2)', 'UserData',2);
uimenu(f,'Label',getString(message('stats:nlintool:label_Curve')), 'Separator','on', ...
       'Callback','nlintool(''conf'',3)', 'UserData',3);
uimenu(f,'Label',getString(message('stats:nlintool:label_Observation')), ...
       'Callback','nlintool(''conf'',4)', 'UserData',4);
uimenu(f,'Label',getString(message('stats:nlintool:label_NoBounds')), 'Separator','on', ...
       'Callback','nlintool(''conf'',5)', 'UserData',5);

% -----------------------
function setconf(ud)
%SETCONF Update menus to reflect confidence bound settings

ma = get(findobj(gcf, 'Type','uimenu', 'UserData','conf'), 'Children');
set(ma, 'Checked', 'off');          % uncheck all menu items

% Check item 1 for simultaneous, 2 for non-simultaneous
hh = findobj(ma, 'Type', 'uimenu', 'UserData', 2-ud.simflag);
if (length(hh) == 1), set(hh, 'Checked', 'on'); end

% Check item 3 for curve, 4 for new observation
hh = findobj(ma, 'Type', 'uimenu', 'UserData', 3+ud.obsflag);
if (length(hh) == 1), set(hh, 'Checked', 'on'); end

% Check item 5 if we are omitting confidence bounds.  The remaining 
hh = findobj(ma, 'Type', 'uimenu', 'UserData', 5);
if ud.confflag
    set(ma, 'Enable','on');
else
    set(ma, 'Enable','off');
    set(hh, 'Checked', 'on', 'Enable','on');
end

% -----------------------
function ud = updateplot(nlin_fig, ud)
% Update plot after change in X value or confidence bound setting

x_field       = ud.x_field;
nlin_axes     = ud.nlin_axes;
xsettings     = ud.xsettings;
y_field1      = ud.y_field(1);
y_field2      = ud.y_field(2);
fitline       = ud.fitline;
reference_line= ud.reference_line;
last_axes     = ud.last_axes;  
n             = size(x_field,2);
xfit          = ud.xfit;

% Get current X values, updating from edit box if necessary
xrow  = xsettings(1,:);
lk = last_axes;
if lk~=0
    cx    = str2double(get(x_field(lk),'String'));
    xrow(lk) = cx;
end

[cy, dcy] = PredictionsPlusError(xrow,ud);

% If we need to update reference lines, lk will be non-empty
if lk~=0
   xsettings(:,lk) = cx(ones(41,1));
   set(reference_line(lk,1),'XData',cx*ones(2,1));
   set(reference_line(lk,2),'YData',[cy cy]);
end

ymax = zeros(n,1);
ymin = zeros(n,1);
ud.xsettings = xsettings;            

for idx = 1:n
   ith_x      = xsettings;
   ith_x(:,idx) = xfit(:,idx);
   [yfit, deltay] = PredictionsPlusError(ith_x,ud);

   if ~ud.confflag
      % No conf bounds wanted, so set their y data to a vector of NaNs
      % so they will not plot but the lines will be around for future use
      deltay = NaN;
   end
   set(nlin_axes(idx),'YlimMode','auto');
   set(fitline(1,idx),'Ydata',yfit);
   set(fitline(2,idx),'Ydata',yfit-deltay);
   set(fitline(3,idx),'Ydata',yfit+deltay);
   ylim = get(nlin_axes(idx),'YLim');
   ymin(idx) = ylim(1);
   ymax(idx) = ylim(2);
end         

ylims = [min(ymin), max(ymax)];
ud.last_axes = 0;
for ix = 1:n      
   set(nlin_axes(ix),'Ylim',ylims);
   set(reference_line(ix,1),'Ydata',ylims);
   set(reference_line(ix,2),'YData',[cy cy]);
end
set(nlin_fig,'Userdata',ud);       

set(y_field1,'String',num2str(double(cy)));
set(y_field2,'String',num2str(double(dcy)));


% ----------- helper to locate the axes under the cursor
function k = findaxes(fig, allaxes, ~, eventData)

k = [];
h = eventData.HitObject;
if h==fig
    return
end
h = ancestor(h,'axes');
if isempty(h)
    return
end
k = find(allaxes==h,1);

% -------------------------
function downFun(varargin)
nlin_fig = gcbf;
ud = get(nlin_fig,'Userdata');
set(nlin_fig,'WindowButtonMotionFcn',@(varargin) motionFun(1,varargin{:}));

k = findaxes(nlin_fig, ud.nlin_axes,varargin{:});
if isempty(k)
    return
end
ud.last_axes = k;
set(nlin_fig,'Pointer','crosshair');
cp = get(ud.nlin_axes(k),'CurrentPoint');
cx=cp(1,1);
xl = get(ud.nlin_axes(k),'Xlim');
maxx = xl(2);
minx = xl(1);
cx = max(minx, min(maxx, cx));

xrow  = ud.xsettings(1,:);
xrow(k) = cx;
ud.xsettings(:,k) = cx(ones(41,1));
[cy, dcy] = PredictionsPlusError(xrow,ud);

set(nlin_fig,'Userdata',ud);       

set(ud.x_field(k),'String',num2str(double(cx)));
set(ud.reference_line(k,1),'XData',cx*ones(2,1));
set(ud.reference_line(k,2),'YData',[cy cy]);

set(ud.y_field(1),'String',num2str(double(cy)));
set(ud.y_field(2),'String',num2str(double(dcy)));

set(nlin_fig,'WindowButtonUpFcn',@upFun);

% -------
function motionFun(flag,varargin)
nlin_fig = gcbf;
ud = get(nlin_fig,'Userdata');
k = findaxes(nlin_fig, ud.nlin_axes,varargin{:});
if isempty(k)
    return
end
xrange = get(ud.nlin_axes(k),'XLim');
newx(k) = str2double(get(ud.x_field(k),'String'));
maxx = xrange(2);
minx = xrange(1);
n = size(ud.x_field,2);

 if flag == 0,
    % Mouse was up, update cursor if over ref line
    % Check for data consistency in each axis
    if n > 1,
       yn = zeros(n,1);
       for idx = 1:n
          y = get(ud.reference_line(idx,2),'Ydata');
          yn(idx) = y(1);
       end
       % If data is inconsistent update all the plots.
       if diff(yn) ~= zeros(n-1,1)
          upFun(varargin{:});
       end
    end
    % Change cursor to plus sign when mouse is on top of vertical line.
     cursorstate = get(nlin_fig,'Pointer');
     cp = get(ud.nlin_axes(k),'CurrentPoint');
     cx = cp(1,1);
     fuzz = 0.02 * (maxx - minx);
     online = cx > newx(k) - fuzz & cx < newx(k) + fuzz;
     if online && strcmp(cursorstate,'arrow'),
         cursorstate = 'crosshair';
     elseif ~online && strcmp(cursorstate,'crosshair'),
         cursorstate = 'arrow';
     end
     set(nlin_fig,'Pointer',cursorstate);
     return
   elseif flag == 1,
     % Mouse was down, move reference line
     if ud.last_axes~=k
         return;
     end
     cp = get(ud.nlin_axes(k),'CurrentPoint');
     cx=cp(1,1);
     cx = max(minx, min(maxx, cx));

     xrow  = ud.xsettings(1,:);
     xrow(k) = cx;
     [cy, dcy] = PredictionsPlusError(xrow,ud);

     set(ud.x_field(k),'String',num2str(double(cx)));
     set(ud.reference_line(k,1),'XData', repmat(cx,2,1));
     set(ud.reference_line(k,2),'YData',[cy cy]);

     set(ud.y_field(1),'String',num2str(double(cy)));
     set(ud.y_field(2),'String',num2str(double(dcy)));
  end

% ------
function upFun(varargin)
nlin_fig = gcbf;
set(nlin_fig,'WindowButtonMotionFcn',@(varargin) motionFun(0,varargin{:}));

ud = get(nlin_fig,'Userdata');
n = size(ud.x_field,2);
p = get(nlin_fig,'CurrentPoint');
k = floor(1+n*(p(1)-0.18)/.80);

lk = ud.last_axes;
if lk==0
   return
end

if lk~=k
    % We are no longer over the axes where we started
    xrange = get(ud.nlin_axes(lk),'XLim');
    if k < lk
        set(ud.x_field(lk),'String',num2str(double(xrange(1))));
    elseif k > lk
        set(ud.x_field(lk),'String',num2str(double(xrange(2))));
    end
end
    
updateplot(nlin_fig, ud);


