function rstool(x,y,model,alpha,xname,yname)
%RSTOOL Multidimensional response surface fitting and visualization (RSM).
%   RSTOOL(X,Y,MODEL) opens an interactive GUI for fitting and visualizing a
%   polynomial response surface for a response variable Y as a function of the
%   multiple predictor variables in X.  Distinct predictor variables should
%   appear in different columns of X.  Y can be a vector, corresponding to a
%   single response, or a matrix, with columns corresponding to multiple
%   responses.  Y must have as many elements (or rows, if it is a matrix) as X
%   has rows.  RSTOOL displays a family of plots, one for each combination of
%   columns in X and Y.  Each plot represents the response surface as a
%   function of one predictor variable, with the remaining predictors held
%   fixed.  RSTOOL plots a 95% simultaneous confidence band for the fitted
%   response surface as two red curves on each plot.
%
%   The optional input MODEL controls the regression model.  By default, 
%   RSTOOL uses a linear additive model with a constant term.  MODEL can be 
%   any one of the following strings:
%
%     'linear'        Constant and linear terms (the default)
%     'interaction'   Constant, linear, and interaction terms
%     'quadratic'     Constant, linear, interaction, and squared terms
%     'purequadratic' Constant, linear, and squared terms
%
%   Alternatively, MODEL can be a matrix of model terms accepted by the 
%   X2FX function.  See X2FX for a description of this matrix and for
%   a description of the order in which terms appear.  You can use this
%   matrix to specify other models including ones without a constant term.
%
%   Drag the dashed blue reference line to examine predicted values.  
%   Specify a predictor by typing its value into the editable text field.
%
%   Use the pop-up menu to change the model.
%
%   Use the Export push button to export fitted coefficients and regression 
%   statistics to the base workspace.  Exported coefficients appear in the 
%   order defined by the X2FX function.
%
%   RSTOOL(X,Y,MODEL,ALPHA) plots 100(1-ALPHA)% confidence intervals for 
%   predictions.
%
%   RSTOOL(X,Y,MODEL,ALPHA,XNAME,YNAME) labels the axes using the names in
%   the strings XNAME and YNAME.  To label each subplot differently, XNAME 
%   and YNAME can be cell arrays of strings.
%
%   See also X2FX, NLINTOOL.
   
%   Copyright 1993-2014 The MathWorks, Inc.


if nargin > 2
    model = convertStringsToChars(model);
end

if nargin > 4
    xname = convertStringsToChars(xname);
end

if nargin > 5
    yname = convertStringsToChars(yname);
end

if nargin==0
    % No data given, so launch rstool using sample data
    msg = getString(message('stats:rstool:dialogText_NoArgs'));
    title = getString(message('stats:rstool:dialogTitle_NoArgs'));
    uiwait(msgbox(msg,title,'modal'));

    s = load('hald');
    x = s.ingredients;
    y = s.heat;
    model = 'interaction';
    alpha = 0.05;
    xname = {'X1' 'X2' 'X3' 'X4'};
    yname = 'Heat';
    nargs = 6;
else
    nargs = nargin;
end

if ~ischar(x) 
    action = 'start';
else
    action = x;
end

%On recursive calls get all necessary handles and data.
if ~strcmp(action,'start')
   if nargs == 2,
      flag = y;
   end
   lin_fig = gcbf;
   if ~isequal(get(lin_fig,'Tag'),'linfig');
      return
   end
   ud = get(lin_fig,'Userdata');
   
   rmse          = ud.rmse;
   residuals     = ud.residuals;
   alpha         = ud.alpha;
   modelpopup    = ud.modelpopup;
      
   beta          = ud.beta;
   model         = ud.model;
   x_field       = ud.x_field;
   lin_axes      = ud.lin_axes;
   xsettings     = ud.xsettings;
   y_field1      = ud.y_field(:,1);
   x             = ud.x;
   y             = ud.y;  
   n             = size(x,2);
   ny            = size(y,2);

   xrange         = zeros(n,2);
   for k = 1:n         
      xrange(k,1:2) = get(lin_axes(k,1),'XLim');
   end

   newy = zeros(numel(y_field1),1);
   for j = 1:length(newy)
      newy(j) = str2double(get(y_field1(j),'String'));
   end
end

switch action

case 'start'

if (nargs<2)
   error(message('stats:rstool:TooFewInputs'));
end
 
% Remove Nan if necessary
if (size(y,1) == 1), y = y(:); end
good = (sum(isnan(y),2) == 0) & (sum(isnan(x),2) == 0);
if (any(~good))
   y = y(good,:);
   x = x(good,:);
end
 
n = size(x,2);
ny = size(y,2);

if nargs < 4 || isempty(alpha)
   alpha = 0.05;
end

if nargs < 5 || isempty(xname)
   xname = cell(n,1);
   for i = 1:n
       xname{i} = getString(message('stats:rstool:xaxisLabel_Default', ...
                                    strjust(int2str(i),'left')));
   end
elseif ischar(xname) && (size(xname,1)==n)
   xname = cellstr(xname);
elseif iscellstr(xname) && (numel(xname)==n)
   % ok
else
   error(message('stats:rstool:BadXName'));
end
if nargs < 6 || isempty(yname)
   yname = cell(ny,1);
   for i = 1:ny
       yname{i} = getString(message('stats:rstool:yaxisLabel_Default', ...
                                    strjust(int2str(i),'left')));
   end
elseif ischar(yname) && (size(yname,1)==ny)
   yname = cellstr(yname);
elseif iscellstr(yname) && (numel(yname)==ny)
   % ok
else
   error(message('stats:rstool:BadYName'));
end

ud.usermodel = [];
if nargs < 3, model = ''; end
if isempty(model) || strcmp(model,'linear') || strcmp(model,'l')
   model = 'linear';
   mval  = 1;
   modelNameUI  = getString(message('stats:rstool:userModel_Linear'));
elseif strcmp(model,'purequadratic') || strcmp(model,'p'),
   mval  = 2;
   modelNameUI  = getString(message('stats:rstool:userModel_PureQuadratic'));
elseif strcmp(model,'interaction') || strcmp(model,'i'),
   mval  = 3;
   modelNameUI  = getString(message('stats:rstool:userModel_Interactions'));
elseif strcmp(model,'quadratic') || strcmp(model,'q'),
   mval  = 4;
   modelNameUI  = getString(message('stats:rstool:userModel_FullQuadratic'));
elseif ~ischar(model)
   mval  = 5;
   modelNameUI  = getString(message('stats:rstool:userModel_UserSpecified'));
   ud.usermodel = model;
else
   if any(exist(model,'file') == [2 3])
      % x2fx supports model as a function name
      mval  = 5;
      modelNameUI  = model;
      ud.usermodel = model;
   else
      error(message('stats:rstool:BadModel', model));
   end
end

design = x2fx(x,model);

% Fit response surface model design
if size(design,2) > size(x,1)
    error(message('stats:rstool:NotEnoughData', modelNameUI));
end

[Q, R] = qr(design,0);
if rcond(R) < 1E-12
    error(message('stats:rstool:NotEnoughData', modelNameUI));
end

[beta,~,residuals,~,~,rmse,crit] = endfit(Q,R,y,design,alpha);

ud.beta = beta;
ud.crit = crit;
ud.R    = R;
ud.model = model;
ud.good = good;

% Set positions of graphic objects
maxx = max(x);
minx = min(x);
xrange = maxx - minx;
maxx = maxx + 0.025 .* xrange;
minx = minx - 0.025 .* xrange;
xrange = 1.05*xrange;

lin_axes       = zeros(n,ny);
fitline        = zeros(3,n,ny);
reference_line = zeros(2,n,ny);

xfit      = xrange(ones(41,1),:)./40;
xfit(1,:) = minx;
xfit      = cumsum(xfit);

avgx      = mean(x);

xsettings = avgx(ones(42,1),:);

ud.xfit = xfit;
ud.xsettings = xsettings;


lin_fig = figure('Units','Normalized','Interruptible','on',...
             'Position',[0.05 0.35 0.90 0.5],...
             'NumberTitle','off', 'IntegerHandle','off', ...
             'Name', ...
             getString(message('stats:rstool:figureTitle',modelNameUI)), ...
             'Tag','linfig','ToolBar','none');
set(0,'CurrentFigure',lin_fig);

% Remove brushing/linking tools that conflict with this GUI
delete(uigettool(lin_fig,'Exploration.Brushing'))
delete(findall(lin_fig,'Tag','figDataManagerBrushTools'))
delete(findall(lin_fig,'Tag','figBrush'))
delete(uigettool(lin_fig,'DataManager.Linking'))
delete(findall(lin_fig,'Tag','figLinked'))

% Set up axes
xtmp = 0:1;
for k = 1:n
   for j = 1:ny
      % Create an axis for each pair of input (x) and output (y) variables
      axisp   = [.18+(k-1)*.80/n  .22+(j-1)*.76/ny  .80/n  .76/ny];
      lin_axes(k,j) = axes;
      set(lin_axes(k,j),'XLim',[minx(k) maxx(k)],'Box','on',...
                        'NextPlot','add',...
                        'Position',axisp,'Gridlinestyle','none');
      if k>1, set(lin_axes(k,j),'Yticklabel',[]); end
      if j>1, set(lin_axes(k,j),'Xticklabel',[]); end

      % Add curves
      fitline(1:3,k,j) = plot(xtmp,xtmp,'g-', xtmp,xtmp,'r--', xtmp,xtmp,'r--');
      
      % Add reference Lines
      reference_line(1,k,j) = plot(xtmp, xtmp,'--');
      reference_line(2,k,j) = plot(xtmp, xtmp,':');
   end
end

ud.fitline = fitline;
ud.lin_axes = lin_axes;
ud.reference_line = reference_line;

uihandles = MakeUIcontrols(xname,lin_fig,yname,avgx,mval,ud.usermodel);

ud.export = uihandles.export;
ud.modelpopup = uihandles.modelpopup;
ud.x_field = uihandles.x_field; 
ud.y_field = uihandles.y_field;

% Update curves and reference lines on all graphs
calcy(lin_fig,n,ny,ud);

ud.texthandle = [];
ud.residuals = residuals;
ud.reference_line = reference_line;
ud.last_axes = 0;
ud.x = x;
ud.y = y;
ud.rmse = rmse;
ud.alpha = alpha;
set(lin_fig,'UserData',ud,'HandleVisibility','callback',...
            'BusyAction','queue', ...
            'WindowButtonMotionFcn',@(varargin) motionFun(0,varargin{:}), ...
            'WindowButtonDownFcn',@downFun,...
            'WindowButtonUpFcn',@upFun,'Interruptible','on');
% Finished with plot startup function.

case 'edittext',
   cx    = str2double(get(x_field(flag),'String'));  
   if isnan(cx)
       set(x_field(flag),'String',num2str(xsettings(1,flag)));
       % Create Bad Settings Warning Dialog.
       warndlg(getString(message('stats:rstool:warndlg_NumbersOnly')),'RSTOOL','modal');
       return
   end  
   
   xl = get(lin_axes(flag,1),'Xlim');
   if cx < xl(1) || cx > xl(2)
       % Create Bad Settings Warning Dialog.
       warndlg(getString(message('stats:rstool:warndlg_Range')),'RSTOOL','modal');
       set(x_field(flag),'String',num2str(xsettings(1,flag)));
       return
   end
   
   xsettings(:,flag) = cx(ones(42,1));
   ud.xsettings = xsettings;            

   % Update graph
   calcy(lin_fig, n, ny, ud, flag, cx);

   set(lin_fig,'Userdata',ud);       

case 'output',
     bmf = get(lin_fig,'WindowButtonMotionFcn');
     bdf = get(lin_fig,'WindowButtonDownFcn');
     set(lin_fig,'WindowButtonMotionFcn','');
     set(lin_fig,'WindowButtonDownFcn','');
 
    checkLabels = {getString(message('stats:rstool:checkbox_Coefficients')), ...
                   getString(message('stats:rstool:checkbox_RMSE')), ...
                   getString(message('stats:rstool:checkbox_Residuals'))};
    defaultVarNames = {'beta', 'rmse', 'residuals'};
 
    if (all(ud.good))
        fullresid = residuals;
    else
        fullresid = NaN(length(ud.good), size(residuals,2));
        fullresid(ud.good,:) = residuals;
    end
    
    items = {beta, rmse, fullresid};
     
    export2wsdlg(checkLabels, defaultVarNames, items, getString(message('stats:rstool:dialogTitle_Export')));
    
     set(lin_fig,'WindowButtonMotionFcn',bmf);
     set(lin_fig,'WindowButtonDownFcn',bdf);

case 'changemodel',
   cases = get(modelpopup,'Value');
   if cases == 1
      model = 'linear';
      modelNameUI  = getString(message('stats:rstool:userModel_Linear'));
   elseif cases == 2
      model = 'purequadratic';
      modelNameUI  = getString(message('stats:rstool:userModel_PureQuadratic'));
   elseif cases == 3
      model = 'interaction';
      modelNameUI  = getString(message('stats:rstool:userModel_Interactions'));
   elseif cases == 4
      model = 'quadratic';
      modelNameUI  = getString(message('stats:rstool:userModel_FullQuadratic'));
   elseif cases == 5
      if ischar(model)
         if isempty(ud.usermodel)
            disp(getString(message('stats:rstool:commandWindowText_CallRSTOOL')));
            disp(getString(message('stats:rstool:commandWindowText_HelpX2FX')));         
            disp(getString(message('stats:rstool:commandWindowText_FittingLinear')));
            model = 'linear';
            modelNameUI  = getString(message('stats:rstool:userModel_Linear'));
            set(modelpopup,'Value',1);
         else
            model = ud.usermodel;
			if ischar(ud.usermodel)
               % This is a file name holding a design matrix.
               modelNameUI = model;
			   set(lin_fig,'Name', getString(message('stats:rstool:figureTitle',modelNameUI)));
            else
               modelNameUI  = getString(message('stats:rstool:userModel_UserSpecified'));
	           set(lin_fig,'Name', getString(message('stats:rstool:figureTitle',modelNameUI)));
			end
         end
      end
   end
    
   % Fit response surface model design
   design = x2fx(x,model);
   [Q, R] = qr(design,0);
 
   if (size(R,1) < size(R,2)) || rcond(R) < 1E-12
       % Create Model Warning Figure.
       s = getString(message('stats:rstool:warndlg_InsufficientData',modelNameUI));
       warndlg(s,'RSTOOL','modal');
       mval = get(modelpopup,'Userdata');
       set(modelpopup,'Value',mval);

       if mval == 1
          set(lin_fig,'Name',getString(message('stats:rstool:figureTitle', ...
              getString(message('stats:rstool:userModel_Linear')))));
       elseif mval == 2
          set(lin_fig,'Name',getString(message('stats:rstool:figureTitle', ...
              getString(message('stats:rstool:userModel_PureQuadratic')))));
       elseif mval == 3
          set(lin_fig,'Name',getString(message('stats:rstool:figureTitle', ...
              getString(message('stats:rstool:userModel_Interactions')))));
       elseif mval == 4
          set(lin_fig,'Name',getString(message('stats:rstool:figureTitle', ...
              getString(message('stats:rstool:userModel_FullQuadratic')))));
       end
       drawnow;
       return;
    end
   
    [beta,~,residuals,~,~,rmse,crit] = endfit(Q,R,y,design,alpha);
    
    ud.crit      = crit;
    ud.rmse      = rmse;
    ud.R         = R;
    ud.beta      = beta;
    ud.model     = model;
    ud.residuals = residuals;
    set(lin_fig,'Userdata',ud);   

    % Update graph
    calcy(lin_fig, n, ny, ud);

   if cases == 1
      set(lin_fig,'Name',getString(message('stats:rstool:figureTitle', ...
         getString(message('stats:rstool:userModel_Linear')))));
   elseif cases == 2
      set(lin_fig,'Name',getString(message('stats:rstool:figureTitle', ...
         getString(message('stats:rstool:userModel_PureQuadratic')))));
   elseif cases == 3
      set(lin_fig,'Name',getString(message('stats:rstool:figureTitle', ...
         getString(message('stats:rstool:userModel_Interactions')))));
   elseif cases == 4
      set(lin_fig,'Name',getString(message('stats:rstool:figureTitle', ...
         getString(message('stats:rstool:userModel_FullQuadratic')))));
   end
   set(modelpopup,'Userdata',cases);
end

function uihandles = MakeUIcontrols(xstr,lin_fig,ystr,avgx,mval,usermodel)
% Local function for Creating uicontrols for rstool.
n = length(xstr);
fcolor = get(lin_fig,'Color');
for k = 1:n
   xfieldp = [.18 + (k-0.5)*.80/n - 0.5*min(.5/n,.15) .09 min(.5/n,.15) .07];
   xtextp  = [.18 + (k-0.5)*.80/n - 0.5*min(.5/n,.18) .02 min(.5/n,.18) .05];
   uicontrol(lin_fig,'Style','text','Units','normalized',...
        'Position',xtextp,'BackgroundColor',fcolor,...
        'ForegroundColor','k','String',xstr{k});

   uihandles.x_field(k)  = uicontrol(lin_fig,'Style','edit',...
         'Units','normalized','Position',xfieldp,'String',num2str(avgx(k)),...
         'BackgroundColor','white',...
         'CallBack',['rstool(''edittext'',',int2str(k),')']);
end

ny = length(ystr);
oldu = get(uihandles.x_field(1),'FontUnits');
set(uihandles.x_field(1), 'FontUnits', 'normalized');
halfsize = 0.6 * get(uihandles.x_field(1), 'FontSize') * .07;
set(uihandles.x_field(1), 'FontUnits', oldu);
yfieldp = [.01 0 .10 2*halfsize];

for j = 1:ny
   ycenter = .22+(j-.5)*.76/ny;

   % y value field
   yfieldp(2) = ycenter + halfsize;
   uihandles.y_field(j,1) =uicontrol(lin_fig,'Style','text',...
         'Units','normalized', 'Position',yfieldp,...
         'String','', 'ForegroundColor','k', 'BackgroundColor',fcolor);

   % y delta field
   yfieldp(2) = ycenter - 3 * halfsize;
   uihandles.y_field(j,2) =uicontrol(lin_fig,'Style','text',...
         'Units','normalized', 'Position',yfieldp,...
         'String','', 'ForegroundColor','k', 'BackgroundColor',fcolor);

   % plus or minus field
   yfieldp(2) = ycenter - halfsize;
   uicontrol(lin_fig,'Style','text','Units','normalized',...
          'Position',yfieldp, 'String',' +/-',...
          'ForegroundColor','k','BackgroundColor',fcolor);

   % y name field
   yfieldp(2) = ycenter + 3 * halfsize;
   uicontrol(lin_fig,'Style','text','Units','normalized',...
          'Position',yfieldp, 'BackgroundColor',fcolor,...
          'ForegroundColor','k','String',ystr(j));
end
   
uihandles.export     = uicontrol(lin_fig,'Style','Pushbutton', ...
              'String',getString(message('stats:rstool:button_Export')), ...
              'Units','normalized','Position', ...
              [0.01 0.14 0.13 0.05], 'CallBack','rstool(''output'')');

modellist = { getString(message('stats:rstool:userModel_Linear')), ...
              getString(message('stats:rstool:userModel_PureQuadratic')), ...
              getString(message('stats:rstool:userModel_Interactions')), ...
              getString(message('stats:rstool:userModel_FullQuadratic')) };
if ~isempty(usermodel)
   modellist = {modellist{:}, getString(message('stats:rstool:userModel_UserSpecified'))};
end

uihandles.modelpopup = uicontrol(lin_fig,'Style','Popup','String',...
        modellist, ...
          'Value',mval,'BackgroundColor','w',...
              'Units','normalized','Position',[0.01 0.08 0.13 0.05],...
              'CallBack','rstool(''changemodel'')','Userdata',mval);

uicontrol('Style','Pushbutton','Units','normalized',...
         'Position',[0.01 0.02 0.13 0.05],'Callback','close', ...
         'String',getString(message('stats:rstool:button_Close')));

% ---- helper to create or update prediction curves
function deltay = calcy(lin_fig,n,ny,ud, xk, xval)

crit      = ud.crit;
beta      = ud.beta;
R         = ud.R;
model     = ud.model;
xsettings = ud.xsettings;
xfit      = ud.xfit;
fitline   = ud.fitline;
lin_axes  = ud.lin_axes;
reference_line = ud.reference_line;

% Get info stored in figure
yextremes = zeros(n,ny,2);

for k = 1:n
   % Calculate y values for fitted line plot.
   ith_x      = xsettings;
   xi         = xfit(:,k);
   ith_x(1:41,k) = xi;
   if (nargin > 4)
      ith_x(42, xk) = xval;  % otherwise row 42 stays at current setting
   end
   
   % Calculate y values for confidence interval lines.
   xpred = x2fx(ith_x,model);
   yfit  = xpred(1:41,:)*beta;
   newy  = xpred(42,:)*beta;
   E     = xpred(1:41,:)/R;
   tmp   = sqrt(sum(E.*E,2));
   dy    = repmat(tmp,1,length(crit)) .* repmat(crit, length(tmp), 1);

   for j = 1:ny
      % Plot prediction line with confidence intervals
      yfitj = yfit(:,j);
      dyj = dy(:,j);
      t1 = yfitj - dyj;
      t2 = yfitj + dyj;
      set(lin_fig,'CurrentAxes',lin_axes(k,j));
      set(fitline(1,k,j),'Xdata',xi,'Ydata',yfitj);
      set(fitline(2,k,j),'Xdata',xi,'Ydata',t1);
      set(fitline(3,k,j),'Xdata',xi,'Ydata',t2);

      % No x ticks right near the end, they might overlap the next axis
      xl = get(lin_axes(k,j),'Xlim');
      xr = diff(xl);
      xt = get(lin_axes(k,j),'XTick');
      lowtick = xl(1) + .1*xr;
      hitick  = xl(2) - .1*xr;
      xt = xt(xt>lowtick & xt<hitick);
      set(lin_axes(k,j),'Xtick',xt);
      
      % Calculate data for vertical reference lines, allow for dyj=NaN
      yextremes(k,j,1) = min(min(yfitj), min(t1));
      yextremes(k,j,2) = max(max(yfitj), max(t2));

      if (k == 1)
         E = xpred(42,:)/R;
         deltay = sqrt(E*E')*crit;
      end
   end

end  % End of the plotting loop over all the axes.

ymin = min(yextremes(:,:,1), [], 1);
ymax = max(yextremes(:,:,2), [], 1);

for k = 1:n
   t1 = xsettings([1 1], k);
   t2 = get(lin_axes(k,1),'XLim');
   for j = 1:ny
      set(lin_axes(k,j), 'XLim', t2);
      set(reference_line(1,k,j),'Xdata', t1);
      set(reference_line(2,k,j),'XData', t2);
   end
end
for j = 1:ny
   t1 = [ymin(j) ymax(j)];
   t2 = newy([j j]);
   for k = 1:n
      set(lin_axes(k,j), 'YLim', t1);
      set(reference_line(1,k,j),'Ydata', t1);
      set(reference_line(2,k,j),'YData', t2);
   end
end

% Update labels
for j = 1:ny
   set(ud.y_field(j,1),'String',num2str(newy(j)));
   set(ud.y_field(j,2),'String',num2str(deltay(j)));
end


% ------- helper to compute fit results
function [b,yhat,res,p,df,rmse,crit]=endfit(Q,R,y,design,alpha)
% complete fit after the QR decomposition 

b = R\(Q'*y);

yhat = design*b;
res = y - yhat;

p = length(b);
df = max(size(y,1)) - p;

if (df > 0)
   rmse = sqrt(sum(res.*res)./df);
   % Not a prediction interval - Scheffe param = p
   crit = sqrt(p*finv(1 - alpha,p,df))*rmse;
else
   rmse = NaN;
   crit = NaN;
end

% -------- helper to update reference lines
function updateref(lin_fig, k, j, ud, maxx, minx)
% Update reference line on axis (k,j) using current settings

cp = get(ud.lin_axes(k,j),'CurrentPoint');
cx = min(maxx, max(minx, cp(1,1)));

xrow  = ud.xsettings(1,:);
xrow(k) = cx;
      
drow = x2fx(xrow, ud.model);
yrow = drow * ud.beta;
E = drow / ud.R;
deltay = sqrt(E*E') * ud.crit;

ud.xsettings(:,k) = cx(ones(42,1));
set(lin_fig,'Userdata',ud);       
set(ud.x_field(k),'String',num2str(cx));
set(ud.y_field(j,1),'String',num2str(yrow(j)));
set(ud.y_field(j,2),'String',num2str(deltay(j)));
set(ud.reference_line(1,k,j),'XData',cx*ones(2,1));
set(ud.reference_line(2,k,j),'YData',[yrow(j); yrow(j)]);

% ----------- helper to locate the axes under the cursor
function [k,j] = findaxes(fig, allaxes, ~, eventData)

k = [];
j = [];
h = eventData.HitObject;
 
if h==fig
    return
end
h = ancestor(h,'axes');
if isempty(h)
    return
end
[k,j] = find(allaxes==h,1);

% ---------- Callbacks in response to mouse movement and clicks
function downFun(varargin)
lin_fig = gcbf;
if ~isequal(get(lin_fig,'Tag'),'linfig');
   return
end
ud = get(lin_fig,'Userdata');

[k,j] = findaxes(lin_fig, ud.lin_axes, varargin{:});
if isempty(k)
    return
end

% Click over axes k. Change the motion function and move the ref line.
ud.last_axes = k;
set(lin_fig,'Pointer','crosshair',...
    'WindowButtonMotionFcn',@(varargin) motionFun(1,varargin{:}));
xrange = get(ud.lin_axes(k,j),'Xlim');
maxx = xrange(2);
minx = xrange(1);
updateref(lin_fig, k, j, ud, maxx, minx);

% ----
function motionFun(flag,varargin)
lin_fig = gcbf;
if ~isequal(get(lin_fig,'Tag'),'linfig');
   return
end
ud = get(lin_fig,'Userdata');

[k,j] = findaxes(lin_fig, ud.lin_axes, varargin{:});
if isempty(k)
    return
end

% Moving over axes k. Set this as current axes.
newx = str2double(get(ud.x_field(k),'String'));
xlim = get(ud.lin_axes(k,1),'XLim');
maxx = xlim(2);
minx = xlim(1);
set(lin_fig,'CurrentAxes',ud.lin_axes(k,j));

if flag == 0 
    % Button is up, set cursor differently if we are on a ref line
    cursorstate = get(lin_fig,'Pointer');
    cp = get(ud.lin_axes(k,j),'CurrentPoint');
    cx = cp(1,1);
    fuzz = 0.02 * (maxx - minx);
    online = cx > newx - fuzz & cx < newx + fuzz;
    if online && strcmp(cursorstate,'arrow'),
        cursorstate = 'crosshair';
    elseif ~online && strcmp(cursorstate,'crosshair'),
        cursorstate = 'arrow';
    end
    set(lin_fig,'Pointer',cursorstate);
    return
    
elseif flag == 1
    % Button is down, move ref line to new location
    if ud.last_axes~=k, return; end
    updateref(lin_fig, k, j, ud, maxx, minx);
    
end

% ----
function upFun(varargin)
lin_fig = gcbf;
if ~isequal(get(lin_fig,'Tag'),'linfig');
    return
end

% Mouse up, reset motion function
set(lin_fig,'WindowButtonMotionFcn',@(varargin) motionFun(0,varargin{:}), ...
    'Pointer','arrow');

ud = get(lin_fig,'Userdata');
lk = ud.last_axes;
if lk==0
    return % unexpected, we have no original axes to update
end

% Find out where we are
n = size(ud.x,2);
ny = size(ud.y,2);
[k,~] = findaxes(lin_fig, ud.lin_axes, varargin{:});
if isempty(k)
    % If we slipped above/below axes, figure out from x location
    p = get(lin_fig,'CurrentPoint');
    k = floor(1+n*(p(1)-0.18)/.80);
end

if k~=lk   % we are not over the axes we started from
    xlim = get(ud.lin_axes(lk,1),'XLim');
    
    if k < lk
        set(ud.x_field(lk),'String',num2str(xlim(1)));
    elseif k > lk
        set(ud.x_field(lk),'String',num2str(xlim(2)));
    end
end

% Update x field and graph
cx    = str2double(get(ud.x_field(lk),'String'));
ud.xsettings(:,lk) = cx(ones(42,1));

calcy(lin_fig, n, ny, ud, lk, cx);

ud.last_axes = 0;
set(lin_fig,'Userdata',ud);
