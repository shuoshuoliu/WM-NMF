function randtool(action,flag)
%RANDTOOL Demonstration of many random number generators.
%   RANDTOOL creates a histogram of random samples from many distributions.
%   This is an example that displays histograms of random samples from
%   the distributions in Statistics and Machine Learning Toolbox.
%
%   Change the parameters of the distribution by typing in a new
%   value or by moving a slider.
%
%   Export the current sample to a variable in the workspace by pressing
%   the export button.
%
%   Change the sample size by typing any positive integer in the
%   sample size edit box.
%
%   Change the distribution type using the popup menu.

%   Copyright 1993-2015 The MathWorks, Inc.


if nargin > 0
    action = convertStringsToChars(action);
end

if nargin < 1
    action = 'start';
end

%On recursive calls get all necessary handles and data.
if ~strcmp(action,'start')   
   childList = allchild(0);
   set(0,'CurrentFigure',childList(strcmp({childList.Tag},'randfig'))); 
   rand_fig = childList(childList == gcf);
   ud = get(rand_fig,'UserData');

   switch action
      case {'setpfield', 'setpslider', 'setphi', 'setplo'}
         ud = feval(action,flag,ud);
      case 'changedistribution', 
         ud = changedistribution(ud);
      case 'setsample',
         ud=updategui(ud,[]); 
      case 'output',
         outrndnum(ud); 
   end
   set(rand_fig,'UserData',ud);
end

% Initialize all GUI objects. Plot Normal CDF with sliders for parameters.
if strcmp(action,'start'),
   % Set positions of graphic objects
   axisp   = [.11 .35 .86 .58];

   pos = cell(7,3);
   pos{5,1} = [.07 .135 .12 .04];   % text
   pos{1,1} = [.19 .145 .12 .04];   % parameter
   pos{3,1} = [.19 .085 .12 .04];   % lower bound
   pos{4,1} = [.32 .085 .03 .16];   % slider
   pos{2,1} = [.19 .205 .12 .04];   % upper bound
   pos{7,1} = [.07 .075 .12 .04];   % lower bound label
   pos{6,1} = [.07 .195 .12 .04];   % upper bound label
   
   pos{5,2} = [.38 .135 .12 .04];   % text
   pos{1,2} = [.50 .145 .12 .04];   % parameter
   pos{3,2} = [.50 .085 .12 .04];   % lower bound
   pos{4,2} = [.63 .085 .03 .16];   % slider
   pos{2,2} = [.50 .205 .12 .04];   % upper bound
   pos{7,2} = [.38 .075 .12 .04];   % lower bound label
   pos{6,2} = [.38 .195 .12 .04];   % upper bound label
   
   pos{5,3} = [.69 .135 .12 .04];   % text
   pos{1,3} = [.81 .145 .12 .04];   % parameter
   pos{3,3} = [.81 .085 .12 .04];   % lower bound
   pos{4,3} = [.94 .085 .03 .16];   % slider
   pos{2,3} = [.81 .205 .12 .04];   % upper bound
   pos{7,3} = [.69 .075 .12 .04];   % lower bound label
   pos{6,3} = [.69 .195 .12 .04];   % upper bound label

   ud.dists = statguidists;
   dfltDist = find(strcmpi('Normal',{ud.dists.code}));
   
   % Set axis limits and data
   rand_fig = figure('Tag', 'randfig', 'Position', [304 380 625 605], ...
                'Visible', 'off','NumberTitle','off','IntegerHandle','off',...
                'Name',getString(message('stats:randtool:TitleRandtool')),...
                'ToolBar','none');
   figcolor  = get(rand_fig,'Color');
   set(rand_fig,'Units','Normalized');
   rand_axes = axes;
   xrange = getxdata(dfltDist,ud);
   set(rand_axes,...
      'Position',axisp,'XLim',xrange,'Box','on','Tag','randaxes');
   set(rand_fig,'UserData',ud);

   % Define graphics objects
   for idx = 1:3
       nstr = int2str(idx);
       ud.pfhndl(idx) = uicontrol('Style','edit','Units','normalized',...
           'Position',pos{1,idx},...
           'String',num2str(ud.dists(dfltDist).parameters(2-rem(idx,2))),...
           'BackgroundColor','white',...
           'Callback',['randtool(''setpfield'',',nstr,')'],...
           'Tag',['pfhndl' num2str(idx)]);
         
       ud.hihndl(idx) = uicontrol('Style','edit','Units','normalized',...
           'Position',pos{2,idx},...
           'String',num2str(ud.dists(dfltDist).phi(2-rem(idx,2))),...
           'BackgroundColor','white',...
           'Callback',['randtool(''setphi'',',nstr,')'],...
           'Tag',['hihndl' num2str(idx)]);
         
       ud.lohndl(idx) = uicontrol('Style','edit','Units','normalized',...
           'Position',pos{3,idx},...
           'String',num2str(ud.dists(dfltDist).plo(2-rem(idx,2))),...
           'BackgroundColor','white',... 
           'Callback',['randtool(''setplo'',',nstr,')'],...
           'Tag',['lohndl' num2str(idx)]);

       ud.pslider(idx) = uicontrol('Style','slider','Units','normalized',...
           'Position',pos{4,idx},...
           'Value',ud.dists(dfltDist).parameters(2-rem(idx,2)),...
           'Max',ud.dists(dfltDist).phi(2-rem(idx,2)),...
           'Min',ud.dists(dfltDist).plo(2-rem(idx,2)),...
           'Callback',['randtool(''setpslider'',',nstr,')'],...
           'Tag',['pslider' num2str(idx)]);

       ud.ptext(idx) = uicontrol('Style','text','Units','normalized',...
           'Position',pos{5,idx},...
           'BackgroundColor',figcolor,'ForegroundColor','k',...
           'String',ud.dists(dfltDist).paramnames{2-rem(idx,2)},...
           'Tag',['ptext' num2str(idx)]);
       
       ud.lowerboundtext(idx) = uicontrol('Style','text','Units','normalized', ...
         'Position', pos{7,idx},...
         'ForegroundColor','k','BackgroundColor',figcolor,...
         'String', getString(message('stats:randtool:LabelLowerBound')),...
         'Tag',['lowerboundtext' num2str(idx)]); 
   
       ud.upperboundtext(idx) = uicontrol('Style','text','Units','normalized', ...
         'Position', pos{6,idx},'ForegroundColor','k',...
         'BackgroundColor',figcolor,...
         'String', getString(message('stats:randtool:LabelUpperBound')),...
         'Tag',['upperboundtext' num2str(idx)]); 

       setincrement(ud.pslider(idx), ud.dists(dfltDist).intparam(2-rem(idx,2)));
   end      
   
   enableParams(ud, 3, 'off');

   distNameList = {ud.dists.name};
   ud.popup=uicontrol('Style','Popup','String',distNameList,...
        'Units','normalized','Position',[.28 .945 .25 .04],...
        'UserData','popup','Value',dfltDist,...
        'BackgroundColor','white',...
        'Callback','randtool(''changedistribution'')','Tag','popup');

   ud.samples_field = uicontrol('Style','edit','Units','normalized',...
         'Position',[.71 .945 .15 .04], ...
         'String',num2str(100),...
         'BackgroundColor','white',...
         'Callback','randtool(''setsample'',1)','Tag','samplesfield');

   uicontrol('Style','Pushbutton',...
         'Units','normalized','Position',[.70 .02 .13 .04],...
         'Callback','randtool(''setsample'',1)',...
         'String',getString(message('stats:randtool:ButtonResample')),...
         'Tag','resamplebutton');

   uicontrol('Style','Pushbutton','Units','normalized',...
         'Position',[.84 .02 .13 .04],...
         'Callback','randtool(''output'',2);',...
         'String',getString(message('stats:randtool:ButtonExport')),...
         'Tag','exportbutton');

   ud = updategui(ud,[]);
   set(rand_fig,'Visible','on'); % need to do this before calling placetitlebar
   placetitlebar(rand_fig);
   set(rand_fig,'UserData',ud,'HandleVisibility','callback',...
       'InvertHardcopy', 'on', 'PaperPositionMode', 'auto')
   drawnow;
end % End of initialization.

% END OF randtool MAIN FUNCTION.

% ==== Begin helper functions ====

%------------------------------------------------------------------------------
% Supply x-axis range for each distribution. GETXDATA
function xrange = getxdata(popupvalue,ud)
phi = ud.dists(popupvalue).phi;
plo = ud.dists(popupvalue).plo;
switch ud.dists(popupvalue).rvname
    case 'betarv', % Beta 
       xrange  = [0 1];
    case 'binorv', % Binomial 
       xrange  = [-.5 phi(1)+.5];
    case 'burrrv', % Burr
       if phi(2) < 1
           xrange = icdf('burr', [0.001 0.95], (phi(1)+plo(1))/2, (phi(2)+plo(2))/2, (phi(3)+plo(3))/2 );
       else
           xrange = icdf('burr', [0.001 0.995], (phi(1)+plo(1))/2, (phi(2)+plo(2))/2, (phi(3)+plo(3))/2 );
       end
	case 'chi2rv', % Chi-square
       xrange  = [0 phi + 4 * sqrt(2 * phi)];
    case 'unidrv', % Discrete Uniform
       xrange  = [0.5 phi+.5];
    case 'exprv', % Exponential
       xrange  = [0 4*phi];
    case 'evrv', % Extreme Value
       xrange = [plo(1)-5*phi(2), phi(1)+2*phi(2)];
    case 'frv', % F 
       xrange  = [0 finv(0.995,plo(1),plo(1))];
    case 'gamrv', % Gamma
       hixvalue = phi(1) * phi(2) + 4*sqrt(phi(1) * phi(2) * phi(2));
       xrange  = [0 hixvalue];
    case'gevrv', % Generalized Extreme Value
       loxvalue = gevinv(0.01,plo(1),phi(2),plo(3));
       hixvalue = gevinv(0.99,phi(1),phi(2),phi(3));
       xrange  = [loxvalue-0.5 hixvalue+.5];       
    case 'gprv', % Generalized Pareto
       hixvalue = gpinv(0.99,phi(1),phi(2),phi(3));
       xrange  = [plo(3)-0.5 hixvalue+.5];       
    case 'georv', % Geometric
       hixvalue = geoinv(0.99,plo(1));
       xrange  = [-0.5 hixvalue+.5];   
    case 'hnrv', % Half Normal
       hixvalue = icdf('halfnorm',0.99,phi(1),phi(2));
       xrange  = [plo(1)-0.5 hixvalue+.5];  
    case 'hygerv', % Hypergeometric
       xrange  = [-0.5 phi(1)+.5];
    case 'lognrv', % Lognormal
       xrange = [0 logninv(0.99,phi(1),phi(2))];
    case 'nbinrv', % Negative Binomial
       xrange = [-0.5 nbininv(0.99,phi(1),plo(2))+.5];
    case 'ncfrv', % Noncentral F
       xrange = [0 phi(3)+30];
    case 'nctrv', % Noncentral T
       xrange = [phi(2)-14 phi(2)+14];
    case 'ncx2rv', % Noncentral Chi-square
       xrange = [0 phi(2)+30];
    case 'normrv', % Normal
       xrange   = [plo(1) - 3 * phi(2) phi(1) + 3 * phi(2)];
    case 'poissrv', % Poisson
      xrange  = [-0.5 4*phi(1)+.5];
    case 'raylrv', % Rayleigh
       xrange = [0 raylinv(0.995,phi(1))];
    case 'trv', % T
       lowxvalue = tinv(0.005,plo(1));
       xrange  = [lowxvalue -lowxvalue];
    case 'unifrv', % Uniform
       xrange  = [plo(1) phi(2)];
    case 'weibrv', % Weibull
       xrange  = [0 wblinv(0.995,plo(1),plo(2))];
end


%-----------------------------------------------------------------------------
% Determine validity of value with respect to min and max
function valid = okwithminmax(cv, pmin, pmax, popupvalue, fieldnum, ud)

if ~isreal(cv) || ~isreal(pmin) || ~isreal(pmax)
    valid = false;

% All parameters may be in the open interval (pmin, pmax)
elseif (pmin < cv) && (cv < pmax)
    valid = true;
    
else
    valid = false;
    rvname = ud.dists(popupvalue).rvname;
    paramname = lower(ud.dists(popupvalue).paramnames{fieldnum});
    
    % Binomial p may also be in the closed interval [pmin, pmax]
    if  (isequal(rvname, 'binorv') && isequal(paramname,'prob')) 
        if (pmin <= cv) &&  (cv <= pmax)
            valid = true;
        end
        
    % NC Chi-sq and F delta may also be in half-open interval [pmin, pmax)
    elseif  (isequal(rvname, 'ncx2rv') && isequal(paramname,'delta')) ...
         || (isequal(rvname, 'ncfrv')  && isequal(paramname,'delta'))
        if (pmin <= cv) &&  (cv < pmax)
            valid = true;
        end
        
    % Hypergeometric may also be in half-open interval [pmin, pmax)
    elseif  (isequal(rvname, 'hygerv') && ...
                (isequal(paramname,'k') || isequal(paramname,'n')))
        if (pmin <= cv) &&  (cv < pmax)
            valid = true;
        end

    end
end


%-----------------------------------------------------------------------------
% Check that new slider positions are legal
function ok = slidervaluesok(popupvalue, cv, fieldno, ud)
ok = true;

% For uniform, pvalue "max" must be greater than pvalue "min"
if isequal(ud.dists(popupvalue).rvname, 'unifrv')
    paramname = lower(ud.dists(popupvalue).paramnames{fieldno});
    if isequal(paramname, 'min')
        pv = ud.dists(popupvalue).parameters(2);
        if cv >= pv
            ok = false;
        end
    else %if isequal(paramname, 'max')
        pv = ud.dists(popupvalue).parameters(1);
        if cv <= pv
            ok = false;
        end
    end
end

% For hypergeometric, N and K must each be no more than M
if isequal(ud.dists(popupvalue).rvname, 'hygerv')
    paramname = lower(ud.dists(popupvalue).paramnames{fieldno});
    switch paramname
        case 'm',
            pv1 = ud.dists(popupvalue).parameters(2);
            pv2 = ud.dists(popupvalue).parameters(3);
            if cv < pv1 || cv < pv2
                ok = false;
            end
        case {'k','n'}
            pv = ud.dists(popupvalue).parameters(1);
            if cv > pv
                ok = false;
            end
    end
end



%------------------------------------------------------------------------------
% set sliders to use integer or continuous values as appropriate
function setincrement(pslider, intparam)
ss = [0.01 0.1];       % MATLAB default
if (intparam)
   d = max(1, get(pslider,'Max') - get(pslider,'Min'));

   ss = max(1, round(ss * d));
   if (ss(2) <= ss(1)), ss(2) = ss(1) + 1; end
   ss = ss ./ d;
end

set(pslider, 'SliderStep', ss);


%------------------------------------------------------------------------------
% Enable or disable a GUI object
function enableParams(ud, p, state)
if strcmp(state, 'off')
    color =  [0.831373 0.815686 0.784314];
    set(ud.pfhndl(p),'String', '');
    set(ud.hihndl(p),'String', '');
    set(ud.lohndl(p),'String', '');
    set(ud.ptext(p),'String', '');
else
    color = 'white';
end
set(ud.pfhndl(p),'Enable', state, 'BackgroundColor', color);
set(ud.hihndl(p),'Enable', state, 'BackgroundColor', color);
set(ud.lohndl(p),'Enable',state, 'BackgroundColor', color);
set(ud.pslider(p),'Enable', state);
set(ud.ptext(p),'Enable',state);
set(ud.lowerboundtext(p),'Enable',state);
set(ud.upperboundtext(p),'Enable',state);

% End helper functions

% BEGIN CALLBACK FUNCTIONS.

%------------------------------------------------------------------------------
% Callback for changing probability distribution function. CHANGEDISTRIBUTION
function ud = changedistribution(ud)

popupvalue = get(ud.popup,'Value');
parameters = ud.dists(popupvalue).parameters;
phi        = ud.dists(popupvalue).phi;
plo        = ud.dists(popupvalue).plo;
paramnames = ud.dists(popupvalue).paramnames;
intparam   = ud.dists(popupvalue).intparam;

nparams = length(parameters);
enableParams(ud, 1:nparams, 'on');
if nparams < 3
    enableParams(ud, (nparams+1):3, 'off');
end

for idx = 1:nparams
    set(ud.ptext(idx),'String',paramnames{idx});
    set(ud.pfhndl(idx),'String',num2str(parameters(idx)));
    set(ud.lohndl(idx),'String',num2str(plo(idx)));
    set(ud.hihndl(idx),'String',num2str(phi(idx)));
    set(ud.pslider(idx),'Min',plo(idx),'Max',phi(idx),'Value',parameters(idx));
    setincrement(ud.pslider(idx), intparam(idx));
end

xrange = getxdata(popupvalue,ud);
set(gca,'XLim',xrange);

ud=updategui(ud,xrange);


%------------------------------------------------------------------------------
% Callback for controlling lower bound of the parameters using editable text boxes.
function ud = setplo(fieldno,ud)
popupvalue = get(ud.popup,'Value');
intparam = ud.dists(popupvalue).intparam(fieldno);
fieldentry = get(ud.lohndl(fieldno),'String');
cv   = str2double(fieldentry);

pv = str2double(get(ud.pfhndl(fieldno),'String'));
cmax = str2double(get(ud.hihndl(fieldno),'String'));

if intparam
    cv = round(cv);
    set(ud.lohndl(fieldno),'String',num2str(cv));
end
badval = false;

% if the proposed lower limit is larger then the current upper limit, no good
if cv >= cmax
  badval = true;
  
% if the proposed lower limit is smaller than the current upper limit but
% larger then the current value, it must be ok, except for the cross check
% needed for some distributions
elseif cv > pv
  if slidervaluesok(popupvalue, cv, fieldno, ud)
      set(ud.pslider(fieldno),'Min',cv);
      ud.dists(popupvalue).plo(fieldno) = cv;
      set(ud.pfhndl(fieldno),'String',num2str(cv));
      ud = setpfield(fieldno,ud);
  else
      badval = true;
  end
  
% else we need to check if it's larger then the minimum allowed value
else
  pmin = ud.dists(popupvalue).pmin(fieldno);
  pmax = ud.dists(popupvalue).pmax(fieldno);
  if okwithminmax(cv, pmin, pmax, popupvalue, fieldno, ud)
    set(ud.pslider(fieldno),'Min',cv);
    ud.dists(popupvalue).plo(fieldno) = cv;
  else
    badval = true;
  end
end
setincrement(ud.pslider(fieldno), intparam);
xrange = getxdata(popupvalue,ud);
ud = updategui(ud,xrange);

if badval
    preventry = num2str(ud.dists(popupvalue).plo(fieldno));
    wmsg = getString(message('stats:randtool:WarnDlgBadLBValue',...
        fieldentry, preventry));
    uiwait(warndlg(wmsg,...
        getString(message('stats:randtool:WarnDlgTitle')),...
        'modal'))
    set(ud.lohndl(fieldno),'String',preventry);
end


%------------------------------------------------------------------------------
% Callback for controlling upper bound of the parameters using editable text boxes.
function ud = setphi(fieldno,ud)
popupvalue = get(ud.popup,'Value');
intparam = ud.dists(popupvalue).intparam(fieldno);
fieldentry = get(ud.hihndl(fieldno),'String');
cv   = str2double(fieldentry);
pv = str2double(get(ud.pfhndl(fieldno),'String'));
cmin = str2double(get(ud.lohndl(fieldno),'String'));

if intparam
    cv = round(cv);
    set(ud.hihndl(fieldno),'String',num2str(cv));
end

badval = false;

% if the proposed upper limit is samller then the current lower limit, no good
if cv <= cmin
  badval = true;
  
% if the proposed upper limit is larger than the current lower limit but
% smaller then the current value, it must be ok, except for the cross check
% needed for some distributions
elseif cv < pv
  if slidervaluesok(popupvalue, cv, fieldno, ud)
      set(ud.pslider(fieldno),'Max',cv);
      ud.dists(popupvalue).phi(fieldno) = cv;
      set(ud.pfhndl(fieldno),'String',num2str(cv));
      ud = setpfield(fieldno,ud);
  else
      badval = true;
  end
  
% else we need to check if it's smaller then the maximum allowed value
else
  pmin = ud.dists(popupvalue).pmin(fieldno);
  pmax = ud.dists(popupvalue).pmax(fieldno);
  if okwithminmax(cv, pmin, pmax, popupvalue, fieldno, ud)
    set(ud.pslider(fieldno),'Max',cv);
    ud.dists(popupvalue).phi(fieldno) = cv;
  else
    badval = true;
  end
end
setincrement(ud.pslider(fieldno), intparam);
xrange = getxdata(popupvalue,ud);
ud=updategui(ud,xrange);

if badval
    preventry = num2str(ud.dists(popupvalue).phi(fieldno));
    wmsg = getString(message('stats:randtool:WarnDlgBadUBValue',...
        fieldentry, preventry));
    uiwait(warndlg(wmsg,...
        getString(message('stats:randtool:WarnDlgTitle')),...
        'modal'))
    set(ud.hihndl(fieldno),'String',preventry);
end


%------------------------------------------------------------------------------
% Callback for controlling the parameter values using sliders.
function ud = setpslider(sliderno,ud)

% turn off this callback in case we have to put up a warning dialog
cbstr = get(ud.pslider(sliderno),'Callback');
set(ud.pslider(sliderno),'Callback',[]);


try
    popupvalue = get(ud.popup,'Value');
    intparam = ud.dists(popupvalue).intparam(sliderno);

    cv = get(ud.pslider(sliderno),'Value');
    if intparam
        cv = round(cv);
    end

    set(ud.pfhndl(sliderno),'String',num2str(cv));

    if slidervaluesok(popupvalue, cv, sliderno, ud)
        ud.dists(popupvalue).parameters(sliderno) = cv;
    else % handle a conflict between parameters for certain distributions
        pv = ud.dists(popupvalue).parameters(sliderno);
        preventry = num2str(pv);
        wmsg = getString(message('stats:randtool:WarnDlgBadParamValue',...
            num2str(cv), preventry));
        uiwait(warndlg(wmsg,...
            getString(message('stats:randtool:WarnDlgTitle')),...
            'modal'))
        set(ud.pslider(sliderno),'Value',pv);
        set(ud.pfhndl(sliderno),'String',preventry);
    end
    ud=updategui(ud,[]);
catch ME
    set(ud.pslider(sliderno),'Callback',cbstr);
    rethrow(ME)
end

% turn this callback back on after warning dialog dismissed
set(ud.pslider(sliderno),'Callback',cbstr);



%------------------------------------------------------------------------------
% Callback for the parameter value editable text boxes.
function ud = setpfield(fieldno,ud)
popupvalue = get(ud.popup,'Value');
intparam = ud.dists(popupvalue).intparam(fieldno);
fieldentry = get(ud.pfhndl(fieldno),'String');
cv = str2double(fieldentry);
if isnan(cv) || ~isreal(cv)
    goodval = false;
else
    goodval = true;
    pmin = ud.dists(popupvalue).pmin(fieldno);
    pmax = ud.dists(popupvalue).pmax(fieldno);
    phivalue = ud.dists(popupvalue).phi(fieldno);
    plovalue = ud.dists(popupvalue).plo(fieldno);
end

if goodval && intparam
    cv = round(cv);
    set(ud.pfhndl(fieldno),'String',num2str(cv));
end
if goodval && slidervaluesok(popupvalue, cv, fieldno, ud) && ... 
              okwithminmax(cv, pmin, pmax, popupvalue, fieldno, ud)
    set(ud.pslider(fieldno),'Value',cv);
    ud.dists(popupvalue).parameters(fieldno) = cv;
    if (cv >= phivalue), 
        set(ud.hihndl(fieldno),'String',num2str(cv));
        ud = setphi(fieldno,ud); 
        set(ud.pslider(fieldno),'Max',cv);
        setincrement(ud.pslider(fieldno), intparam);
        return; % this return is to avoid using updategui twice.
    end
    if (cv <= plovalue), 
        set(ud.lohndl(fieldno),'String',num2str(cv));
        ud = setplo(fieldno,ud); 
        set(ud.pslider(fieldno),'Min',cv);
        setincrement(ud.pslider(fieldno), intparam);
        return; 
    end
else
    preventry = num2str(ud.dists(popupvalue).parameters(fieldno));
    set(ud.pfhndl(fieldno),'String',preventry);
    wmsg = getString(message('stats:randtool:WarnDlgBadParamValue',...
        fieldentry, preventry));
    uiwait(warndlg(wmsg,...
        getString(message('stats:randtool:WarnDlgTitle')),...
        'modal'))
end
if goodval
   xrange = getxdata(popupvalue,ud);
   ud= updategui(ud, xrange);
end


%------------------------------------------------------------------------------
% Callback to update graphic objects in GUI.
function ud=updategui(ud,xrange)

if isempty(xrange)
   xrange = get(gca,'XLim');
end

popupvalue = get(ud.popup,'Value');

code = ud.dists(popupvalue).code;
if strcmpi(code,'weibull')  % use new name to avoid warning
   code = 'wbl';
end
samples = str2double(get(ud.samples_field,'String'));
if isnan(samples) || samples<=0 || samples~=floor(samples) || ~isreal(samples)
    wmsg = getString(message('stats:randtool:WarnDlgBadSampleSize'));
    uiwait(warndlg(wmsg,...
        getString(message('stats:randtool:WarnDlgTitle')),...
        'modal'))
    samples = 100;
    set(ud.samples_field,'String','100');
end

pval = num2cell(ud.dists(popupvalue).parameters);

set(gcf,'Pointer','watch');
try
   random_numbers = random(code,pval{:},samples,1);
   ud.random_numbers = random_numbers;

   % Create Histogram
   minrn = min(random_numbers);
   maxrn = max(random_numbers);
   if ud.dists(popupvalue).discrete
       edges = minrn-0.5:1:(maxrn+.5); % bins centered on the integers
   else
       bins = floor(sqrt(samples));
       if bins > 1
           edges = linspace(minrn,maxrn+eps(maxrn),bins+1);
       else
           edges = [minrn, minrn+1]; % single bin one unit wide
       end
   end
   counts = histc(random_numbers,edges);
   hbar = bar(edges,counts,'histc');
   set(hbar,'Tag','dbar');

   set(gca,'XLim',xrange, 'XTickMode','auto', 'XTickLabelMode','auto');

   strCounts = getString(message('stats:randtool:yLabelCounts'));
   strValues = getString(message('stats:randtool:xLabelValues'));
   strSamples = getString(message('stats:randtool:LabelSamples'));
   strDistribution = getString(message('stats:randtool:LabelDistribution'));
   
   text(-0.08, 0.45,strCounts,'Units','Normalized', 'Rotation', 90);
   text(0.45,-0.10,strValues,'Units','Normalized');
   text(0.55, 1.06,strSamples,'Units','Normalized');
   text(0.01, 1.06,strDistribution,'Units','Normalized');
   
catch ME
    % Ensure the pointer is set back to arrow.
    set(gcf,'Pointer','arrow');
    m = message('stats:randtool:GeneratorError');
    throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
end
% Set the pointer back to arrow.
set(gcf,'Pointer','arrow');



%------------------------------------------------------------------------------
% Callback to output data
function outrndnum(ud)
popupvalue = get(ud.popup,'Value');
labelnames = get(ud.popup,'String');
def = {ud.dists(popupvalue).rvname};
label = {deblank(labelnames(popupvalue,:))};
item = {ud.random_numbers};
export2wsdlg(label, def, item,getString(message('stats:randtool:TitleExport')));

% END CALLBACK FUNCTIONS.
