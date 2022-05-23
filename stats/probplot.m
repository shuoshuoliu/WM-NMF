function h = probplot(varargin)
%PROBPLOT Probability plot
%   PROBPLOT(Y) produces a normal probability plot comparing the distribution
%   of the data Y to the normal distribution.  Y can be a single vector, or
%   a matrix with a separate sample in each column. The plot includes a
%   reference line that passes through the lower and upper quartiles of
%   Y, and that is useful for judging whether the data follow a normal
%   distribution. PROBPLOT uses midpoint probability plotting positions.
%
%   PROBPLOT('DISTNAME',Y) creates a probability plot for the specified
%   distribution.  DISTNAME is a character string chosen from the following
%   list of distributions:
%
%       'exponential'                      Exponential
%       'extreme value' or 'ev'            Extreme value
%       'half normal' or 'hn'              Half normal
%       'lognormal'                        Lognormal
%       'normal'                           Normal
%       'rayleigh'                         Rayleigh
%       'weibull' or 'wbl'                 Weibull
%
%   PROBPLOT(Y,CENS,FREQ) or PROBPLOT('DISTNAME',Y,CENS,FREQ) requires
%   a vector Y.  CENS is a vector of the same size as Y and contains
%   1 for observations that are right-censored and 0 for observations
%   that are observed exactly.  FREQ is a vector of the same size as Y,
%   containing integer frequencies for the corresponding elements in Y.
%
%   PROBPLOT(AX,Y) takes a handle AX to an existing probability plot, and
%   adds additional lines for the samples in Y.  AX is a handle for a
%   set of axes or uiaxes.
%
%   Use the data cursor to read precise values from the plot, as well as
%   observation numbers from the points when CENS and FREQ are left
%   unspecified.
%
%   PROBPLOT(...,'noref') omits the reference line.
%
%   PROBPLOT(AX,PD) takes a probability distribution object PD, and adds a
%   fitted line to the axes specified by AX to represent the probability
%   distribution specified by PD.
%
%   PROBPLOT(AX,FUN,PARAMS) takes a function FUN and a set of parameters
%   PARAMS, and adds a fitted line to the axes specified by AX.  FUN is
%   a function to compute a cdf function, and is specified with @
%   (such as @wblcdf).  PARAMS is the set of parameters required to
%   evaluate FUN, and is specified as a cell array or vector.  The
%   function must accept a vector of X values as its first argument,
%   then the optional parameters, and must return a vector of cdf
%   values evaluated at X.
%
%   H=PROBPLOT(...) returns handles to the plotted lines.
%
%   The y axis scale is based on the selected probability distribution. The
%   x axis has a log scale for the Weibull and lognormal distributions, and
%   a linear scale for the others.  The Ith sorted value from a sample of
%   size N is plotted against the midpoint in the jump of the empirical CDF
%   on the y axis.  With uncensored data, that midpoint is (I-0.5)/N.  With
%   censored data, the y value is more complicated to compute.
%
%   Example:  Generate exponential data.  A normal probability plot does
%             not show a good fit.  A Weibull probability plot looks
%             better, because the exponential distribution is part of the
%             Weibull family of distributions.
%       y = exprnd(5,200,1);
%       probplot(y);                      % normal probability plot
%       figure; probplot('weibull',y);    % Weibull probability plot
%
%   See also NORMPLOT, WBLPLOT, ECDF.

%   Copyright 2003-2020 The MathWorks, Inc.

% At least one argument required
if nargin > 0
    [varargin{:}] = convertStringsToChars(varargin{:});
end

narginchk(1,Inf);

% Check for a special flag at the end to skip reference lines
if nargin>0 && isequal(varargin{end},'noref')
    narginchk(2,Inf);
    addrefline = false;
    varargin(end) = [];
else
    addrefline = true;
end

% Gather any gpuArray inputs
[varargin{:}] = gather(varargin{:});

% Get probability distribution info from input name, input handle, or default
[hAxes,newdist,dist,param,varargin] = getdistname(varargin);
if length(dist)~=1
    error(message('stats:probplot:InvalidDistribution'));
end

% Get some properties of this distribution
distname = dist.name;
distcode = dist.code;
if dist.uselogpp
    invcdffunc = @(x,varargin) log(feval(dist.invfunc,x,varargin{:}));
    cdffunc = @(x,varargin) feval(dist.cdffunc,exp(x),varargin{:});
else
    invcdffunc = dist.invfunc;
    cdffunc = dist.cdffunc;
end

% Plot either a data sample or a user-supplied fit
if isempty(varargin)
    % Just setting up empty plot
    x = [];
    addfit = false;
elseif isnumeric(varargin{1})
    % Plot the sample or samples in X
    [x,cens,freq,originds] = checkdata(dist,varargin{:});
    addfit = false;
elseif isequal(class(varargin{1}),'function_handle') ...
        || ischar(varargin{1}) ...
        || isa(varargin{1},'ProbDist') || isa(varargin{1},'prob.ProbabilityDistribution')
    % Must be plotting a fit rather than a data sample
    x = [];
    addfit = true;
    if isa(varargin{1},'ProbDist') || isa(varargin{1},'prob.ProbabilityDistribution')
        fitcdffunc = @(x) cdf(varargin{1},x);
        param = {};
    else
        fitcdffunc = varargin{1};
        if length(varargin)>=2
            if iscell(varargin{2})
                param = varargin{2};
            else
                param = num2cell(varargin{2});
            end
        end
    end
else
    error(message('stats:probplot:InvalidData'))
end

% Create probability plot of the proper type
if newdist
    hAxes = createprobplot(hAxes,distcode,distname,param,dist,...
                                 cdffunc,invcdffunc);
end

isUIAxes = false;
if isa(hAxes,'matlab.ui.control.UIAxes')
    isUIAxes = true;
    % Need to set hold on on UIAxes because each call to fplot overwrites
    % the previous one. Since these axes are always passed in (never gca),
    % the syntax definition means that we should not clear them.
    isHoldOn = ishold(hAxes);
    hold(hAxes,'on');
end

% Plot a data sample
hfit = [];
href = [];
hdata = [];
if ~isempty(x)
    % Get plotting positions and expanded x data vector for plot
    [xData,ppos] = getppos(x,cens,freq);

    % Draw fit lines first, so they will be below the data points
    if addrefline
        % Add lines representing fits
        nsamples = size(xData,2);
        samplenum = zeros(nsamples,1);
        for j=1:nsamples
           % Get some points on the reference line
           if isempty(cens) && isempty(freq)
              probpoints = [];
           else
              probpoints = ppos;
           end
           
           refxy = getrefxy(xData(:,j),probpoints,...
                            invcdffunc,param,dist.uselogpp,distname);

           if ~isempty(refxy)
              % Create a function that will put a line through these points
              %      A = y1 - x1*(y2-y1)/(x2-x1)
              %      B = (y2-y1)/(x2-x1)
              B = (refxy(2,2)-refxy(1,2)) ./ (refxy(2,1)-refxy(1,1));
              A = refxy(1,2) - refxy(1,1) .* B;
              if dist.uselogpp
                 % The reference x1 and x2 values are already on the log scale
                 linefun = @(x)(A+B*log(x));
              else
                 linefun = @(x)(A+B*x);
              end

              % Create a reference line based on this function
              if isUIAxes
                  % primitve.FunctionLine is not defined for uiaxes, call
                  % fplot instead
                  hrefj = fplot(hAxes,linefun,[min(xData(:,j)),max(xData(:,j))]);
              else
                hrefj = matlab.graphics.chart.primitive.FunctionLine(...
                    'Function',linefun,'Parent',hAxes);
              end
              href = [href; hrefj];
              samplenum(numel(href)) = j;
              
              % Add custom data cursor
              hB = hggetbehavior(hrefj,'datacursor');
              set(hB,'UpdateFcn',{@probplotDatatipCallback,hAxes});
              if nsamples>1
                  setappdata(hrefj,'group',j);
              end

              % Start off with reasonable default properties
              set(hrefj,'XLimInclude','off','YLimInclude','off',...
                          'Color','k','LineStyle','--','Tag','reference');
           end
        end
    end
    
    % Plot data points
    hdata = adddatapoints(hAxes,xData,ppos,originds);

    if addrefline && nsamples>1
        % With multiple samples, make sure fit color matches data color
        for j=1:numel(href)
            set(href(j),'Color',get(hdata(samplenum(j)),'Color'));
        end
    end
end

% Add a fit
if addfit
    hfit = ppaddfit(hAxes,fitcdffunc,param);
end

% Now that y limits are established, define good y tick positions and labels
setprobticks(hAxes,distname,hdata);

% Reset hold state if UIAxes
if isUIAxes && ~isHoldOn
    hold(hAxes,'off');
end

if nargout>0
   h = [hdata(:); href(:); hfit(:)];
end

% -------------------------------------
function hAxes = createprobplot(hAxes,dcode,dname,param,dist,cdffunc,invcdffunc)
%CREATEPROBPLOT Create an empty probability plot of a particular type

% Create a set of axes if none supplied
if isempty(hAxes)
   hAxes = cla('reset');
end

% Store information about the reference distribution
setappdata(hAxes,'ReferenceDistribution',dcode);
setappdata(hAxes,'CdfFunction',cdffunc);
setappdata(hAxes,'InverseCdfFunction',invcdffunc);
setappdata(hAxes,'DistributionParameters',param);
setappdata(hAxes,'LogScale',dist.uselogpp);
setappdata(hAxes,'DistSpec',dist);

% Add labels and a title
set(get(hAxes,'XLabel'),'String',getString(message('stats:probplot:Data')));
set(get(hAxes,'YLabel'),'String',getString(message('stats:probplot:Probability')));
if isempty(param)
    paramtext = '';
else
    fmt = repmat('%g,',1,length(param));
    fmt(end) = [];
    paramtext = sprintf(['(' fmt ')'],param{:});
end
set(get(hAxes,'Title'),'String',...
                     sprintf('%s',getString(message('stats:probplot:ProbPlotForDistribution',...
                             dname,paramtext))));

% Set X axis to log scale if this distribution uses that scale
if dist.uselogpp
   set(hAxes,'XScale','log');
end

% ------------------------------------
function hdata = adddatapoints(hAxes,x,ppos,originds)
%ADDDATAPOINTS Add points representing data to a probability plot
%   x is already sorted

% Get information about the reference distribution
invcdffunc = getappdata(hAxes,'InverseCdfFunction');
param = getappdata(hAxes,'DistributionParameters');

% Compute y values for this distribution
q = feval(invcdffunc,ppos,param{:});

% Add to plot
hdata = line(x,q,'Parent',hAxes);

% Attach custom data cursor
for i=1:length(hdata)
    hB = hggetbehavior(hdata(i),'datacursor');
    set(hB,'UpdateFcn',{@probplotDatatipCallback,hAxes});
    if ~isempty(originds)
        setappdata(hdata(i),'originds',originds(:,i));
    end
    if numel(hdata)>1
        setappdata(hdata(i),'group',i);
    end
end

% Use log X scale if distribution requires it
xmin = min(x(1,:));
xmax = max(x(end,:));
if isequal(get(hAxes,'XScale'),'log')
    xmin = log(xmin);
    xmax = log(xmax);
    dx = 0.02 * (xmax - xmin);
    if (dx==0)
       dx = 1;
    end
    xmin = exp(xmin - dx);
    xmax = exp(xmax + dx);
    newxlim = [xmin xmax];
else
    dx = 0.02 * (xmax - xmin);
    if (dx==0)
       dx = 1;
    end
    newxlim = [xmin-dx, xmax+dx];
end
oldxlim = get(hAxes,'XLim');
set(hAxes,'XLim',[min(newxlim(1),oldxlim(1)), max(newxlim(2),oldxlim(2))]);

% Make sure they have different markers and no connecting line
markerlist = {'x','o','+','*','s','d','v','^','<','>','p','h'}';
set(hdata,'LineStyle','none','Tag','data',...
          {'Marker'},markerlist(1+mod((1:length(hdata))-1,length(markerlist))));

% -------------------------------------------
function setprobticks(hAxes,distname,hndl)
%SETPROBTICKS Set the y axis tick marks to use a probability scale

invcdffunc = getappdata(hAxes,'InverseCdfFunction');
param = getappdata(hAxes,'DistributionParameters');

% Define tick locations
if strcmp(distname,'Half Normal')
    % Make small probabilities not bunched up
    ticklevels = [.1 .25 .5 .75 .9 .95 .99 .995 .999 .9995 .9999];
else
    ticklevels = [.0001 .0005 .001 .005 .01 .05 .1 .25 .5 ...
              .75 .9 .95 .99 .995 .999 .9995 .9999];
end
tickloc = feval(invcdffunc,ticklevels,param{:});

% Place probability-scale ticks at the right locations, making sure to use
% blank labels in case of overlap between adjacent ticks.
placeTicks(hAxes,tickloc,ticklevels,hndl);

% Repeat this operation whenever the axes size changes or the limits change
addlistener(hAxes,'SizeChanged',@(varargin)placeTicks(hAxes,tickloc,ticklevels,hndl));
try   % if ~isa(hAxes,'matlab.ui.control.UIAxes')
    addlistener(hAxes,'YLim','PostSet',@(varargin)placeTicks(hAxes,tickloc,ticklevels,hndl));
catch
end

function placeTicks(hAxes,tickloc,ticklevels,hndl)
% Remove ticks outside axes range
if ~ishandle(hAxes)
    % Axes are gone, this is no longer relevant
    return
end
if ~any(ishandle(hndl),'all')
    % Probplot contents are gone, possibly after newplot, this is no longer
    % relevant
    return
end
ylim = get(hAxes,'YLim');
t = tickloc>=ylim(1) | tickloc<=ylim(2);
tickloc = tickloc(t);
ticklevels = ticklevels(t);

% Remove ticks that are too close together
ttt = text(hAxes.XLim(1),hAxes.YLim(1),'M','Parent',hAxes,'Visible','off','Units','data');
yheight = ttt.Extent(end);
delete(ttt);

temploc = tickloc;
prevtick = temploc(1);
lasttick = temploc(end);
for j=2:floor(length(temploc)/2)
    % March through the ticks from each end. Remove the label from any tick
    % that would overlap its non-removed neighbor.
    if temploc(j)-prevtick < yheight
        temploc(j) = NaN;
    else
        prevtick = temploc(j);
    end
    if lasttick-temploc(end+1-j) < yheight
        temploc(end+1-j) = NaN;
    else
        lasttick = temploc(end+1-j);
    end
end

remove = isnan(temploc);
if ~isequal(hAxes.YTick,tickloc)
    % Apparently we are doing this for the first time, so set both tick
    % positions and labels. Update the labels afterward, after they are
    % converted to text.
    doupdate = true;
else
    % Tick locations are not changed. Figure out if we need to update the
    % ticks with new labels having different blank values.
    oldlevels = hAxes.YTickLabel;
    wasnan = all(oldlevels==' ',2);
    doupdate = ~isequal(remove,wasnan);
end
if doupdate
    set(hAxes,'YTick',tickloc,'YTickLabel',ticklevels);
    hAxes.YTickLabel(remove,:) = ' ';
end

% -------------------------------------------
function [x,ppos] = getppos(x,cens,freq)
%GETPPOS Get plotting positions for probability plot

if isempty(cens) && isempty(freq)
   % Use a simple definition compatible with the censored calculation.
   % Plot sorted x(i) against (i-1/2)/n.
   % This code supports X as a matrix or a vector.
   if ~any(isnan(x(:)))
       n = size(x,1);
       ppos = ((1:n)' - 0.5) / n;
   else
       % Any NaNs in X were sorted to the bottom, so we will fill in
       % ppos values starting at the top
       ppos = nan(size(x));
       for j=1:size(x,2)
           n = sum(~isnan(x(:,j)));
           ppos(1:n,j) = ((1:n)' - 0.5) / n;
       end
   end
else
   % Compute the empirical cdf
   [fecdf,xecdf,temp1,temp2,D] = ecdf(x,'cens',cens,'freq',freq); %#ok<ASGLU>
   
   % Create outputs with one row for each observed failure
   N = sum(min(100,D));
   xout = zeros(N,1);
   ppos = zeros(N,1);
   
   % Fill in with points equally spaced at each jump of the cdf
   % If there are M failures at x(i), plot x(i) against the M values
   % equally spaced between F(x(i)-) and F(x(i)+)
   xbase = 0;
   for j=1:length(xecdf)-1
      Dj = D(j);
      Npts = min(100,Dj);
      rownums = xbase+(1:Npts);
      xout(rownums) = xecdf(j+1);
      ppos(rownums) = fecdf(j) + ((1:Npts)-0.5) * (fecdf(j+1)-fecdf(j)) / Npts;
      xbase = xbase + Npts;
   end
   
   % Replace old data with the new version
   x = xout;
end

% -------------------------------------------
function refxy = getrefxy(x,probpoints,invcdffunc,param,uselogpp,distname)
%GETREFXY Get two points on a reference line

if isempty(probpoints)
    % If there is no censoring, we get the first and third quartiles. For
    % half normal plot, we uses zero and median.
    if strcmp(distname,'Half Normal')
        pt1 = 0;
        pt2 = .5;
    else
        pt1 = .25;
        pt2 = .75;
    end
    DataQ = prctile(x,100*[pt1; pt2]);
    DataQ1 = DataQ(1);
    DataQ3 = DataQ(2);
else
    % Get upper quartile, or as high up as we can go
    pt2 = min(.75,probpoints(end));
    DataQ3 = interp1(probpoints,x,pt2);
    
    % Get lower quartile or a substitute for it
    pt1 = max(pt2/3,probpoints(1));
    DataQ1 = interp1(probpoints,x,pt1);
end

% Use log scale if necessary
if uselogpp
    DataQ1 = log(DataQ1);
    DataQ3 = log(DataQ3);
end

% Get the y values for these x values
DistrQ = feval(invcdffunc, [pt1 pt2], param{:});

% Package up the points and return them
if DataQ3 > DataQ1 % make sure we have distinct points
    refxy = [DataQ1 DistrQ(1)
        DataQ3 DistrQ(2)];
else
    refxy = [];
end

% ------------------------------------------
function hfit = ppaddfit(hAxes,cdffunc,params)
%PPADDFIT Add fit to probability plot

if nargin<3 || isempty(params)
   params = {};
elseif isnumeric(params)
   params = num2cell(params);
end

% Define function using local function handle
isUIAxes = isa(hAxes,'matlab.ui.control.UIAxes');

if isequal(get(hAxes,'XScale'),'log')
    xlim = get(hAxes,'XLim');
    xtemp = logspace(log10(xlim(1)),log10(xlim(2)),300)';
    hfit = line(xtemp,calcfit(xtemp,hAxes,cdffunc,params),'Parent',hAxes);
elseif isUIAxes
    hfit = fplot(hAxes,@(x)calcfit(x,hAxes,cdffunc,params));
else
    hfit = matlab.graphics.chart.primitive.FunctionLine(...
             'Function',@calcfit,'UserArgs',{hAxes,cdffunc,params},...
             'Parent',hAxes);
end

% Add custom data cursor
hB = hggetbehavior(hfit,'datacursor');
set(hB,'UpdateFcn',{@probplotDatatipCallback,hAxes});

% Start off with reasonable default properties
set(hfit,'XLimInclude','off','YLimInclude','off',...
         'Color','k','LineStyle','--','Tag','fit');

% ------------------------------------------
function fx = calcfit(x,hAxes,cdffunc,fitparams)
%CALCFIT Calculated values of function for plot

% Get y values in units of the tick mark labels
p = feval(cdffunc,x,fitparams{:});

% Get y values in units of the real y axis
invcdffunc = getappdata(hAxes,'InverseCdfFunction');
plotparams = getappdata(hAxes,'DistributionParameters');
fx = feval(invcdffunc,p,plotparams{:});

% ------------------------------------------
function [hAxes,newdist,dist,param,vin] = getdistname(vin)
%GETDISTNAME Get probability distribution info from user input

% Check for probability plot axes handle as a first argument
if numel(vin{1})==1 && ishghandle(vin{1})
   hAxes = vin{1};
   if ~isequal(get(hAxes,'Type'),'axes')
      error(message('stats:probplot:BadHandle'));
   end
   vin(1) = [];
else
   hAxes = [];

   % Try to respect "hold on" if the target is a probability plot
   H = get(groot,'CurrentFigure');
   if ~isempty(H) && ishghandle(H) && strcmp(H.Type,'figure')
       H = H.CurrentAxes;
   end
   if ~isempty(H) && ishold(H)
        dname = getappdata(H,'ReferenceDistribution');
        param = getappdata(H,'DistributionParameters');
        isprobplot = ~isempty(getappdata(H,'InverseCdfFunction')) && iscell(param);
        if isprobplot
            hAxes = H;
        end
   end
end

% Now get the distribution name and parameters for it
if ~isempty(vin) && (isa(vin{1},'ProbDist') || isa(vin{1},'prob.ProbabilityDistribution')) && isempty(hAxes)
   % ProbDist object followed by data
   pd = vin{1};
   
   try
       dist = dfgetdistributions(pd.DistName);
   catch
       if isa(pd,'prob.TriangularDistribution') || isa(pd,'prob.PiecewiseLinearDistribution') || isa(pd,'prob.MultinomialDistribution')
           % These three distributions are not supported and they do not have 
           % property 'DistName'; if they are passed in, error out the same 
           % message as using probplot('DISTNAME',Y)
           error(message('stats:probplot:UnrecognizedDistribution', pd.DistributionName));
       elseif  isa(pd,'prob.UniformDistribution')
           % Uniform distribution does not have property 'DistName'
           % neither, but it differs from the three distributions above:
           % once the parameter estimates are provided, probplot works for
           % Uniform distirbution
           dist = dfgetdistributions('uniform');
       end
   end
   
   newdist = true;
   vin(1) = [];
   param = num2cell(pd.ParameterValues);
elseif ~isempty(vin) && (ischar(vin{1}) || iscell(vin{1}))
   % Preference is to use passed-in name, if any
   dname = vin{1};
   vin(1) = [];
   newdist = true;
   
   % May be passed in as a name or a {dist,param} array
   if iscell(dname)
      dist = dname{1};
      param = num2cell(dname{2});
   elseif ischar(dname)
      if strncmpi(dname,'ev',length(dname))
          dname = 'extreme value';
      elseif strncmpi(dname,'wbl',length(dname))
          dname = 'weibull';
      end
      dist = dfgetdistributions(dname);
      if ~isscalar(dist)
         error(message('stats:probplot:UnrecognizedDistribution', dname));
      elseif ~dist.islocscale
         error(message('stats:probplot:InappropriateDistribution', dname));
      end
      param = {};
   else
      error(message('stats:probplot:BadDistribution'));
   end   
else
    if ~isempty(hAxes)
        dname = getappdata(hAxes,'ReferenceDistribution');
        param = getappdata(hAxes,'DistributionParameters');
        notpp = isempty(getappdata(hAxes,'InverseCdfFunction')) || ~iscell(param);
    end
    if isempty(hAxes) || (notpp && ~ishold(hAxes))
        % Use default if no axes passed in
        dname = 'normal';
        param = {};
        newdist = true;
    else
        % Otherwise use current distribution
        if notpp
            error(message('stats:probplot:NotProbPlotHandle'))
        end
        newdist = false;
    end
    dist = dfgetdistributions(dname);
end

% ------------------------------------------
function [x,cens,freq,originds] = checkdata(dist,varargin)
%CHECKDATA Get data and check that it works with this distribution

x = varargin{1};
if length(varargin)<2
   cens = [];
else
   cens = varargin{2};
end
if length(varargin)<3
   freq = [];
else
   freq = varargin{3};
end

if ~isempty(cens) || ~isempty(freq)
    % Remove NaNs now if we have to maintain x, cens, and freq in parallel.
    % Otherwise if we don't have cens and freq, we'll deal with NaNs in x
    % as required.  X must be a vector in this case.
    [~,~,x,cens,freq] = internal.stats.removenan(x,cens,freq);
end

% Follow the usual convention by treating a row vector as a column
if ~ismatrix(x)
   error(message('stats:probplot:YVectorMatrix'))
end
if size(x,1)==1
   x = x';
end
[x,sortidx] = sort(x);

[nobs,nsamples] = size(x);
if ~isempty(cens)
   if length(cens)~=nobs || ~(isnumeric(cens) || islogical(cens))
      error(message('stats:probplot:BadCensSize'))
   end
   cens = cens(sortidx);
end
if ~isempty(freq)
   if length(freq)~=nobs || ~(isnumeric(freq) || islogical(freq))
      error(message('stats:probplot:BadFreqSize'))
   end
   freq = freq(sortidx);
end
if isempty(freq) && isempty(cens)
    originds = sortidx;
else
    originds = [];
end

% Check match between data and distribution.
xmin = min(x(1,:));
xmax = max(x(end,:));
if xmin==dist.support(1) && ~dist.closedbound(1)
    error(message('stats:probplot:OpenLowerSupport', dist.name, sprintf( '%g', dist.support( 1 ) )));
elseif xmin<dist.support(1)
    error(message('stats:probplot:LowerSupport', dist.name, sprintf( '%g', dist.support( 1 ) )));
elseif xmax==dist.support(2) && ~dist.closedbound(2)
    error(message('stats:probplot:OpenUpperSupport', dist.name, sprintf( '%g', dist.support( 2 ) )));
elseif xmax>dist.support(2)
    error(message('stats:probplot:UpperSupport', dist.name, sprintf( '%g', dist.support( 2 ) )));
elseif (nsamples>1) && ~dist.islocscale
    error(message('stats:probplot:SingleSampleDistribution', dist.name));
elseif (~isempty(cens) || ~isempty(freq)) && (nsamples>1)
    error(message('stats:probplot:SingleSample'));
end


% ------------------------------------------
function datatipTxt = probplotDatatipCallback(~,evt,hAxes)

target = get(evt,'Target');
ind = get(evt,'DataIndex');
pos = get(evt,'Position');

x = pos(1);
y = pos(2);

cdffunc = getappdata(hAxes,'CdfFunction');
param = getappdata(hAxes,'DistributionParameters');

% Compute position to display.
yper = cdffunc(y,param{:});

% Get the original row number of the selected point, which is set if 
% freq and cens are left unspecified. Also empty for reference lines.
originds = getappdata(target,'originds');

% Get the group number, which is set on points and reference lines 
% if more than one series.
group = getappdata(target,'group');

% Generate text to display.
datatipTxt = {
    sprintf('%s: %s',getString(message('stats:probplot:Data')),num2str(x)),...
    sprintf('%s: %s',getString(message('stats:probplot:Probability')),num2str(yper)),...
    };
if ~isempty(originds) || ~isempty(group)
    datatipTxt{end+1} = '';
end
if ~isempty(originds)
    origind = originds(ind);
    datatipTxt{end+1} = sprintf('%s: %s',...
        getString(message('stats:probplot:Observation')),num2str(origind));
end
if ~isempty(group)
    datatipTxt{end+1} = sprintf('%s: %s',...
        getString(message('stats:probplot:Group')),num2str(group));
end
%datatipTxt = [sprintf('%s\n',datatipTxt{1:end-1}), datatipTxt{end}];
