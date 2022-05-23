function setfunction(dffig,ftype,dtype)
%SETFUNCTION Set value of probability function to be plotted


%   Copyright 2003-2020 The MathWorks, Inc.

% Figure out what we're changing from
oldftype = dfgetset('ftype');
olddtype = dfgetset('dtype');

% Remember these settings
dfgetset('ftype',ftype);
if nargin>=3
   dfgetset('dtype',dtype);
else
   dtype = [];
end

% Remember current x axis limits
ax = findobj(dffig, 'Tag', 'dfMainAxes');
oldxlim = get(ax,'XLim');
oldxscale = get(ax,'XScale');

% Remember legend properties
legh = legend('-find',ax);
leginfo = stats.internal.dfit.getlegendinfo(legh);

% Get rid of axis controls
if isequal(dfgetset('showaxlimctrl'),'on')
   stats.internal.dfit.toggleaxlimctrl(dffig);
end

% Determine and remember if this function type supports bounds
oktypes = {'cdf' 'survivor' 'cumhazard' 'icdf'};
dobounds = ismember(ftype, oktypes);
dfgetset('dobounds',dobounds);

% Get the array of data sets
dsdb = stats.internal.dfit.getdsdb;
ds = down(dsdb);

if ~isequal(oldftype,ftype)
   % Change the function type for each one
   while(~isempty(ds))
      setftype(ds,ftype);
      ds = right(ds);
   end
elseif isequal(ftype,'probplot') && ~isequal(olddtype,dtype)
   while(~isempty(ds))
      clearplot(ds);
      ds = right(ds);
   end
end   

% Get the array of fits
fitdb = stats.internal.dfit.getfitdb;
ft = down(fitdb);
referenceFit = [];

if ~isequal(oldftype,ftype)
   % Change the function type for each one
   while(~isempty(ft))
      setftype(ft,ftype);
      ft = right(ft);
   end
end

% Determine if a specific set of parameters (reference fit) is required
if isequal(ftype,'probplot')
   if ishandle(dtype)
      referenceFit = dtype;
   else
      referenceFit = [];
   end
end

% Update title, labels, appdata, and (for probability plots) axes
newplot(ax);
setappdata(ax,'ReferenceDistribution','');
setappdata(ax,'CdfFunction','');
setappdata(ax,'InverseCdfFunction','');
setappdata(ax,'DistributionParameters','');
setappdata(ax,'LogScale','');

% Define the colors to be used here
a = [3 0 2 1 3 3 3 2 2 0 2 3 0 1 2 1 0 1 0 1 1
     0 0 1 1 0 3 2 2 1 2 0 1 3 2 3 0 1 3 0 2 0
     0 3 0 1 3 0 1 2 3 1 1 2 0 3 1 2 2 2 0 0 2]'/3;
set(ax,'ColorOrder',a);

% Turn on grid if requested
if isequal(dfgetset('showgrid'), 'on')
   set(ax,'xgrid','on','ygrid','on')
end

if isequal(ftype,'probplot')
   if isempty(referenceFit)
      probplot(ax,dtype);
   else
      probplot(ax,{referenceFit.distspec, referenceFit.params})
   end
   title(ax,'');
elseif isequal(ftype,'icdf')
   set(get(ax,'XLabel'),'String',getString(message('stats:dfstrings:errormsg_Probability')));
   set(get(ax,'YLabel'),'String',getString(message('stats:dfstrings:label_Quantile')));
else
   othertypes =  {'pdf'               'cdf'                    ...
                  'survivor'          'cumhazard'};
   otherlabels = {getString(message('stats:dfstrings:errormsg_Density')) ...
                  getString(message('stats:dfstrings:label_CumulativeProbability')) ...
                  getString(message('stats:dfstrings:dropdown_SurvivorFunction')) ...
                  getString(message('stats:dfstrings:dropdown_CumulativeHazard'))};
   jtype = find(strcmpi(ftype,othertypes));
   if isempty(jtype)   % should never happen
      ylab = ftype;
   else
      ylab = otherlabels{jtype};
   end
   set(get(ax,'XLabel'),'String',getString(message('stats:dfstrings:label_Data')));
   set(get(ax,'YLabel'),'String',ylab);
end
set(ax, 'box','on','Tag','dfMainAxes',...
        'XLimMode','manual','YLimMode','manual','ZLimMode','manual',...
        'CLimMode','manual','AlimMode','manual');

samexaxis = {'cdf' 'survivor' 'cumhazard' 'probplot'};

% Reset the x limits, update plotted curves, and set y limits
stats.internal.dfit.updateallplots(true,false);    % update data sets
newxscale = get(ax,'XScale');
if ismember(ftype,samexaxis) && isequal(oldxscale,newxscale) && ...
     (ismember(oldftype,samexaxis) || isequal(oldftype,'pdf'))
   stats.internal.dfit.updatexlim(oldxlim,false);     % use old x limits
else
   stats.internal.dfit.updatexlim([],false);          % get new x limits
end
stats.internal.dfit.updateallplots(false,true);    % update fits
stats.internal.dfit.updateylim;                    % now compute y limits

% Update the legend
stats.internal.dfit.updatelegend(dffig,false,leginfo);

% Update the Data Sets Table
dataObj = stats.internal.dfit.Data.getInstance();
dataObj.updateManageDataSetsTable(); 
