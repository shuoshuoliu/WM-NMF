function h = initdistfitApp(h,name,fitEditor)
%INITDISTFIT Helper function to allow object construction outside the constructor

% Copyright 2003-2020 The MathWorks, Inc.

% This is the meat of the constructor.  It is here, rather than in distfit.m, 
% because there is no way to have another method, specifically copyfit, call
% the constructor.  Since it is a method of the distfit object, it only
% calls the builtin.  An example of this workaround in action is in copyfit.

% Initialize data members
h.support = 'unbounded';
h.numevalpoints = 100;
h.x = [];
h.y = [];
h.ftype = 'cdf';
h.conflev = 0.95;
h.exclusionrulename= '';
h.plot=0;
h.linehandle=[];
h.ColorMarkerLine = [];


% Create a default fit name if none was supplied
if nargin==1
	name = stats.internal.dfit.getfitname();
end

h.name=name;

if nargin<3
   h.fitframe = [];
else
   h.fitframe=fitEditor; % an instance of the Fit.m class
end

list(9) = handle.listener(h,findprop(h,'kernel'),'PropertyPostSet',...
                           {@handleTableInfoChanged});
list(8) = handle.listener(h,findprop(h,'distname'),'PropertyPostSet',...
                           {@handleTableInfoChanged});
list(7) = handle.listener(h,findprop(h,'dataset'),'PropertyPostSet',...
                          {@handleTableInfoChanged});
list(6) = handle.listener(h,findprop(h,'name'),'PropertyPostSet',...
                          {@updatename,h});
list(5) = handle.listener(h,findprop(h,'conflev'),'PropertyPostSet',...
                          {@updateconflev,h});
list(4) = handle.listener(h,findprop(h,'showbounds'),'PropertyPostSet', {@changebounds,h});
list(3) = handle.listener(h,findprop(h,'plot'),'PropertyPostSet', ...
                          {@localupdate,h});
list(2) = handle.listener(h,'ObjectBeingDestroyed', {@cleanup,h});
list(1) = handle.listener(h,findprop(h,'dshandle'),'PropertyPostSet',...
                          {@updatelim,h});

h.listeners=list;
dfgetset('dirty',true);   % session has changed since last save

%=============================================================================
function handleTableInfoChanged(~, ~)
dfgetset('dirty',true);   % session has changed since last save

%=============================================================================
function updatename(~,~,fit)

if fit.plot && ishghandle(fit.linehandle)
   stats.internal.dfit.updatelegend(fit.linehandle);
end
dfgetset('dirty',true);   % session has changed since last save

%=============================================================================
function updateconflev(~,~,fit)

if ~isempty(fit.ybounds) && isequal(fit.fittype,'param')
   updateplot(fit);
end
%=============================================================================
function updatelim(~,~,fit)

fit.xlim = xlim(fit);

%=============================================================================
function localupdate(~,~,fit)

% Update plotted curve
updateplot(fit);

% Update plot limits
stats.internal.dfit.updatexlim;
stats.internal.dfit.updateylim;

%Update Fits table to show current state
fitManagerObj = stats.internal.dfit.FitsManager.getInstance();
fitManagerObj.updateFitsTable();
       
dfgetset('dirty',true);   % session has changed since last save

%=============================================================================
function changebounds(~,~,fit)

dist = fit.distspec;
if isempty(dist)
   hasbounds = false;
else
   hasbounds = dist.hasconfbounds;
end

if fit.showbounds && ~hasbounds
   % Undo requests for bounds if this distribution can't handle that
   fit.showbounds = hasbounds;
else
   % Otherwise the change is successful, so update the plot
   updateplot(fit);
   dfgetset('dirty',true);   % session has changed since last save
end


%=============================================================================
function cleanup(~,~,fit)

if ~isempty(fit.linehandle)
   if ishghandle(fit.linehandle)
      list = fit.listeners;
      list(3).enable = 'off';
      fit.plot = 0;
      updateplot(fit);
      % Update plot limits
     stats.internal.dfit.updatexlim;
     stats.internal.dfit.updateylim;
   end
end

dfgetset('dirty',true);   % session has changed since last save

% Update list of probability plot distributions to remove this fit
stats.internal.dfit.updateppdists;
