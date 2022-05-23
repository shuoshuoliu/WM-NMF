function varargout=noJavaDFApp(varargin)
%distributionFitter Open Distribution Fitter non-Java version
%
%   distributionFitter opens an app for fitting distributions to data.
%   This app enables you to import data from the workspace, fit a variety
%   of distributions, and use plots to compare the fitted distributions to
%   each other and to the data.
%
%   distributionFitter(Y) displays the Distribution Fitter App and
%   creates a data set with data specified by the vector y.
%
%   distributionFitter(Y,CENS) uses the vector CENS to specify whether the
%   observation Y(j) is censored, (CENS(j)==1) or observed exactly
%   (CENS(j)==0).  If CENS is omitted or empty, no Y values are censored.
%
%   distributionFitter(Y,CENS,FREQ) uses the vector FREQ to specify the
%   frequency of each element of Y.  If FREQ is omitted or empty, all Y
%   values have a frequency of 1.
%
%   distributionFitter(Y,CENS,FREQ,'DSNAME') creates a data set with the
%   name 'dsname' using the data vector Y, censoring indicator CENS, and
%   frequency vector FREQ.
%
%   See also MLE, DISTTOOL, RANDTOOL.

%   Copyright 2019-2020 The MathWorks, Inc.

if nargin > 0
    [varargin{:}] = convertStringsToChars(varargin{:});
end

% Handle call-back
if (nargin > 0 && ischar(varargin{1}))
    out = switchyard(varargin{:});
    if nargout>0
       varargout = {out};
    end
    return
end

% Get a reference to the Distribution Fitting App
dft = stats.internal.dfit.DistributionFitting.getInstance;

% make sure there are instances of the udd databases for datasets, fits and outliers 
stats.internal.dfit.getdsdb;
stats.internal.dfit.getfitdb;
stats.internal.dfit.getoutlierdb;

% Send the gui information about the distributions we can fit
dfgetset('alldistributions','');  % clear definitions left-over from before
[dists,emsg] = dfgetdistributions;
if ~isempty(emsg)
   edlg = errordlg(getString(message('stats:dfstrings:dlg_ErrorImportCustom',...
                                     emsg)),...
                   getString(message('stats:dfstrings:dlg_DistributionFittingTool')));
else
   edlg = [];
end
stats.internal.dfit.setdistributions(dft,dists);

dfgetset('dft',dft);

% Try to get old figure
dffig = dfgetset('dffig');

% If the handle is empty, create the object and save the handle
makefig = (isempty(dffig) || ~ishghandle(dffig));
if makefig
   mainPanel = stats.internal.dfit.MainPanel();
   dffig = mainPanel.DistributionFitterUIFigure;
   stats.internal.dfit.session('clear');
   stats.internal.dfit.setfunction(dffig,'pdf');
   
   % Initialize default bin width rules information
   initdefaultbinwidthrules;
   
   dfgetset('dirty',true);   % session has changed since last save
else
   figure(dffig);
end
 
% Start with input data, or put up message about importing data
ds = [];
if nargin>0
   % If data were passed in, set up argument list for that case
   dsargs = {[] [] [] ''};
   n = min(4,nargin);
   dsargs(1:n) = varargin(1:n);
   for j=1:min(3,n)
      dsargs{4+j} = inputname(j);   % get data names if possible
   end
   [ds,err] = stats.internal.dfit.createdataset(dsargs{:});
   
   if ~isempty(err)
      err = getString(message('stats:dfstrings:sprintf_ErrorImportData',err));
      errordlg(err,getString(message('stats:dfstrings:dlg_BadInputData')),'modal');
   end
   dfgetset('dirty',true);   % session has changed since last save
   if ~makefig
      delete(findall(dffig,'Tag','dfstarthint')); % remove any old message
   end
end

ax = findobj(dffig, 'Tag', 'dfMainAxes');
if makefig && (nargin==0 || isempty(ds))
   text(ax, .5, .5, getString(message('stats:dfstrings:gui_StartUpMessage')),...
        'Tag','dfstarthint',...
        'HorizontalAlignment','center', ...
        'Units', 'normalized');
end

if nargout==1
   varargout={dft};
end

if ~isempty(edlg)
    figure(edlg);
end

% --------------------------------------------
function out = switchyard(action,varargin)
%SWITCHYARD Dispatch menu call-backs and other actions to private functions

dffig = dfgetset('dffig');
out = [];

switch(action)
    % Fitting actions
    case 'addsmoothfit'
         fit = stats.internal.dfit.addsmoothfit(varargin{:});
         if ~isempty(fit)
              setPanelName(fit.fitframe, getString(message('stats:dfittool:title_editFit')));
         end
         
         stats.internal.dfit.updatelegend(dffig);
         stats.internal.dfit.updateylim;
         if ~isempty(fit)
             out = fit;
         end
         stats.internal.dfit.updateppdists(dffig);
    case 'addparamfit'
         fit = stats.internal.dfit.addparamfit(varargin{:});
         if ~isempty(fit)
            setPanelName(fit.fitframe, getString(message('stats:dfittool:title_editFit')));
         end

         stats.internal.dfit.updatelegend(dffig);
         stats.internal.dfit.updateylim;
         if ~isempty(fit)
             out = fit;
          end
         stats.internal.dfit.updateppdists(dffig);
    % Various graph manipulation actions
    case 'adjustlayout'
         dfgetset('oldposition',get(dffig,'Position'));
         dfgetset('oldunits',get(dffig,'Units'));
    case 'defaultaxes'
         stats.internal.dfit.updatexlim([],true,true);
         stats.internal.dfit.updateylim(true);
end

% --------------------------------------------
function initdefaultbinwidthrules()
binDlgInfo = struct('rule', 1, 'nbinsExpr', '', 'nbins', [], 'widthExpr', '', ...
                    'width', [], 'placementRule', 1, 'anchorExpr', '', ...
                    'anchor', [], 'applyToAll', false, 'setDefault', false);
dfgetset('binDlgInfo', binDlgInfo);


