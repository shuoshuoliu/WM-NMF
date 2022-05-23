function dfupdateallplots(dods,dofit,force)
%DFUPDATEALLPLOTS Call update methods for all data sets and fits


%   Copyright 2003-2011 The MathWorks, Inc.

msg = '';
if nargin<3
   force = false;
end

% Supply defaults if missing or if called as a listener (args not logical)
if nargin<1 || ~islogical(dods)
   dods = true;
end
if nargin<2 || ~islogical(dofit)
   dofit = true;
end

% Get the array of data sets and update each one
if dods
   dsdb = getdsdb;
   ds = down(dsdb);
   while(~isempty(ds))
      try
         if (force)
             clearplot(ds);
         end    
         updateplot(ds);
      catch ME
         msg = appendmsg(msg,ds.name,ME.message);
      end
      ds = right(ds);
   end
end

% Get the array of fits and update each one
if dofit
   fitdb = getfitdb;
   ft = down(fitdb);
   while(~isempty(ft))
      try
         updateplot(ft);
      catch ME
         msg = appendmsg(msg,ft.name,ME.message);
      end
      ft = right(ft);
   end
end

if ~isempty(msg)
   errordlg(msg,getString(message('stats:dfstrings:dlg_ErrorUpdatingPlot')),'modal');
end


%--------------------------------------------
function msg = appendmsg(msg,objname,newmsg)
%APPENDMSG Append a new section to an existing set of error messages

if isempty(msg)
   msg = getString(message('stats:dfstrings:sprintf_ErrorPlotting0',objname,newmsg));
elseif ~isempty(newmsg)
   msg = getString(message('stats:dfstrings:sprintf_ErrorPlotting',msg,objname,newmsg));
end
