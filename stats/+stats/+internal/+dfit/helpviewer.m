function helpviewer(topic, errorname)
% HELPVIEWER  is a helper file for Distribution Fitter 
% HELPVIEWER Displays help for Distribution Fitter TOPIC. If the map file 
% cannot be found, an error is displayed using ERRORNAME

%   Copyright 2003-2020 The MathWorks, Inc.

mapfilename = [docroot '/toolbox/stats/stats.map'];
try
    helpview(mapfilename, topic);
catch
    errordlg(getString(message('stats:dfstrings:sprintf_UnableToDisplayHelp', errorname)));
end
