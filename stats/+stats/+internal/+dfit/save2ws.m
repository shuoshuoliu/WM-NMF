function save2ws(fitobj)
%SAVE2WS Utility to save probability distribution object from dfittool


%   Copyright 2003-2011 The MathWorks, Inc.

if ischar(fitobj)
    % Fit name passed in, so get fit object
    fitdb = stats.internal.dfit.getfitdb;
    fitobj = find(fitdb, 'name', fitobj);
end

if ~isa(fitobj,'stats.dffit')
    error(message('stats:dfsave2w:BadFit'));
end

% Bring up dialog to get variable name for this fit
export2wsdlg({getString(message('stats:dfstrings:dlg_SaveFittedDistributionAs'))},...
             {'pd'},{fitobj.probdist},...
             getString(message('stats:dfstrings:dlg_SaveFit2WS')));

end
