function ok = dfasksavesession(~)
%DFASKSAVESESSION Ask whether current session should be saved


%   Copyright 2003-2011 The MathWorks, Inc.

dsdb = getdsdb;
fitdb = getfitdb;

yesText = getString(message('stats:dfstrings:button_Yes'));
noText = getString(message('stats:dfstrings:button_No'));
cancelText = getString(message('stats:dfstrings:button_Cancel'));

% Offer to save session unless there's nothing to save or the
% session has not changed
if ~dfgetset('dirty') || (isempty(down(dsdb)) && isempty(down(fitdb)))
   resp = noText;
else
   resp = questdlg(getString(message('stats:dfstrings:dlg_SaveThisDistributionFittingSession')), ...
                   getString(message('stats:dfstrings:dlg_DistributionFitting')), ...
                   yesText, noText, cancelText,    yesText);
end

if isempty(resp)
	resp = cancelText;
end

if isequal(resp,yesText)
   ok = dfsession('save');
   if ~ok
      resp = cancelText;
   end
end

ok = ~isequal(resp,cancelText);
