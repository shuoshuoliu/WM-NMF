function savelineproperties(fit)
%SAVELINEPROPERTIES Save current line properties for later recovery


% Copyright 2003-2008 The MathWorks, Inc.

% Get line properties previously saved
oldcml = fit.ColorMarkerLine;

cml = cell(4,1);
lineproperties = {'Color' 'Marker' 'LineStyle' 'LineWidth'};

if ~isempty(fit.linehandle) && ishghandle(fit.linehandle)
   cml = get(fit.linehandle, lineproperties);
end

% Only the discrete pdf uses markers, so don't overwrite the marker
% unless that's what we're plotting.
if length(oldcml) >= 4
   if ~isequal(fit.ftype,'pdf') || iscontinuous(fit)
      cml{2} = oldcml{2};
   end
end

% Save current properties
fit.ColorMarkerLine = cml;
