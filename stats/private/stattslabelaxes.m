function stattslabelaxes(ax,ts,yname,charttype)
%STATTSLABELAXES Label axes with time information.

% Copyright 2011 The MathWorks, Inc. 


% Get variable name and units, if any
if isempty(ts) || isequal(ts.Name,'unnamed') || isempty(ts.Name)
    varname = '';
else
    varname = ts.Name;
end
if isempty(ts) || isempty(ts.DataInfo.units)
    units = '';
else
    units = ts.DataInfo.units;
end

% Derive y axis label from varname and units
if isempty(yname)
    if isempty(varname)
        yname = getString(message('stats:stattslabelaxes:yaxisLabelMeasurement'));
    else
        yname = varname;
    end
end
if isempty(units)
    ylabel(ax,yname);
else
    ylabel(ax,sprintf('%s (%s)',yname,units));
end

% Create title with chart type
if isempty(varname)
    title(ax,charttype);
else
    title(ax,sprintf('%s for %s',charttype,varname));
end

% Label x axes with time or sample number
if isempty(ts)
    xlabel(getString(message('stats:stattslabelaxes:xaxisLabelSamples')));
elseif isempty(ts.TimeInfo.StartDate)
    if isempty(ts.TimeInfo.units)
        xlabel(ax,getString(message('stats:stattslabelaxes:xaxisLabelTime')));
    else
        xlabel(ax,getString(message('stats:stattslabelaxes:xaxisLabelTimeUnits',ts.TimeInfo.units)));
    end
else
    xlabel('');
    if isempty(ts.TimeInfo.Format)
        datetick(ax,'x');
    else
        datetick(ax,'x',ts.TimeInfo.Format);
    end
end
