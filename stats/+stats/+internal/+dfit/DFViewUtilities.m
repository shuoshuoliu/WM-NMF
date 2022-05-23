classdef DFViewUtilities < handle
    % DFViewUtilities   utility functions that help create UI components
    % for ViewData and ViewExclusionRule classes
    
    %   Copyright 2019-2020 The MathWorks, Inc.
    
    methods(Static)
        function dataSetTable = createDataSetTable(grid, layoutRow, layoutColumn, dataSet, exclusionRule)
            dataSetTable = uitable(grid, 'BackgroundColor', iWhite());
            dataSetTable.Data = iGetData(dataSet);
            dataSetTable.RowName = [];
            dataSetTable.ColumnFormat = {'numeric', 'numeric','numeric', 'numeric'};
            dataSetTable.ColumnName = {...
                iGetMessageString('stats:dfittool:columnLabel_index'), ...
                iGetMessageString('stats:dfittool:columnLabel_data'), ...
                iGetMessageString('stats:dfittool:columnLabel_censoring'), ...
                iGetMessageString('stats:dfittool:columnLabel_frequency')};
            if ~isempty(exclusionRule)
                excludedIdx = iIsExcluded(dataSet.y,exclusionRule.YLowLessEqual, ...
                    exclusionRule.getlowerbound, exclusionRule.YHighGreaterEqual, ...
                    exclusionRule.getupperbound);
                dataSetTable.BackgroundColor = iGetTableBackgroundColor(excludedIdx);
            end
            dataSetTable.Layout.Row = layoutRow;
            dataSetTable.Layout.Column = layoutColumn;
        end
        
        function previewAxes = createExclusionRulePreview(panelObj, dataSet, exclusionRule)
            if isempty(dataSet)
                previewAxes = iCreateExclusionPreviewAxesWithoutData(panelObj, exclusionRule);
            else
                previewAxes = iCreateExclusionPreviewAxesWithData(panelObj, dataSet, exclusionRule);
            end
            previewAxes.Toolbar = [];
            disableDefaultInteractivity(previewAxes);
        end
        
        function previewAxes = createDataSetDistributionPreview(panelObj, dataSet, binInfo)
            if nargin < 3
                binInfo = dataSet.binDlgInfo;
            end
            previewAxes = iCreateDataSetPreviewAxes(panelObj, dataSet, binInfo);
            previewAxes.Toolbar = [];
            disableDefaultInteractivity(previewAxes);
        end
        
        function isValidDataSet = createDataSetDistributionPreviewFromExpressions(panelObj, dataExpr, censoringExpr, frequencyExpr)            
            [errMessage, data, censoring, frequency] = stats.internal.dfit.checkselections(dataExpr, censoringExpr, frequencyExpr);
            
            if isempty(errMessage)
                dataSet.y = data;
                dataSet.censored = censoring;
                dataSet.frequency = frequency;
                
                % Use default bin rules
                defaultBinRule = dfgetset('binDlgInfo');
                
                stats.internal.dfit.DFViewUtilities.createDataSetDistributionPreview(panelObj, dataSet, defaultBinRule);
                isValidDataSet = true;
            else
                errMessage = [getString(message('stats:dfittool:title_invalidSelections')), ...
                    newline, ...
                    newline, ...
                    errMessage];

                stats.internal.dfit.DFViewUtilities.createMessagePreview(panelObj, errMessage);
                isValidDataSet = false;
            end
        end
        
        function createMessagePreview(panelObj, messageStr)
            delete(panelObj.Children);
            
            scrollableGrid = uigridlayout(panelObj, ...
                'ColumnWidth', {'fit'}, ...
                'RowHeight', {'1x'}, ...
                'Scrollable', 'on');
            
            uilabel(scrollableGrid, ...
                'Text', messageStr, ...
                'VerticalAlignment', 'Top');
            
            scroll(scrollableGrid,'left', 'top');
        end
    end
end

% helpers
function outputStr = iGetMessageString(inputStr)
outputStr = getString(message(inputStr));
end

function d = iGetData(dataSetObj)
d = cell(dataSetObj.ylength, 4);
for i = 1:dataSetObj.ylength
    d{i, 1} = i;
    d{i, 2} = dataSetObj.y(i);
    if isempty(dataSetObj.censored)
        d{i, 3} = [];
    else
        d{i, 3} = dataSetObj.censored(i);
    end
    if isempty(dataSetObj.frequency)
        d{i, 4} = [];
    else
        d{i, 4} = dataSetObj.frequency(i);
    end
end
end

function boolArray = iIsExcluded(dataVector, isLessEqual, lowerLimit, isGreaterEqual, upperLimit)
if isLessEqual
    excludedLowerBound = dataVector <= lowerLimit;
else
    excludedLowerBound = dataVector < lowerLimit;
end

if isGreaterEqual
    excludedUpperBound = dataVector >= upperLimit;
else
    excludedUpperBound = dataVector > upperLimit;
end
boolArray = excludedLowerBound | excludedUpperBound;
end

function colorArray = iGetTableBackgroundColor(boolArray)
colorArray = ones(numel(boolArray), 3);
colorArray(boolArray,:) = repmat(iDarkGrey(), sum(boolArray), 1);
end

function ax = iCreateExclusionPreviewAxesWithData(panelObj, dataSet, exclusionRule)
% Get data w/o NaNs but with no exclusion rule applied
[ydata, cens, freq] = getincludeddata(dataSet, []);
ydata = real(ydata);
if isempty(cens)
    cens = zeros(size(ydata));
end
if isempty(freq)
    freq = ones(size(ydata));
end

% Sort y and carry along the rest
[ydata,i] = sort(ydata);
cens = cens(i);
freq = freq(i);

% Create x and y vectors to plot
n = sum(freq);
x = zeros(n,1);
y = zeros(n,1);
g = zeros(n,1);
x(1:freq(1)) = ydata(1);
y(1:freq(1)) = (1:freq(1))';
g(1:freq(1)) = cens(1);
i = freq(1)+1;
for k=2:length(ydata)
    for j=1:freq(k)
        x(i) = ydata(k);
        g(i) = cens(k);
        if (i>1) && (x(i)==x(i-1))
            y(i) = y(i-1) + 1;
        else
            y(i) = 1;
        end
        i = i+1;
    end
end

ylo = exclusionRule.YLow;
if isempty(ylo)
    ylo = -Inf;
else
    ylo = str2double(ylo);
end

yhi = exclusionRule.YHigh;
if isempty(yhi)
    yhi = Inf;
else
    yhi = str2double(yhi);
end

ylotest = exclusionRule.YLowLessEqual;
yhitest = exclusionRule.YLowLessEqual;

xlim = [min(x) max(x)];
xlim = xlim + .05 * [-1 1] * diff(xlim);
ylim = [min(y) max(y)];
ylim = ylim + .05 * [-1 1] * diff(ylim);
if ylim(1) == ylim(2)
    ylim = [0 2];
end

ax = axes(...
    'Parent', panelObj, ...
    'Xtick',[],...
    'Ytick',[], ...
    'XLim', xlim, ...
    'YLim', ylim, ...
    'Box', 'on', ...
    'Color', iWhite(), ...
    'Units', 'normalized', ...
    'Position', iAxesPosition(), ...
    'Visible','off');

if ylotest==0
    inbounds = x>=ylo;
else
    inbounds = x>ylo;
end

if yhitest==0
    inbounds = inbounds & x<=yhi;
else
    inbounds = inbounds & x<yhi;
end

t = inbounds;
l1 = line('XData',x(t), 'YData',y(t),...
    'Color', 'b', 'Marker', '.','LineStyle', 'none', 'Parent', ax);
t = ~inbounds;
l2 = line('XData',x(t), 'YData',y(t),...
    'Color', iDefaultFigureColor()/2, 'Marker', '.', 'LineStyle', 'none', 'Parent', ax);
alllines = [l1 l2];

% Create patches to show the area outside the domain and range
allpatches = [];

if (ylo ~= -Inf)
    xlo = max(ylo, xlim(1));
    xp = [xlim(1) xlo xlo xlim(1)];
    yp = [ylim(1) ylim(1) ylim(2) ylim(2)];
    p1 = patch(xp, yp, iLightGrey(), 'LineStyle', 'none', 'Parent', ax);
    if isempty(allpatches)
        allpatches = p1;
    else
        allpatches(end+1)= p1;
    end
end

if (yhi ~= Inf)
    xhi = min(yhi,xlim(2));
    xp = [xlim(2) xhi xhi xlim(2)];
    yp = [ylim(1) ylim(1) ylim(2) ylim(2)];
    p2 = patch(xp, yp, iLightGrey(), 'LineStyle', 'none', 'Parent', ax);
    if isempty(allpatches)
        allpatches = p2;
    else
        allpatches(end+1)= p2;
    end
end

set(ax,'Children',[alllines allpatches]);
end

function ax = iCreateExclusionPreviewAxesWithoutData(panelObj, exclusionRule)
xlim = [0 4];
ylim = [0 4];

ax = axes(...
    'Parent', panelObj, ...
    'XTick',[], ...
    'YTick',[], ...
    'XLim',xlim, ...
    'YLim',ylim, ...
    'Box','on', ...
    'Color', iWhite(), ...
    'Units', 'normalized', ...
    'Position', [0.05, 0.05, 0.9, 0.9], ...
    'Visible','off');

xlo = exclusionRule.YLow;
xhi = exclusionRule.YHigh;

if ~isempty(xlo)
    patch([0 1 1 0], [0 0 4 4], iLightGrey(), 'LineStyle', 'none', 'Parent', ax);
end
if ~isempty(xhi)
    patch([3 4 4 3], [0 0 4 4], iLightGrey(), 'LineStyle', 'none','Parent', ax);
end
end

function ax = iCreateDataSetPreviewAxes(panelObj, dataSet, binInfo)

data = dataSet.y;
censoring = dataSet.censored;
frequency = dataSet.frequency;

delete(panelObj.Children);

ax = axes( ...
    'Parent', panelObj, ...
    'Box','on', ...
    'Color', iWhite(), ...
    'Units', 'normalized', ...
    'Position', iAxesPosition(), ...
    'Visible','off');

% If data has a complex part, it will spit a warning to the command line, so
% turn off warnings before plotting.
warnstate=warning('off', 'all');

% If we're working on expressions rather than data in an existing data set,
% we may need to remove NaNs
[~,~,data,censoring,frequency] = internal.stats.removenan(data,censoring,frequency);

% Compute the bin centers using the ecdf
% to allow a quartile computation even when there is censoring.
[fstep, xstep] = ecdf(data, 'censoring', censoring, 'frequency', frequency);
[~,binEdges] = internal.stats.histbins(data,censoring,frequency,binInfo,fstep,xstep);

% Plot a histogram from ecdf using the computed number of bins
ecdfhist(ax, fstep, xstep, 'edges', binEdges);
set(ax, 'xtick',[],'ytick',[]);
axis(ax,'tight');
allchildren = get(ax, 'children');
patchchildren = findobj(allchildren,'flat','Type','patch');
set(patchchildren, 'facecolor', iLightGrey());
warning(warnstate);

end

function rgb = iDefaultFigureColor()
rgb = get(0,'defaultuicontrolbackgroundcolor');
end

function rgb = iWhite()
rgb = [1, 1, 1];
end

function rgb = iLightGrey()
rgb = [0.9, 0.9, 0.9];
end

function rgb = iDarkGrey()
rgb = [0.65, 0.65, 0.65];
end

function position = iAxesPosition()
position = [0.05, 0.05, 0.9, 0.9];
end
