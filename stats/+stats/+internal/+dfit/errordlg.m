function errorDialogFigure = errordlg(messageStr, titleStr)
% errordlg   Creates modal uifigure-based error dialog box with specified
% error message and title

%   Copyright 2019-2020 The MathWorks, Inc.

errorDialogFigure = uifigure( ...
    'Name', titleStr, ...
    'Visible', 'off', ...
    'Resize', 'on', ... 
    'Position', iGetDefaultPanelPosition(), ...
    'Tag', 'dfErrorDialogUIFigure', ...
    'WindowStyle', 'modal');

mainGrid = uigridlayout(errorDialogFigure);
mainGrid.ColumnWidth = {'1x'};
mainGrid.RowHeight = {'1x', 'fit'};
mainGrid.Padding = [5 10 0 10];

topGrid = uigridlayout(mainGrid);
topGrid.ColumnWidth = {'fit', '1x'};
topGrid.RowHeight = {'1x'};
topGrid.Padding = [0 0 0 0];
topGrid.Layout.Row = 1;
topGrid.Layout.Column = 1;

errorIconGrid = uigridlayout(topGrid);
errorIconGrid.ColumnWidth = {42};
errorIconGrid.RowHeight = {'1x', 41, '1x'};
errorIconGrid.Padding = [0 0 0 0];
errorIconGrid.RowSpacing = 0;
errorIconGrid.Layout.Row = 1;
errorIconGrid.Layout.Column = 1;

p = uipanel(errorIconGrid, ...
    'BorderType', 'none');
p.Layout.Row = 2;
p.Layout.Column = 1;

createErrorIcon(p);

messageLabel = uilabel(topGrid, ...
    'Text', messageStr, ...
    'HorizontalAlignment', 'left', ...
    'FontColor', [0.149 0.149 0.149], ...
    'WordWrap', 'on', ...
    'Tag', 'dfErrorDialogMessageLabel');
messageLabel.Layout.Row = 1;
messageLabel.Layout.Column = 2;

bottomGrid = uigridlayout(mainGrid);
bottomGrid.ColumnWidth = {'1x', 'fit', '1x'};
bottomGrid.RowHeight = {'fit'};
bottomGrid.Padding = [0 0 0 0];
bottomGrid.Layout.Row = 2;
bottomGrid.Layout.Column = 1;

iCreateDummyMinSizeButton(bottomGrid, 1, 2);

OKButton = uibutton(bottomGrid, ...
    'Text', getString((message('stats:dfittool:button_OK'))), ...
    'ButtonPushedFcn', @(~,~) delete(errorDialogFigure), ...
    'Tag', 'dfErrorDialogOKButton');
OKButton.Layout.Row = 1;
OKButton.Layout.Column = 2;

errorDialogFigure.Visible = 'on';
end

% helpers
function createErrorIcon(panelObj)
ax = axes(...
    'Parent', panelObj, ...
    'Units', 'normalized', ...
    'XTick',[], ...
    'YTick',[], ...
    'Color', iGetDefaultUIControlBackgroundColor(), ...
    'Tag', 'dfErrorDialogErrorIconAxes');


[iconData, alphaData] = matlab.ui.internal.dialog.DialogUtils.imreadDefaultIcon('error');
Img = image( ...
    'Parent',ax, ...
    'CData',iconData, ...
    'AlphaData', alphaData);

set(ax, ...
    'XLim', get(Img,'XData')+[-0.5 0.5], ...
    'YLim', get(Img,'YData')+[-0.5 0.5],  ...
    'Visible', 'off', ...
    'YDir', 'reverse');

ax.Toolbar = [];
disableDefaultInteractivity(ax);
end

function position = iGetDefaultPanelPosition()
width = 400;
height = 150;

position = get(gcbf, 'Position');
if isempty(position)
    position = get(0,'defaultfigureposition');
end

x = position(1) + (position(3) - width)/2;
y = position(2) + (position(4) - height)/2 - 50;

position = [x, y, width, height];
end

function rgb = iGetDefaultUIControlBackgroundColor()
rgb = get(0, 'defaultuicontrolbackgroundcolor');
end

function iCreateDummyMinSizeButton(grid, row, column)
dummyText = repmat('X',1,8);
minButton = uibutton(grid, 'Text', dummyText);
minButton.Layout.Row = row;
minButton.Layout.Column = column;
minButton.Visible = 'off';
end
