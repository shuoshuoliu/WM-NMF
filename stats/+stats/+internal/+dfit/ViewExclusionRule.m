classdef ViewExclusionRule < handle
    % ViewExclusionRule   Dialog to see the data points excluded by an
    % exclusion rule
    
    %   Copyright 2019 The MathWorks, Inc.
    
    properties(Access = private)
        % Panel   (uifigure) The main window of the dialog
        Panel
        
        % UI Components
        ExclusionRuleLabel
        DataLabel
        ExclusionRulePreviewPanel
        LowerLimitLabel
        UpperLimitLabel
        DataSetTable
        CloseButton
        
        % ExclusionRuleName   (string) Name of exclusion rule shown in the
        % dialog
        ExclusionRuleName
        
        % DataSetName   (string) Name of data set shown in the dialog
        DataSetName
        
        % Listeners   (cell of listener) Listeners to various events
        Listeners = {}
        
        % SessionClearingListener (listener) Listener to SessionClearing
        % event in DistributionFitting object
        SessionClearingListener
    end
    
    methods
        function this = ViewExclusionRule(exclusionRule)
            this.ExclusionRuleName = exclusionRule.name;
            this.DataSetName = exclusionRule.dataSet;
            dataSet = find(stats.internal.dfit.getdsdb, 'name', this.DataSetName);
            
            this.createGUIComponents(exclusionRule, dataSet);
            
            this.addListener(stats.internal.dfit.Exclude.getInstance, 'ExclusionRuleRenamed', @this.exclusionRuleRenamedCallback);
            this.addListener(stats.internal.dfit.Data.getInstance, 'DataSetRenamed', @this.dataSetRenamedCallback);
            
            this.SessionClearingListener = event.listener(stats.internal.dfit.DistributionFitting.getInstance(), 'SessionClearing', @(~,~) this.delete);
        end
        
        function showPanel(this)
            if strcmp(this.Panel.Visible, 'off')
                this.resetPanelPositionAndSize();
            end
            figure(this.Panel);
        end
        
        function delete(this)
            delete(this.Panel);
        end
    end
    
    methods(Access = private)
        function createGUIComponents(this, exclusionRule, dataSet)
            hasData = ~isempty(dataSet);
            
            this.Panel = uifigure(...
                'Name', iGetMessageString('stats:dfittool:title_viewExclusionRule'), ...
                'Position', iGetDefaultPanelPosition(hasData), ...
                'Visible', 'off', ...
                'CloseRequestFcn', @this.closeButtonClickedCallback, ...
                'Tag', 'dfViewExclusionRuleUIFigure');
            
            mainGrid = uigridlayout(this.Panel);
            mainGrid.ColumnWidth = {'1x'};
            mainGrid.RowHeight = {'1x', 'fit'};
            
            upperGrid = uigridlayout(mainGrid);
            upperGrid.Padding = [0 0 0 0];
            upperGrid.ColumnWidth = {'fit', '1x'};
            upperGrid.RowHeight = {'1x'};
            upperGrid.Layout.Row = 1;
            upperGrid.Layout.Column = 1;
            
            leftUpperGrid = uigridlayout(upperGrid);
            leftUpperGrid.Padding = [0 0 0 0];
            leftUpperGrid.RowSpacing = 7;
            leftUpperGrid.ColumnWidth = {180};
            leftUpperGrid.RowHeight = {'fit', 'fit', 180, 'fit', 'fit', 'fit', '1x'};
            leftUpperGrid.Layout.Row = 1;
            leftUpperGrid.Layout.Column = 1;
            
            lowerGrid = uigridlayout(mainGrid);
            lowerGrid.Padding = [1 1 1 1];
            lowerGrid.ColumnWidth = {'1x', 'fit'};
            lowerGrid.RowHeight = {'fit'};
            lowerGrid.Layout.Row = 2;
            lowerGrid.Layout.Column = 1;
            
            this.createExclusionRuleLabel(leftUpperGrid, this.ExclusionRuleName);
            this.createDataLabel(leftUpperGrid, this.DataSetName);
            this.createExclusionRulePreviewPanel(leftUpperGrid, exclusionRule, dataSet);
            this.createExcludedSectionsLabel(leftUpperGrid);
            this.createLowerLimitLabel(leftUpperGrid, exclusionRule.YLowLessEqual, exclusionRule.YLow);
            this.createUpperLimitLabel(leftUpperGrid, exclusionRule.YHighGreaterEqual, exclusionRule.YHigh);
            if hasData
                this.createDataSetTable(upperGrid, exclusionRule, dataSet);
            end
            this.createCloseButton(lowerGrid);
            
            stats.internal.dfit.centergui(this.Panel);
        end
        
        function resetPanelPositionAndSize(this)
            hasData = ~isempty(this.DataSetTable);
            this.Panel.Position = iGetDefaultPanelPosition(hasData);
            stats.internal.dfit.centergui(this.Panel);
        end
        
        % UI components
        function createExclusionRuleLabel(this, grid, exclusionRuleName)
            this.ExclusionRuleLabel = uilabel(grid, ...
                'Text', iCreateExclusionRuleLabelString(exclusionRuleName), ...
                'Tag', 'dfViewExclusionRuleExclusionRuleLabel');
            this.ExclusionRuleLabel.Layout.Row = 1;
            this.ExclusionRuleLabel.Layout.Column = 1;
        end
        
        function createDataLabel(this, grid, dataSetName)
            this.DataLabel = uilabel(grid, ...
                'Text', iCreateDataLabelString(dataSetName), ...
                'Tag', 'dfViewExclusionRuleDataLabel');
            this.DataLabel.Layout.Row = 2;
            this.DataLabel.Layout.Column = 1;
        end
        
        function createExclusionRulePreviewPanel(this, grid, exclusionRule, dataSet)
            this.ExclusionRulePreviewPanel = uipanel(grid, ...
                'BackgroundColor', [1,1,1], ...
                'Tag', 'dfViewExclusionRuleExclusionRulePreviewPanel');
            this.ExclusionRulePreviewPanel.Layout.Row = 3;
            this.ExclusionRulePreviewPanel.Layout.Column = 1;
            stats.internal.dfit.DFViewUtilities.createExclusionRulePreview(this.ExclusionRulePreviewPanel, dataSet, exclusionRule);
        end
        
        function createExcludedSectionsLabel(~, grid)
            excludedSectionsLabel = uilabel(grid, ...
                'Text', strcat(iGetMessageString('stats:dfittool:label_excludeSections'),":"));
            excludedSectionsLabel.Layout.Row = 4;
            excludedSectionsLabel.Layout.Column = 1;
        end
        
        function createLowerLimitLabel(this, grid, hasEqualSign, limitValue)
            if hasEqualSign
                limitSign = '<=';
            else
                limitSign = '<';
            end
            this.LowerLimitLabel = uilabel(grid, ...
                'Text', iCreateExcludeDataString(limitSign, limitValue), ...
                'Tag', 'dfViewExclusionRuleLowerLimitLabel');
            this.LowerLimitLabel.Layout.Row = 5;
            this.LowerLimitLabel.Layout.Column = 1;
        end
        
        function createUpperLimitLabel(this, grid, hasEqualSign, limitValue)
            if hasEqualSign
                limitSign = '>=';
            else
                limitSign = '>';
            end
            this.LowerLimitLabel = uilabel(grid, ...
                'Text', iCreateExcludeDataString(limitSign, limitValue), ...
                'Tag', 'dfViewExclusionRuleUpperLimitLabel');
            this.LowerLimitLabel.Layout.Row = 6;
            this.LowerLimitLabel.Layout.Column = 1;
        end
        
        function createDataSetTable(this, grid, exclusionRule, dataSet)
            this.DataSetTable = stats.internal.dfit.DFViewUtilities.createDataSetTable(grid, 1, 2, dataSet, exclusionRule);
            this.DataSetTable.Tag = 'dfViewExclusionRuleDataSetTable';
        end
        
        function createCloseButton(this, grid)
            iCreateDummyMinSizeButton(grid, 1, 2);
            
            this.CloseButton = uibutton(grid, ...
                'Text', iGetMessageString('stats:dfittool:button_close'), ...
                'ButtonPushedFcn', @(~, ~) this.closeButtonClickedCallback);
            this.CloseButton.Layout.Row = 1;
            this.CloseButton.Layout.Column = 2;
        end
        
        % callbacks
        function exclusionRuleRenamedCallback(this, ~, renameEventData)
            oldName = renameEventData.OldName;
            newName = renameEventData.NewName;
            if strcmp(this.ExclusionRuleName, oldName)
                this.ExclusionRuleName = newName;
                this.ExclusionRuleLabel.Text  = iCreateExclusionRuleLabelString(newName);
            end
        end
        
        function dataSetRenamedCallback(this, ~, renameEventData)
            oldName = renameEventData.OldName;
            newName = renameEventData.NewName;
            if strcmp(this.DataSetName, oldName)
                this.DataSetName = newName;
                this.DataLabel.Text = iCreateDataLabelString(newName);
            end
        end
        
        function closeButtonClickedCallback(this, ~, ~)
            this.Panel.Visible = 'off';
        end
        
        function addListener(this, obj, eventName, callback)
            this.Listeners{end+1} = event.listener(obj, eventName, callback);
        end
    end
end

% helper functions
function outputStr = iGetMessageString(inputStr)
outputStr = getString(message(inputStr));
end

function str = iCreateExclusionRuleLabelString(nameStr)
exclusionRuleLabel = iGetMessageString('stats:dfittool:label_exclusionRule');
str = strcat(exclusionRuleLabel, " ", nameStr);
end

function str = iCreateDataLabelString(nameStr)
dataLabel = iGetMessageString('stats:dfittool:label_data');
if isempty(nameStr)
    nameStr = iNone();
end
str = strcat(dataLabel, " ", nameStr);
end

function str = iCreateExcludeDataString(signStr, valueStr)
excludeDataLabel = iGetMessageString('stats:dfittool:label_restrictY');
if isempty(valueStr)
    str = "";
else
    str = strcat(iIndentationPadding(), excludeDataLabel, " ", signStr, " ", valueStr);
end
end

function str = iIndentationPadding()
str = "      ";
end

function str = iNone()
str = strcat('(',iGetMessageString('stats:dfittool:label_none'),')');
end

function position = iGetDefaultPanelPosition(hasData)
position = get(0,'defaultfigureposition');
if hasData
    position([3,4]) = [530 365];
else
    position([3,4]) = [230 365];
end
end

function iCreateDummyMinSizeButton(grid, row, column)
dummyText = repmat('X',1,8);
minButton = uibutton(grid, 'Text', dummyText);
minButton.Layout.Row = row;
minButton.Layout.Column = column;
minButton.Visible = 'off';
end




