classdef ViewData < handle
    % ViewData   Dialog to display data set in a table
    
    %   Copyright 2019 The MathWorks, Inc.
    
    properties(Access = private)
        % Panel   (uifigure) The main window of the dialog
        Panel
        TopGrid
        LeftTopGrid
        
        % UI Components
        DataLabel
        CensoringLabel
        FrequencyLabel
        PreviewPanel
        ExclusionRuleDropDown
        DataSetTable
        CloseButton
        
        % DataSetObj   (UDD) dataSet object associated with the ViewData
        % object
        DataSetObj
        
        % Listeners   (cell of listener) Listeners to various events
        Listeners = {};
        
        % SessionClearingListener (listener) Listener to
        % DistributionFitting object's SessionClearing event
        SessionClearingListener
    end
    
    methods
        function this = ViewData(dataSetObj)
            this.DataSetObj = dataSetObj;
            
            this.createAndInitializePanel();
            
            this.SessionClearingListener = event.listener(stats.internal.dfit.DistributionFitting.getInstance(), 'SessionClearing', @(~,~) this.delete);
        end
        
        function showPanel(this)
            if ~isvalid(this.Panel) || isempty(this.Panel)
                this.createAndInitializePanel();
            end
            figure(this.Panel);
        end
        
        function delete(this)
            delete(this.Panel);
        end
                
        function updateDistributionPreviewAxes(this, ~, ~)
            if ~isvalid(this.Panel)
                return;
            else
                selectedExclusionRuleObj = find(iOutlierDatabase(), 'name', this.ExclusionRuleDropDown.Value);
                if isempty(selectedExclusionRuleObj)
                    stats.internal.dfit.DFViewUtilities.createDataSetDistributionPreview(this.PreviewPanel, this.DataSetObj);
                end
            end
        end
    end
    
    methods(Access = private)
        function createAndInitializePanel(this)
            this.createGUIComponents();
            
            excludeObj = stats.internal.dfit.Exclude.getInstance;
            binRulesManagerObj = stats.internal.dfit.BinRulesManager.getInstance;
            this.addListener(excludeObj, 'ExclusionRulesDeleted', @this.exclusionRuleAddedOrDeletedCallback);
            this.addListener(excludeObj, 'ExclusionRuleAdded', @this.exclusionRuleAddedOrDeletedCallback);
            this.addListener(excludeObj, 'ExclusionRuleRenamed', @this.exclusionRuleRenamedCallback);
            this.addListener(binRulesManagerObj, 'ApplyBinRuleToAllDataSets', @this.updateDistributionPreviewAxes);
            
            this.Panel.Visible = 'on';
        end
        
        function createGUIComponents(this)
            this.Panel = uifigure( ...
                'Name', getString(message('stats:dfittool:title_viewDataSet')), ...
                'Position', iGetDefaultPanelPosition(), ...
                'Visible', 'off',...
                'CloseRequestFcn', @(~, ~)this.closeButtonClickedCallback(), ...
                'Tag', 'dfViewDataUIFigure');
            
            mainGrid = uigridlayout(this.Panel);
            mainGrid.ColumnWidth = {'1x'};
            mainGrid.RowHeight = {'1x', 'fit'};
            
            this.TopGrid = uigridlayout(mainGrid, [1, 2]);
            this.TopGrid.Padding = [0 0 0 0];
            this.TopGrid.ColumnWidth = {'fit', '1x'};
            this.TopGrid.RowHeight = {'1x'};
            
            this.LeftTopGrid = uigridlayout(this.TopGrid, [6, 1]);
            this.LeftTopGrid.Padding = [0 0 0 0];
            this.LeftTopGrid.RowSpacing = 7;
            this.LeftTopGrid.ColumnWidth = {180};
            this.LeftTopGrid.RowHeight = {'fit','fit','fit',180,'fit','fit','1x'};
            
            bottomGrid = uigridlayout(mainGrid, [1, 2]);
            bottomGrid.Padding = [0 0 0 0];
            bottomGrid.Layout.Row = 2;
            bottomGrid.Layout.Column = 1;
            bottomGrid.ColumnWidth = {'1x', 'fit'};
            bottomGrid.RowHeight = {'fit'};
            
            this.createDataLabel();
            this.createCensoringLabel();
            this.createFrequencyLabel();
            this.createPreviewPanel([]);
            this.createExclusionRuleDropDown();
            this.createDataSetTable([]);
            
            this.createCloseButton(bottomGrid);
            
            stats.internal.dfit.centergui(this.Panel);
        end
        
        function updateViewData(this, exclusionRule)
            this.createPreviewPanel(exclusionRule);
            this.createDataSetTable(exclusionRule);
        end
        
        % UI Components
        function createDataLabel(this)
            this.DataLabel = uilabel(this.LeftTopGrid,...
                'Text', [getString(message('stats:dfittool:label_dataset')), ' ', this.DataSetObj.yexp], ...
                'Tag', 'dfViewDataDataLabel');
            this.DataLabel.Layout.Row = 1;
            this.DataLabel.Layout.Column = 1;
        end
        
        function createCensoringLabel(this)
            this.CensoringLabel = uilabel(this.LeftTopGrid,...
                'Text', [getString(message('stats:dfittool:label_censoring')), ' ', this.DataSetObj.censexp], ...
                'Tag', 'dfViewDataCensoringLabel');
            this.CensoringLabel.Layout.Row = 2;
            this.CensoringLabel.Layout.Column = 1;
        end
        
        function createFrequencyLabel(this)
            this.FrequencyLabel = uilabel(this.LeftTopGrid,...
                'Text', [getString(message('stats:dfittool:label_frequency')), ' ', this.DataSetObj.freqexp], ...
                'Tag', 'dfViewDataFrequencyLabel');
            this.FrequencyLabel.Layout.Row = 3;
            this.FrequencyLabel.Layout.Column = 1;
        end
       
        function createPreviewPanel(this, exclusionRule)
            this.PreviewPanel = uipanel(this.LeftTopGrid, ...
                'BackgroundColor', [1,1,1], ...
                'Tag', 'dfViewDataPreviewPanel');
            this.PreviewPanel.Layout.Row = 4;
            this.PreviewPanel.Layout.Column = 1;
            
            if isempty(exclusionRule)
               stats.internal.dfit.DFViewUtilities.createDataSetDistributionPreview(this.PreviewPanel, this.DataSetObj);
            else
               stats.internal.dfit.DFViewUtilities.createExclusionRulePreview(this.PreviewPanel, this.DataSetObj, exclusionRule);
            end
        end
        
        function createExclusionRuleDropDown(this)
            exclusionRuleLabel = uilabel(this.LeftTopGrid,...
                'Text', getString(message('stats:dfittool:label_exclusionRule')), ...
                'Tag', 'dfViewDataExclusionRuleLabel');
            exclusionRuleLabel.Layout.Row = 5;
            exclusionRuleLabel.Layout.Column = 1;
            
            this.ExclusionRuleDropDown = uidropdown(this.LeftTopGrid, ...
                'Items', iExclusionRuleList(), ...
                'ValueChangedFcn', @(~, ~) this.exclusionRuleDropDownValueChanged, ...
                'Tag', 'dfViewDataExclusionRuleDropDown');
            this.ExclusionRuleDropDown.Layout.Row = 6;
            this.ExclusionRuleDropDown.Layout.Column = 1;
        end
       
        function createDataSetTable(this, exclusionRule)
            if ~isempty(this.DataSetTable)
                delete(this.DataSetTable);
            end
            this.DataSetTable = stats.internal.dfit.DFViewUtilities.createDataSetTable(this.TopGrid, 1, 2, this.DataSetObj, exclusionRule);
            this.DataSetTable.Tag = 'dfViewDataDataSetTable';
        end
       
        function createCloseButton(this, grid)
            iCreateDummyMinSizeButton(grid, 1, 2)
            
            this.CloseButton = uibutton(grid,...
                'Text', getString(message('stats:dfittool:button_close')),...
                'ButtonPushedFcn', @(~, ~) this.closeButtonClickedCallback(), ...
                'Tag', 'dfViewDataCloseButton');
            this.CloseButton.Layout.Row = 1;
            this.CloseButton.Layout.Column = 2;
        end
        
        function addListener(this, obj, eventName, callback)
            this.Listeners{end+1} = event.listener(obj, eventName, callback);
        end
        
        function deleteListeners(this)
            cellfun(@(x) delete(x), this.Listeners);
            this.Listeners = {};
        end

        % Callbacks        
        function exclusionRuleDropDownValueChanged(this, ~, ~)
            selectedExclusionRuleObj = find(iOutlierDatabase(), 'name', this.ExclusionRuleDropDown.Value);
            this.updateViewData(selectedExclusionRuleObj);
        end
        
        function exclusionRuleAddedOrDeletedCallback(this, ~, ~)
            selectedExclusionRule = this.ExclusionRuleDropDown.Value;
            newExclusionRuleList = iExclusionRuleList();
            if ~ismember(selectedExclusionRule, newExclusionRuleList)
                this.ExclusionRuleDropDown.Value = iNone();
                this.exclusionRuleDropDownValueChanged();
            end
            this.ExclusionRuleDropDown.Items = newExclusionRuleList;
        end
        
        function exclusionRuleRenamedCallback(this, ~, renameEventData)
            selectedExclusionRule = this.ExclusionRuleDropDown.Value;
            
            oldExclusionRuleName = renameEventData.OldName;
            newExclusionRuleName = renameEventData.NewName;
            
            this.ExclusionRuleDropDown.Items = iExclusionRuleList();
            if strcmp(selectedExclusionRule, oldExclusionRuleName)
                this.ExclusionRuleDropDown.Value = newExclusionRuleName;
            end
        end
        
        function closeButtonClickedCallback(this, ~, ~)
            delete(this.Panel);
            this.deleteListeners();
        end
    end
end

% helpers
function str = iNone()
str = strcat('(',getString(message('stats:dfittool:label_none')),')');
end

function databaseObj = iOutlierDatabase()
databaseObj = stats.internal.dfit.getoutlierdb;
end

function names = iExclusionRuleList()
names = [iNone(),...
    stats.internal.dfit.getDataBaseNames(iOutlierDatabase())];
end

function position = iGetDefaultPanelPosition()
position = get(0,'defaultfigureposition');
position([3,4]) = [530 365];
end

function iCreateDummyMinSizeButton(grid, row, column)
dummyText = repmat('X',1,8);
minButton = uibutton(grid, 'Text', dummyText);
minButton.Layout.Row = row;
minButton.Layout.Column = column;
minButton.Visible = 'off';
end