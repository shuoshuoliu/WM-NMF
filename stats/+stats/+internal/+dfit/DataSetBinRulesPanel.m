classdef DataSetBinRulesPanel < handle
    % DataSetBinRulesPanel    Dialog to set data set bin rules
    
    %   Copyright 2019 The MathWorks, Inc.
    
    properties(Access = private)
        % Panel   (uifigure) The main window of the dialog
        Panel
        TopGrid
        TopLeftGrid
        TopRightGrid
        
        % UI components
        DataLabel
        BinRulesButtonGroup
        FreedmanDiaconisRuleRadioButton
        ScottRuleRadioButton
        NumberOfBinsRadioButton
        BinsCenteredOnIntegersRadioButton
        BinWidthRadioButton
        NumberOfBinsEditField
        BinWidthEditField
        BinWidthButtonGroup
        AutoBinPlacementRadioButton
        BinBoundaryAtRadioButton
        AutoBinPlacementLabel
        BinBoundaryAtLabel
        BinBoundaryAtEditField
        ApplyToAllCheckBox
        SaveAsDefaultCheckBox
        PreviewPanel
        
        % DataSetObj   (stats.dfdata) associated data set object in the
        % dialog
        DataSetObj
        
        % DataSetName   (string) Name of data set shown in the dialog
        DataSetName
        
        % BinRulesManager   (stats.internal.dfit.BinRulesManager)
        BinRulesManager
        
        % Listeners   (cell of listener) Listeners to various events
        Listeners = {};
        
        % SessionClearingListener (listener) Listener to
        % DistributionFitting's SessionClearing event
        SessionClearingListener
    end
    
    methods
        function this = DataSetBinRulesPanel(dataSet)
            this.DataSetObj = dataSet;
            this.DataSetName = dataSet.Name;
            this.BinRulesManager = stats.internal.dfit.BinRulesManager.getInstance();
            
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
        
        function radioButtonObj = getBinRuleRadioButton(this, binRuleIdx)
            switch binRuleIdx
                case 1
                    radioButtonObj = this.FreedmanDiaconisRuleRadioButton;
                case 2
                    radioButtonObj = this.ScottRuleRadioButton;
                case 3
                    radioButtonObj = this.NumberOfBinsRadioButton;
                case 4
                    radioButtonObj = this.BinsCenteredOnIntegersRadioButton;
                case 5
                    radioButtonObj = this.BinWidthRadioButton;
                otherwise
                    radioButtonObj = this.FreedmanDiaconisRuleRadioButton;
            end
        end
        
        function radioButtonObj = getBinWidthRadioButton(this, binWidthIdx)
            switch binWidthIdx
                case 1
                    radioButtonObj = this.AutoBinPlacementRadioButton;
                case 2
                    radioButtonObj = this.BinBoundaryAtRadioButton;
                otherwise
                    radioButtonObj = this.AutoBinPlacementRadioButton;
            end
        end
    end
    
    methods(Access = private)
        function createAndInitializePanel(this)
            this.createGUIComponents();
            this.initializeGUIStates();
            
            this.addListener(this.BinRulesManager, 'ApplyBinRuleToAllDataSets', @this.initializeGUIStates);
            this.addListener(stats.internal.dfit.Data.getInstance, 'DataSetRenamed', @this.dataSetRenamedCallback);
        end
        
        function createGUIComponents(this)
            this.Panel = uifigure( ...
                'Position', iGetDefaultPanelPosition(),...
                'Name', getString(message('stats:dfittool:title_binWidth')), ...
                'CloseRequestFcn', @(~,~) this.cancelButtonClickedCallback, ...
                'Tag', 'dfDataSetBinRulesPanelUIFigure', ...
                'Visible', 'off');
            
            mainGrid = uigridlayout(this.Panel, [5, 1]);
            mainGrid.RowHeight = {'1x', 1, 'fit', 1, 'fit'};
            
            this.TopGrid = uigridlayout(mainGrid, [1, 3]);
            this.TopGrid.Padding = [0 0 0 0];
            this.TopGrid.RowHeight = {'fit'};
            this.TopGrid.ColumnWidth = {'fit','1x',180};
            this.TopGrid.Layout.Row = 1;
            this.TopGrid.Layout.Column = 1;
            
            this.TopLeftGrid = uigridlayout(this.TopGrid, [8, 3]);
            this.TopLeftGrid.Padding = [0 0 0 0];
            this.TopLeftGrid.RowHeight = {'fit', 18, 18, 18, 18, 18, 18, 18};
            this.TopLeftGrid.ColumnWidth = {24, 24, 'fit'};
            this.TopLeftGrid.Layout.Row = 1;
            this.TopLeftGrid.Layout.Column = 1;
            
            this.TopRightGrid = uigridlayout(this.TopGrid, [3,1]);
            this.TopRightGrid.Padding = [0 0 0 0];
            this.TopRightGrid.RowHeight = {'fit',180,'1x'};
            this.TopRightGrid.ColumnWidth = {180};
            this.TopRightGrid.Layout.Row = 1;
            this.TopRightGrid.Layout.Column = 3;
            
            iCreateHorizontalLine(mainGrid, 2, 1);

            middleGrid = uigridlayout(mainGrid, [2, 1]);
            middleGrid.Padding = [0 0 0 0];
            middleGrid.RowHeight = {'fit', 'fit'};
            middleGrid.ColumnWidth = {'fit'};
            middleGrid.Layout.Row = 3;
            middleGrid.Layout.Column = 1;
            
            iCreateHorizontalLine(mainGrid, 4, 1);

            bottomGrid = uigridlayout(mainGrid, [1, 5]);
            bottomGrid.Padding = [0 0 0 0];
            bottomGrid.RowHeight = {'fit'};
            bottomGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit', 'fit'};
            bottomGrid.Layout.Row = 5;
            bottomGrid.Layout.Column = 1;
                        
            this.DataLabel = uilabel(this.TopLeftGrid, ...
                'Text', iCreateDataLabel(this.DataSetObj.name), ...
                'Tag', 'dfDataSetBinRulesPanelDataLabel');
            this.DataLabel.Layout.Row = 1;
            this.DataLabel.Layout.Column = [1 3];
            
            this.BinRulesButtonGroup = uibuttongroup(this.TopLeftGrid, ...
                'BorderType', 'none', ...
                'SelectionChangedFcn', @(~,~) this.radioButtonSelectionChangedCallback, ...
                'Tag', 'dfDataSetBinRulesPanelBinRulesButtonGroup');
            this.BinRulesButtonGroup.Layout.Row = [2 6];
            this.BinRulesButtonGroup.Layout.Column = 1;
            
            this.BinWidthButtonGroup = uibuttongroup(this.TopLeftGrid, ...
                'BorderType', 'none', ...
                'Tag', 'dfDataSetBinRulesPanelBinWidthButtonGroup');
            this.BinWidthButtonGroup.Layout.Row = [7 8];
            this.BinWidthButtonGroup.Layout.Column = 2;
            
            freedManDiaconisRuleLabel = uilabel(this.TopLeftGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_Freedman')));
            freedManDiaconisRuleLabel.Layout.Row = 2;
            freedManDiaconisRuleLabel.Layout.Column = [2 3];
            
            scottRuleLabel = uilabel(this.TopLeftGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_Scott')));
            scottRuleLabel.Layout.Row = 3;
            scottRuleLabel.Layout.Column = [2 3];
            
            numberOfBinsGrid = uigridlayout(this.TopLeftGrid, [1, 2]);
            numberOfBinsGrid.Padding = [0 0 0 0];
            numberOfBinsGrid.ColumnWidth = {'fit', '1x'};
            numberOfBinsGrid.RowHeight = {18};
            numberOfBinsGrid.Layout.Row = 4;
            numberOfBinsGrid.Layout.Column = [2 3];
            
            numberOfBinsLabel = uilabel(numberOfBinsGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_NumBins')));
            numberOfBinsLabel.Layout.Row = 1;
            numberOfBinsLabel.Layout.Column = 1;
            
            this.NumberOfBinsEditField = uieditfield(numberOfBinsGrid, ...
                'ValueChangingFcn', @(~,~) this.numberOfBinsEditFieldValueChangingCallback(), ...
                'Tag', 'dfDataSetBinRulesPanelNumberOfBinsEditField');
            this.NumberOfBinsEditField.Layout.Row = 1;
            this.NumberOfBinsEditField.Layout.Column = 2;
            
            binsCenteredOnIntegersLabel = uilabel(this.TopLeftGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_Integers')));
            binsCenteredOnIntegersLabel.Layout.Row = 5;
            binsCenteredOnIntegersLabel.Layout.Column = [2 3];
            
            binWidthGrid = uigridlayout(this.TopLeftGrid, [1, 2]);
            binWidthGrid.Padding = [0 0 0 0];
            binWidthGrid.ColumnWidth = {'fit', '1x'};
            binWidthGrid.RowHeight = {18};
            binWidthGrid.Layout.Row = 6;
            binWidthGrid.Layout.Column = [2 3];
            
            binWidthLabel = uilabel(binWidthGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_BinWidth')));
            binWidthLabel.Layout.Row = 1;
            binWidthLabel.Layout.Column = 1;
            
            this.BinWidthEditField = uieditfield(binWidthGrid, ...
                'ValueChangingFcn', @(~,~) this.binWidthEditFieldValueChangingCallback(), ...
                'Tag', 'dfDataSetBinRulesPanelBinWidthEditField');
            this.BinWidthEditField.Layout.Row = 1;
            this.BinWidthEditField.Layout.Column = 2;
            
            this.AutoBinPlacementLabel = uilabel(this.TopLeftGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_BinWidthAny')), ...
                'Tag', 'dfDataSetBinRulesPanelAutoBinPlacementlabel');
            this.AutoBinPlacementLabel.Layout.Row = 7;
            this.AutoBinPlacementLabel.Layout.Column = 3;
            
            binBoundaryAtGrid = uigridlayout(this.TopLeftGrid, [1, 2]);
            binBoundaryAtGrid.Padding = [0 0 0 0];
            binBoundaryAtGrid.ColumnWidth = {'fit', '1x'};
            binBoundaryAtGrid.RowHeight = {18};
            binBoundaryAtGrid.Layout.Row = 8;
            binBoundaryAtGrid.Layout.Column = 3;
            
            this.BinBoundaryAtLabel = uilabel(binBoundaryAtGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_BinWidthSpecify')), ...
                'Tag', 'dfDataSetBinRulesPanelBinBoundaryAtLabel');
            this.BinBoundaryAtLabel.Layout.Row = 1;
            this.BinBoundaryAtLabel.Layout.Column = 1;
            
            this.BinBoundaryAtEditField = uieditfield(binBoundaryAtGrid, ...
                'ValueChangingFcn', @(~,~) this.binBoundaryAtEditFieldValueChangingCallback(), ...
                'Tag', 'dfDataSetBinRulesPanelBinBoundaryAtEditField');
            this.BinBoundaryAtEditField.Layout.Row = 1;
            this.BinBoundaryAtEditField.Layout.Column = 2;
            
            this.FreedmanDiaconisRuleRadioButton = uiradiobutton(this.BinRulesButtonGroup, 'Text', '', 'Position', [10 112 150 20], 'Tag', 'FreedmanDiaconis');
            this.ScottRuleRadioButton = uiradiobutton(this.BinRulesButtonGroup, 'Text', '', 'Position', [10 84 150 20], 'Tag', 'ScottRule');
            this.NumberOfBinsRadioButton = uiradiobutton(this.BinRulesButtonGroup, 'Text', '', 'Position', [10 56 150 20], 'Tag', 'NumberOfBins');
            this.BinsCenteredOnIntegersRadioButton = uiradiobutton(this.BinRulesButtonGroup, 'Text', '', 'Position', [10 28 15 20], 'Tag', 'BinsCenteredOnIntegers');
            this.BinWidthRadioButton = uiradiobutton(this.BinRulesButtonGroup, 'Text', '', 'Position', [10 0 15 20], 'Tag', 'BinWidth');
            
            this.AutoBinPlacementRadioButton = uiradiobutton(this.BinWidthButtonGroup, 'Text', '', 'Position', [10 28 15 20], 'Tag', 'BinWidthAutoBinPlacement');
            this.BinBoundaryAtRadioButton = uiradiobutton(this.BinWidthButtonGroup, 'Text', '', 'Position', [10 0 15 20], 'Tag', 'BinWidthBinBoundaryAt');
            
            dataPreviewLabel = uilabel(this.TopRightGrid, ...
                'Text', getString(message('stats:dfittool:label_preview')));
            dataPreviewLabel.Layout.Row = 1;
            dataPreviewLabel.Layout.Column = 1;
            
            this.PreviewPanel = uipanel(this.TopRightGrid, ...
                'BackgroundColor', [1,1,1]);
            this.PreviewPanel.Layout.Row = 2;
            this.PreviewPanel.Layout.Column = 1;
            stats.internal.dfit.DFViewUtilities.createDataSetDistributionPreview(this.PreviewPanel, this.DataSetObj);
            
            this.ApplyToAllCheckBox = uicheckbox(middleGrid, ...
                'Text', getString(message('stats:dfittool:checkbox_applyToAll')), ...
                'Tag', 'dfDataSetBinRulesPanelApplyToAllCheckBox');
            this.ApplyToAllCheckBox.Layout.Row = 1;
            this.ApplyToAllCheckBox.Layout.Column = 1;
            
            this.SaveAsDefaultCheckBox = uicheckbox(middleGrid, ...
                'Text', getString(message('stats:dfittool:checkbox_saveAsDefault')), ...
                'Tag', 'dfDataSetBinRulesPanelSaveAsDefaultCheckBox');
            this.SaveAsDefaultCheckBox.Layout.Row = 2;
            this.SaveAsDefaultCheckBox.Layout.Column = 1;
            
            iCreateDummyMinSizeButton(bottomGrid, 1, 1);
            iCreateDummyMinSizeButton(bottomGrid, 1, 3);
            iCreateDummyMinSizeButton(bottomGrid, 1, 4);
            iCreateDummyMinSizeButton(bottomGrid, 1, 5);
            
            helpButton = uibutton(bottomGrid, ...
                'Text', getString(message('stats:dfittool:button_help')), ...
                'ButtonPushedFcn', @(~,~) this.helpButtonClickedCallback, ...
                'Tag', 'dfDataSetBinRulesPanelHelpButton');
            helpButton.Layout.Row = 1;
            helpButton.Layout.Column = 1;
            
            updatePreviewButton = uibutton(bottomGrid, ...
                'Text', getString(message('stats:dfittool:button_updatePreview')), ...
                'ButtonPushedFcn', @(~,~) this.updatePreviewButtonClickedCallback, ...
                'Tag', 'dfDataSetBinRulesPanelUpdatePreviewButton');
            updatePreviewButton.Layout.Row = 1;
            updatePreviewButton.Layout.Column = 3;
            
            oKButton = uibutton(bottomGrid, ...
                'Text', getString(message('stats:dfittool:button_OK')), ...
                'ButtonPushedFcn', @(~,~) this.oKButtonClickedCallback, ...
                'Tag', 'dfDataSetBinRulesPanelOKButton');
            oKButton.Layout.Row = 1;
            oKButton.Layout.Column = 4;
            
            cancelButton = uibutton(bottomGrid, ...
                'Text', getString(message('stats:dfittool:button_cancel')), ...
                'ButtonPushedFcn', @(~,~) this.cancelButtonClickedCallback, ...
                'Tag', 'dfDataSetBinRulesPanelCancelButton');
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 5;
            
            stats.internal.dfit.centergui(this.Panel);
            this.Panel.Visible = 'on';
        end
        
        function initializeGUIStates(this, ~, ~)
            [binRuleIdx, numBinValueExpr, binWidthValueExpr, binWidthIdx, ...
                binBoundaryAtValueExpr, isApplyBinRuleToAll, isSaveAsDefault] = ...
                iGetBinDialogState(this.DataSetObj);

            selectedBinRuleRadioButton = this.getBinRuleRadioButton(binRuleIdx);
            selectedBinWidthRadioButton = this.getBinWidthRadioButton(binWidthIdx);
            
            this.BinRulesButtonGroup.SelectedObject = selectedBinRuleRadioButton;
            this.BinWidthButtonGroup.SelectedObject = selectedBinWidthRadioButton;
            this.enableBinWidthRadioButton(isBinWidthRadioButton(selectedBinRuleRadioButton));
            this.NumberOfBinsEditField.Value = numBinValueExpr;
            this.BinWidthEditField.Value = binWidthValueExpr;
            this.BinBoundaryAtEditField.Value = binBoundaryAtValueExpr;
            this.ApplyToAllCheckBox.Value = isApplyBinRuleToAll;
            this.SaveAsDefaultCheckBox.Value = isSaveAsDefault;
            
            this.updatePreviewButtonClickedCallback();
        end
        
        function enableBinWidthRadioButton(this, tf)
            this.AutoBinPlacementRadioButton.Enable = tf;
            this.BinBoundaryAtRadioButton.Enable = tf;
            this.AutoBinPlacementLabel.Enable = tf;
            this.BinBoundaryAtLabel.Enable = tf;
            this.BinBoundaryAtEditField.Enable = tf;
        end
        
        function isValid = checkValidityOfBinRuleValues(this)            
            selectedBinRuleRadioButton = this.BinRulesButtonGroup.SelectedObject;
            selectedBinWidthRadioButton = this.BinWidthButtonGroup.SelectedObject;
            numBinValue = char(this.NumberOfBinsEditField.Value);
            binWidthValue = char(this.BinWidthEditField.Value);
            binBoundaryAtValue = char(this.BinBoundaryAtEditField.Value);
            
            [isValid, errMessage] = this.BinRulesManager.checkBinRulesValidity(...
                selectedBinRuleRadioButton.Tag, selectedBinWidthRadioButton.Tag,...
                numBinValue, binWidthValue, binBoundaryAtValue);
            
            if ~isValid
                stats.internal.dfit.errordlg(errMessage, getString(message('stats:dfittool:title_invalidValues')));
            end
        end
        
        % Callbacks
        function dataSetRenamedCallback(this, ~, renameEventData)
            oldName = renameEventData.OldName;
            newName = renameEventData.NewName;
            if strcmp(this.DataSetName, oldName)
                this.DataSetName = newName;
                this.DataLabel.Text = iCreateDataLabel(newName);
            end
        end
        
        function numberOfBinsEditFieldValueChangingCallback(this, ~, ~)
            this.BinRulesButtonGroup.SelectedObject = this.NumberOfBinsRadioButton;
            this.radioButtonSelectionChangedCallback();
        end
        
        function binWidthEditFieldValueChangingCallback(this, ~, ~)
            this.BinRulesButtonGroup.SelectedObject = this.BinWidthRadioButton;
            this.radioButtonSelectionChangedCallback();
        end
        
        function binBoundaryAtEditFieldValueChangingCallback(this, ~, ~)
            this.BinWidthButtonGroup.SelectedObject = this.BinBoundaryAtRadioButton;
        end
        
        function radioButtonSelectionChangedCallback(this, ~, ~)
            tf = isBinWidthRadioButton(this.BinRulesButtonGroup.SelectedObject);
            this.enableBinWidthRadioButton(tf);
        end
        
        function updatePreviewButtonClickedCallback(this, ~, ~)
            isValid = this.checkValidityOfBinRuleValues();
            if ~isValid
                return;
            end
            
            selectedBinRuleRadioButton = this.BinRulesButtonGroup.SelectedObject;
            selectedBinWidthRadioButton = this.BinWidthButtonGroup.SelectedObject;
            
            binInfo.rule = iGetBinRuleRadioButtonIndex(selectedBinRuleRadioButton);
            binInfo.nbins = str2double(this.NumberOfBinsEditField.Value);
            binInfo.width = str2double(this.BinWidthEditField.Value);
            binInfo.placementRule = iGetBinWidthRadioButtonIndex(selectedBinWidthRadioButton);
            binInfo.anchor = str2double(this.BinBoundaryAtEditField.Value);
            
            stats.internal.dfit.DFViewUtilities.createDataSetDistributionPreview(this.PreviewPanel, this.DataSetObj, binInfo);
        end
        
        function oKButtonClickedCallback(this, ~, ~)
            isValid = this.checkValidityOfBinRuleValues();
            if ~isValid
                return;
            end
            
            selectedBinRuleRadioButton = this.BinRulesButtonGroup.SelectedObject;
            selectedBinWidthRadioButton = this.BinWidthButtonGroup.SelectedObject;
            
            binRuleIdx = iGetBinRuleRadioButtonIndex(selectedBinRuleRadioButton);
            numBinValue = char(this.NumberOfBinsEditField.Value);
            binWidthValue = char(this.BinWidthEditField.Value);
            binWidthIdx = iGetBinWidthRadioButtonIndex(selectedBinWidthRadioButton);
            binBoundaryAtValue = char(this.BinBoundaryAtEditField.Value);
            applyToAll = this.ApplyToAllCheckBox.Value;
            saveAsDefault = this.SaveAsDefaultCheckBox.Value;
                        
            stats.internal.dfit.setbinwidthrules(this.DataSetObj, ...
                binRuleIdx, numBinValue, ...
                binWidthValue, ...
                binWidthIdx, binBoundaryAtValue, ...
                applyToAll, saveAsDefault);
            
            if applyToAll
                this.BinRulesManager.applyBinRuleToAllDataSets();
            else
                this.BinRulesManager.applyBinRuleToDataSet(this.DataSetObj);
            end
            
            if saveAsDefault
                this.BinRulesManager.saveBinRuleAsDefault();
            end
            
            this.cancelButtonClickedCallback();
        end
        
        function cancelButtonClickedCallback(this, ~, ~)
            delete(this.Panel);
            this.deleteListeners();
        end
        
        function helpButtonClickedCallback(this, ~, ~) %#ok<INUSD>
            stats.internal.dfit.helpviewer('set_bin_rules','setting bin rules');
        end
        
        function addListener(this, obj, eventName, callback)
            this.Listeners{end+1} = event.listener(obj, eventName, callback);
        end
        
        function deleteListeners(this)
            cellfun(@(x) delete(x), this.Listeners);
            this.Listeners = {};
        end
    end
end

% helper functions
function position = iGetDefaultPanelPosition()
position = get(0,'defaultfigureposition');
position([3,4]) = [500 351];
end

function iCreateHorizontalLine(grid, row, column)
p = uipanel(grid);
p.BackgroundColor = [.75 .75 .75];
p.Layout.Row = row;
p.Layout.Column = column;
end

function iCreateDummyMinSizeButton(grid, row, column)
dummyText = repmat('X',1,8);
minButton = uibutton(grid, 'Text', dummyText);
minButton.Layout.Row = row;
minButton.Layout.Column = column;
minButton.Visible = 'off';
end

function str = iCreateDataLabel(dataSetName)
if isempty(dataSetName)
    str = "";
else
    str = strcat(getString(message('stats:dfittool:label_data')), " ", dataSetName);
end
end

function idx = iGetBinRuleRadioButtonIndex(radioButtonObj)
binRuleTag = radioButtonObj.Tag;
switch binRuleTag
    case 'FreedmanDiaconis'
        idx = 1;
    case 'ScottRule'
        idx = 2;
    case 'NumberOfBins'
        idx = 3;
    case 'BinsCenteredOnIntegers'
        idx = 4;
    case 'BinWidth'
        idx = 5;
    otherwise
        idx = 1;
end
end

function idx = iGetBinWidthRadioButtonIndex(radioButtonObj)
binWidthTag = radioButtonObj.Tag;
switch binWidthTag
    case 'BinWidthAutoBinPlacement'
        idx = 1;
    case 'BinWidthBinBoundaryAt'
        idx = 2;
    otherwise
        idx = 1;
end
end

function [ruleName, numBinsValueExpr, binWidthValueExpr, binWidthIdx, binBoundaryAtValueExpr, isApplyBinRuleToAll, isSaveAsDefault] = iGetBinDialogState(dataSetObj)
if isempty(dataSetObj)
    [ruleName, ...
        numBinsValueExpr, ...
        binWidthValueExpr, ...
        binWidthIdx, ...
        binBoundaryAtValueExpr, ...
        isApplyBinRuleToAll, ...
        isSaveAsDefault] = dfgetbinwidthdefaults();
else
    binDlgInfo = dataSetObj.binDlgInfo;
    ruleName = binDlgInfo.rule;
    numBinsValueExpr = binDlgInfo.nbinsExpr;
    binWidthValueExpr = binDlgInfo.widthExpr;
    binWidthIdx = binDlgInfo.placementRule;
    binBoundaryAtValueExpr = binDlgInfo.anchorExpr;
    isApplyBinRuleToAll = binDlgInfo.applyToAll;
    isSaveAsDefault = binDlgInfo.setDefault;
end
end

function tf = isBinWidthRadioButton(radioButtonObj)
tf = strcmp(radioButtonObj.Tag, 'BinWidth');
end

