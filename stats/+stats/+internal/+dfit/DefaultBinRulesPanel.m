classdef DefaultBinRulesPanel < handle
    % DefaultBinRulesPanel    Dialog to set default bin rules for any new
    % data sets created
    
    %   Copyright 2019-2020 The MathWorks, Inc.
    
    properties(Access = private)
        % Panel   (uifigure) The main window of the dialog
        Panel
        TopGrid
        TopLeftGrid
        
        % UI components
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

        % BinRulesManager   (stats.internal.dfit.BinRulesManager)
        BinRulesManager
        
        % Listeners   (cell of listener) Listeners to various events
        Listeners = {};
        
        % SessionClearingListener (listener) Listener to
        % DistributionFitting's SessionClearing event
        SessionClearingListener
    end
    
    methods(Static)
        function this = getInstance()
            persistent Instance;
            if ~isa(Instance, 'stats.internal.dfit.DefaultBinRulesPanel') || ~isvalid(Instance)
                Instance = stats.internal.dfit.DefaultBinRulesPanel();
            end
            this = Instance;
        end
    end
    
    methods
        function showPanel(this)
            if ~isvalid(this.Panel) || isempty(this.Panel)
                this.createAndInitializeGUIStates();
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
        function this = DefaultBinRulesPanel()
            this.BinRulesManager = stats.internal.dfit.BinRulesManager.getInstance();
            this.createAndInitializeGUIStates();
            
            this.SessionClearingListener = event.listener(stats.internal.dfit.DistributionFitting.getInstance(), 'SessionClearing', @(~,~) this.delete);
        end
        
        function createAndInitializeGUIStates(this)
            this.createGUIComponents();
            this.initializeGUIStates();
            this.addListener(this.BinRulesManager, 'SaveBinRuleAsDefault', @this.initializeGUIStates);
        end
        
        function createGUIComponents(this)
            this.Panel = uifigure('Position', iGetDefaultPanelPosition(),...
                'Name', getString(message('stats:dfittool:title_defaultBinWidth')), ...
                'CloseRequestFcn', @(~,~) this.cancelButtonClickedCallback, ...,
                'Tag', 'dfDefaultBinRulesPanelUIFigure', ...
                'Visible', 'off');
            
            mainGrid = uigridlayout(this.Panel, [5, 1]);
            mainGrid.RowHeight = {'1x', 1, 'fit', 1, 'fit'};
            
            this.TopGrid = uigridlayout(mainGrid, [1, 2]);
            this.TopGrid.Padding = [0 0 0 0];
            this.TopGrid.RowHeight = {'fit'};
            this.TopGrid.ColumnWidth = {'fit','1x'};
            this.TopGrid.Layout.Row = 1;
            this.TopGrid.Layout.Column = 1;
            
            this.TopLeftGrid = uigridlayout(this.TopGrid, [7, 3]);
            this.TopLeftGrid.Padding = [0 0 0 0];
            this.TopLeftGrid.RowHeight = {18, 18, 18, 18, 18, 18, 18};
            this.TopLeftGrid.ColumnWidth = {24, 24, 'fit'};
            this.TopLeftGrid.Layout.Row = 1;
            this.TopLeftGrid.Layout.Column = 1;
            
            iCreateHorizontalLine(mainGrid, 2, 1);

            middleGrid = uigridlayout(mainGrid, [2, 1]);
            middleGrid.Padding = [0 0 0 0];
            middleGrid.RowHeight = {'fit', 'fit'};
            middleGrid.ColumnWidth = {'fit'};
            middleGrid.Layout.Row = 3;
            middleGrid.Layout.Column = 1;
            
            iCreateHorizontalLine(mainGrid, 4, 1);

            bottomGrid = uigridlayout(mainGrid, [1, 4]);
            bottomGrid.Padding = [0 0 0 0];
            bottomGrid.RowHeight = {'fit'};
            bottomGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit'};
            bottomGrid.Layout.Row = 5;
            bottomGrid.Layout.Column = 1;
            
            this.BinRulesButtonGroup = uibuttongroup(this.TopLeftGrid, ...
                'BorderType', 'none', ...
                'SelectionChangedFcn', @(~,~) this.radioButtonSelectionChangedCallback, ...
                'Tag', 'dfDefaultBinRulesPanelBinRulesButtonGroup');
            this.BinRulesButtonGroup.Layout.Row = [1 5];
            this.BinRulesButtonGroup.Layout.Column = 1;
            
            this.BinWidthButtonGroup = uibuttongroup(this.TopLeftGrid, ...
                'BorderType', 'none', ...
                'Tag', 'dfDefaultBinRulesPanelBinWidthButtonGroup');
            this.BinWidthButtonGroup.Layout.Row = [6 7];
            this.BinWidthButtonGroup.Layout.Column = 2;
            
            freedManDiaconisRuleLabel = uilabel(this.TopLeftGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_Freedman')));
            freedManDiaconisRuleLabel.Layout.Row = 1;
            freedManDiaconisRuleLabel.Layout.Column = [2 3];
            
            scottRuleLabel = uilabel(this.TopLeftGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_Scott')));
            scottRuleLabel.Layout.Row = 2;
            scottRuleLabel.Layout.Column = [2 3];
            
            numberOfBinsGrid = uigridlayout(this.TopLeftGrid, [1, 2]);
            numberOfBinsGrid.Padding = [0 0 0 0];
            numberOfBinsGrid.ColumnWidth = {'fit', '1x'};
            numberOfBinsGrid.RowHeight = {18};
            numberOfBinsGrid.Layout.Row = 3;
            numberOfBinsGrid.Layout.Column = [2 3];
            
            numberOfBinsLabel = uilabel(numberOfBinsGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_NumBins')));
            numberOfBinsLabel.Layout.Row = 1;
            numberOfBinsLabel.Layout.Column = 1;
            
            this.NumberOfBinsEditField = uieditfield(numberOfBinsGrid, ...
                'ValueChangingFcn', @(~,~) this.numberOfBinsEditFieldValueChangingCallback(), ...
                'Tag', 'dfDefaultBinRulesPanelNumberOfBinsEditField');
            this.NumberOfBinsEditField.Layout.Row = 1;
            this.NumberOfBinsEditField.Layout.Column = 2;
            
            binsCenteredOnIntegersLabel = uilabel(this.TopLeftGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_Integers')));
            binsCenteredOnIntegersLabel.Layout.Row = 4;
            binsCenteredOnIntegersLabel.Layout.Column = [2 3];
            
            binWidthGrid = uigridlayout(this.TopLeftGrid, [1, 2]);
            binWidthGrid.Padding = [0 0 0 0];
            binWidthGrid.ColumnWidth = {'fit', '1x'};
            binWidthGrid.RowHeight = {18};
            binWidthGrid.Layout.Row = 5;
            binWidthGrid.Layout.Column = [2 3];
            
            binWidthLabel = uilabel(binWidthGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_BinWidth')));
            binWidthLabel.Layout.Row = 1;
            binWidthLabel.Layout.Column = 1;
            
            this.BinWidthEditField = uieditfield(binWidthGrid, ...
                'valueChangingFcn', @(~,~) this.binWidthEditFieldValueChangingCallback(), ...
                'Tag', 'dfDefaultBinRulesPanelBinWidthEditField');
            this.BinWidthEditField.Layout.Row = 1;
            this.BinWidthEditField.Layout.Column = 2;
            
            this.AutoBinPlacementLabel = uilabel(this.TopLeftGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_BinWidthAny')), ...
                'Tag', 'dfDefaultBinRulesPanelAutoBinPlacementlabel');
            this.AutoBinPlacementLabel.Layout.Row = 6;
            this.AutoBinPlacementLabel.Layout.Column = 3;
            
            binBoundaryAtGrid = uigridlayout(this.TopLeftGrid, [1, 2]);
            binBoundaryAtGrid.Padding = [0 0 0 0];
            binBoundaryAtGrid.ColumnWidth = {'fit', '1x'};
            binBoundaryAtGrid.RowHeight = {18};
            binBoundaryAtGrid.Layout.Row = 7;
            binBoundaryAtGrid.Layout.Column = 3;
            
            this.BinBoundaryAtLabel = uilabel(binBoundaryAtGrid, ...
                'Text', getString(message('stats:dfittool:radiobutton_BinWidthSpecify')), ...
                'Tag', 'dfDefaultBinRulesPanelBinBoundaryAtLabel');
            this.BinBoundaryAtLabel.Layout.Row = 1;
            this.BinBoundaryAtLabel.Layout.Column = 1;
            
            this.BinBoundaryAtEditField = uieditfield(binBoundaryAtGrid, ...
                'ValueChangingFcn', @(~,~) this.binBoundaryAtEditFieldValueChangingCallback(), ...
                'Tag', 'dfDefaultBinRulesPanelBinBoundaryAtEditField');
            this.BinBoundaryAtEditField.Layout.Row = 1;
            this.BinBoundaryAtEditField.Layout.Column = 2;
            
            this.FreedmanDiaconisRuleRadioButton = uiradiobutton(this.BinRulesButtonGroup, 'Text', '', 'Position', [10 112 150 20], 'Tag', 'FreedmanDiaconis');
            this.ScottRuleRadioButton = uiradiobutton(this.BinRulesButtonGroup, 'Text', '', 'Position', [10 84 150 20], 'Tag', 'ScottRule');
            this.NumberOfBinsRadioButton = uiradiobutton(this.BinRulesButtonGroup, 'Text', '', 'Position', [10 56 150 20], 'Tag', 'NumberOfBins');
            this.BinsCenteredOnIntegersRadioButton = uiradiobutton(this.BinRulesButtonGroup, 'Text', '', 'Position', [10 28 15 20], 'Tag', 'BinsCenteredOnIntegers');
            this.BinWidthRadioButton = uiradiobutton(this.BinRulesButtonGroup, 'Text', '', 'Position', [10 0 15 20], 'Tag', 'BinWidth');
            
            this.AutoBinPlacementRadioButton = uiradiobutton(this.BinWidthButtonGroup, 'Text', '', 'Position', [10 28 15 20], 'Tag', 'BinWidthAutoBinPlacement');
            this.BinBoundaryAtRadioButton = uiradiobutton(this.BinWidthButtonGroup, 'Text', '', 'Position', [10 0 15 20], 'Tag', 'BinWidthBinBoundaryAt');
            
            this.ApplyToAllCheckBox = uicheckbox(middleGrid, ...
                'Text', getString(message('stats:dfittool:checkbox_applyToAll')), ...
                'Tag', 'dfDefaultBinRulesPanelApplyToAllCheckBox');
            this.ApplyToAllCheckBox.Layout.Row = 1;
            this.ApplyToAllCheckBox.Layout.Column = 1;
            
            iCreateDummyMinSizeButton(bottomGrid, 1, 1);
            iCreateDummyMinSizeButton(bottomGrid, 1, 3);
            iCreateDummyMinSizeButton(bottomGrid, 1, 4);
            
            helpButton = uibutton(bottomGrid, ...
                'Text', getString(message('stats:dfittool:button_help')), ...
                'ButtonPushedFcn', @(~,~) this.helpButtonClickedCallback, ...
                'Tag', 'dfDefaultBinRulesPanelHelpButton');
            helpButton.Layout.Row = 1;
            helpButton.Layout.Column = 1;
            
            oKButton = uibutton(bottomGrid, ...
                'Text', getString(message('stats:dfittool:button_OK')), ...
                'ButtonPushedFcn', @(~,~) this.oKButtonClickedCallback, ...
                'Tag', 'dfDefaultBinRulesPanelOKButton');
            oKButton.Layout.Row = 1;
            oKButton.Layout.Column = 3;
            
            cancelButton = uibutton(bottomGrid, ...
                'Text', getString(message('stats:dfittool:button_cancel')), ...
                'ButtonPushedFcn', @(~,~) this.cancelButtonClickedCallback, ...
                'Tag', 'dfDefaultBinRulesPanelCancelButton');
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 4;
            
            stats.internal.dfit.centergui(this.Panel);
            this.Panel.Visible = 'on';
        end
        
        function initializeGUIStates(this, ~, ~)
            binDlgInfo = dfgetset('binDlgInfo');
            binRuleIdx = binDlgInfo.rule;
            numBinValueExpr = binDlgInfo.nbinsExpr;
            binWidthValueExpr = binDlgInfo.widthExpr;
            binWidthIdx = binDlgInfo.placementRule;
            binBoundaryAtValueExpr = binDlgInfo.anchorExpr;
            isApplyBinRuleToAll = binDlgInfo.applyToAll;
            
            selectedBinRuleRadioButton = this.getBinRuleRadioButton(binRuleIdx);
            selectedBinWidthRadioButton = this.getBinWidthRadioButton(binWidthIdx);
            
            this.BinRulesButtonGroup.SelectedObject = selectedBinRuleRadioButton;
            this.BinWidthButtonGroup.SelectedObject = selectedBinWidthRadioButton;
            this.enableBinWidthRadioButton(isBinWidthRadioButton(selectedBinRuleRadioButton));
            this.NumberOfBinsEditField.Value = numBinValueExpr;
            this.BinWidthEditField.Value = binWidthValueExpr;
            this.BinBoundaryAtEditField.Value = binBoundaryAtValueExpr;
            
            this.ApplyToAllCheckBox.Value = isApplyBinRuleToAll; 
        end
        
        function enableBinWidthRadioButton(this, onoff)
            this.AutoBinPlacementRadioButton.Enable = onoff;
            this.BinBoundaryAtRadioButton.Enable = onoff;
            this.AutoBinPlacementLabel.Enable = onoff;
            this.BinBoundaryAtLabel.Enable = onoff;
            this.BinBoundaryAtEditField.Enable = onoff;
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
            
            stats.internal.dfit.setbinwidthrules([], ...
                binRuleIdx, numBinValue, ...
                binWidthValue, ...
                binWidthIdx, binBoundaryAtValue, ...
                applyToAll, true);
            
            if applyToAll
                this.BinRulesManager.applyBinRuleToAllDataSets();
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
position([3,4]) = [310 300];
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

function tf = isBinWidthRadioButton(radioButtonObj)
tf = strcmp(radioButtonObj.Tag, 'BinWidth');
end