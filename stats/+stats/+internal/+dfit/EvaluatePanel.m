classdef EvaluatePanel < handle
    % EvaluatePanel   View for the Evaluate dialog
    
    % Copyright 2019 The MathWorks, Inc
    
    properties(Access = private)
        % Panel   (uifigure) The main window of the dialog
        Panel
        
        % Various UI components
        FitListBox
        FunctionDropDown
        VectorTextFieldLabel
        VectorTextField
        ConfidenceBoundsCheckBox
        ConfidenceBoundsLevelLabel
        ConfidenceBoundsTextField
        ConfidenceBoundsPercentLabel
        PlotFunctionCheckBox
        ApplyButton
        ExportToWorkspaceButton
        CloseButton
        HelpButton
        ClickApplyMessageLabel
        
        % EvaluateResultsTable   (uitable)
        EvaluateResultsTable
        
        FunctionOptionTags = {...
            'pdf', ...
            'cdf', ...
            'icdf', ...
            'survivor', ...
            'cumhazard', ...
            'hazrate'};
        
        FunctionOptions = {...
            'stats:dfittool:choice_pdf', ...
            'stats:dfittool:choice_cdf', ...
            'stats:dfittool:choice_invCDF', ...
            'stats:dfittool:choice_survivor', ...
            'stats:dfittool:choice_cumHaz', ...
            'stats:dfittool:choice_hazRate'};
    end
    
    properties(Dependent, SetAccess = private)
        FitList
        SelectedFits
        SelectedFunction
        VectorValue
        IsComputeConfidenceBounds
        LevelValue
        IsPlotFunction
    end
    
    events        
        ApplyButtonClicked
        ExportToWorkspaceButtonClicked
        CloseButtonClicked
        HelpButtonClicked
        
        FunctionDropDownValueChanged
        PlotFunctionCheckBoxValueChanged
        ParameterValuesChanged
    end
    
    methods
        function this = EvaluatePanel(fitList, functionType, vectorLabel, vectorValue, isConfidenceBoundEnabled, evaluateTable)
            this.createGUIComponents();
            
            this.FitListBox.Items = fitList;
            
            this.FunctionDropDown.Value = (functionType);
            
            this.updateVectorTextFieldLabel(vectorLabel);
            this.VectorTextField.Value = vectorValue;
            
            this.updateConfidenceBoundComponentStates(isConfidenceBoundEnabled);
            
            this.ConfidenceBoundsCheckBox.Value = false;
            this.ConfidenceBoundsTextField.Value = iDefaultConfidenceBoundValue();
            
            this.PlotFunctionCheckBox.Value = false;
            
            this.updateEvaluateResultsTable(evaluateTable);
        end
        
        function showPanel(this)
            figure(this.Panel)
        end
        
        function showMessageDialog(this, messageStr)
            uialert(this.Panel,messageStr,...
                getString(message('stats:dfittool:title_evaluateError')));
        end
        
        function delete(this)
            delete(this.Panel);
        end
        
        % update UI component states
        function updateFitListBox(this, fitList, selectedFits)
            this.FitListBox.Items = fitList;
            this.FitListBox.Value = selectedFits;
        end
        
        function updateVectorTextFieldLabel(this, label)
            this.VectorTextFieldLabel.Text = label;
        end
        
        function updateConfidenceBoundComponentStates(this, isEnabled)
            if isEnabled
                this.ConfidenceBoundsCheckBox.Enable = 'on';
                this.ConfidenceBoundsTextField.Enable = 'on';
                this.ConfidenceBoundsLevelLabel.Enable = 'on';
                this.ConfidenceBoundsPercentLabel.Enable = 'on';
            else
                this.ConfidenceBoundsCheckBox.Enable = 'off';
                this.ConfidenceBoundsTextField.Enable = 'off';
                this.ConfidenceBoundsLevelLabel.Enable = 'off';
                this.ConfidenceBoundsPercentLabel.Enable = 'off';
            end
        end
        
        function updatePlotCheckBox(this, isChecked)
            this.PlotFunctionCheckBox.Value = isChecked;
        end
        
        function updateEvaluateResultsTable(this, tableData)
            this.EvaluateResultsTable.Data = tableData;
            this.showEvaluateResultsTable(~isempty(tableData));
        end
        
        % Setter/Getter methods
        function value = get.FitList(this)
            value = this.FitListBox.Items;
        end
        
        function value = get.SelectedFits(this)
            value = this.FitListBox.Value;
        end
        
        function value = get.SelectedFunction(this)
            value = this.FunctionDropDown.Value;
        end
        
        function value = get.VectorValue(this)
            value = this.VectorTextField.Value;
        end
        
        function tf = get.IsComputeConfidenceBounds(this)
            tf = this.ConfidenceBoundsCheckBox.Value;
        end
        
        function value = get.LevelValue(this)
            value = this.ConfidenceBoundsTextField.Value;
        end
        
        function tf = get.IsPlotFunction(this)
            tf = this.PlotFunctionCheckBox.Value;
        end
    end
    
    methods(Access = private)        
        function createGUIComponents(this)
            this.Panel = uifigure( ...
                'Name', getString(message('stats:dfittool:title_evaluate')), ...
                'Position', iGetDefaultPanelPosition(), ...
                'CloseRequestFcn', @(~,~) this.closeButtonClickedCallback, ...
                'Tag', 'dfEvaluatePanelUIFigure', ...
                'Visible', 'off');
            
            this.createGridLayout();
            this.showEvaluateResultsTable(false);
            
            stats.internal.dfit.centergui(this.Panel);
        end
        
        function showEvaluateResultsTable(this, isShowTable)
            if isShowTable
                this.EvaluateResultsTable.Visible = 'on';
                this.ClickApplyMessageLabel.Visible = 'off';
                this.ExportToWorkspaceButton.Enable = 'on';
            else
                this.EvaluateResultsTable.Visible = 'off';
                this.ClickApplyMessageLabel.Visible = 'on';
                this.ExportToWorkspaceButton.Enable = 'off';
            end
        end
        
        % Grid layout
        function createGridLayout(this)
            mainGrid = uigridlayout(this.Panel, [3, 1]);
            mainGrid.RowHeight = {'1x', 1, 'fit'};
            mainGrid.ColumnWidth = {'1x'};
 
            this.createTopPanels(mainGrid);
            this.createBottomPanel(mainGrid);
        end
        
        function createTopPanels(this, parentGrid)
            topGrid = uigridlayout(parentGrid, [1 3]);
            topGrid.Padding = [0 0 0 0];
            topGrid.RowHeight = {'1x'};
            topGrid.ColumnWidth = {'fit', 1, '1x'};
            topGrid.Layout.Row = 1;
            topGrid.Layout.Column = 1;
            
            p1 = uipanel(topGrid);
            p1.BackgroundColor = [.80 .80 .80];
            p1.Layout.Row = 1;
            p1.Layout.Column = 2;
            
            p1 = uipanel(parentGrid);
            p1.BackgroundColor = [.80 .80 .80];
            p1.Layout.Row = 2;
            p1.Layout.Column = 1;
            
            this.createTopLeftPanel(topGrid);
            this.createTopRightPanel(topGrid);
        end
        
        function createTopLeftPanel(this, parentGrid)
            topLeftGrid = uigridlayout(parentGrid, [8 1]);
            topLeftGrid.Padding = [0 0 0 0];
            topLeftGrid.ColumnWidth = {250};
            topLeftGrid.RowHeight = {'fit', '1x', 'fit', 'fit', ...
                'fit', 'fit', 1, 'fit'};
            topLeftGrid.Layout.Row = 1;
            topLeftGrid.Layout.Column = 1;
            
            funAtGrid = uigridlayout(topLeftGrid, [2 2]);
            funAtGrid.RowHeight = {'fit', 'fit'};
            funAtGrid.ColumnWidth = {'fit', '1x'};
            funAtGrid.Layout.Row = 3;
            funAtGrid.Layout.Column = 1;
            
            levelGrid = uigridlayout(topLeftGrid, [1, 3]);
            levelGrid.RowHeight = {'fit'};
            levelGrid.ColumnWidth = {'fit', 60, 'fit'};
            levelGrid.Layout.Row = 5;
            levelGrid.Layout.Column = 1;
            
            applyGrid = uigridlayout(topLeftGrid, [1, 2]);
            applyGrid.RowHeight = {'fit'};
            applyGrid.ColumnWidth = {'1x', 'fit'};
            applyGrid.Layout.Row = 6;
            applyGrid.Layout.Column = 1;
            
            p1 = uipanel(topLeftGrid);
            p1.BackgroundColor = [.80 .80 .80];
            p1.Layout.Row = 7;
            p1.Layout.Column = 1;
            
            this.createFitListBox(topLeftGrid);
            this.createConfidenceBoundsCheckBox(topLeftGrid);
            this.createFunctionDropDown(funAtGrid);
            this.createVectorTextField(funAtGrid);
            this.createConfidenceBoundLevelTextField(levelGrid);
            
            this.createApplyButton(applyGrid);
            this.createPlotFunctionCheckBox(topLeftGrid);
        end
        
        function createTopRightPanel(this, parentGrid)
            topRightGrid = uigridlayout(parentGrid, [3 3]);
            topRightGrid.ColumnWidth = {'1x', 'fit', '1x'};
            topRightGrid.RowHeight = {'1x', 'fit', '1x'};
            topRightGrid.Layout.Row = 1;
            topRightGrid.Layout.Column = 3;
            
            this.createClickApplyMessageLabel(topRightGrid);
            this.createEvaluateResultsTable(topRightGrid);
        end
        
        function createBottomPanel(this, parentGrid)
            bottomGrid = uigridlayout(parentGrid, [1 4]);
            bottomGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit'};
            bottomGrid.RowHeight = {'fit'};
            bottomGrid.Padding = [0 0 0 0];
            bottomGrid.Layout.Row = 3;
            bottomGrid.Layout.Column = 1;
            
            this.createExportToWorkspaceButton(bottomGrid, 3);
            this.createHelpButton(bottomGrid, 1);
            this.createCloseButton(bottomGrid, 4);
        end
        
        % UI components
        function createFitListBox(this, parentGrid)
            l = uilabel(parentGrid, ...
                'Text', getString(message('stats:dfittool:label_selectFit')));
            l.Layout.Row = 1;
            l.Layout.Column = 1;
            
            fitNames = stats.internal.dfit.getDataBaseNames(stats.internal.dfit.getfitdb);
            this.FitListBox = uilistbox(parentGrid, ...
                'Items', fitNames, ...
                'Tag', 'dfEvaluatePanelFitListBox');
            this.FitListBox.Multiselect = 'on';
            this.FitListBox.Layout.Row = 2;
            this.FitListBox.Layout.Column = 1;
        end
        
        function createFunctionDropDown(this, parentGrid)
            l = uilabel(parentGrid,...
                'Text', iGetMessageString({'stats:dfittool:label_function'}));
            l.Layout.Row = 1;
            l.Layout.Column = 1;
            
            this.FunctionDropDown = uidropdown(parentGrid,...
                'Items', iGetMessageString(this.FunctionOptions),...
                'ItemsData', this.FunctionOptionTags, ...
                'ValueChangedFcn', @(~,~) this.functionDropDownValueChanged, ...
                'Tag', 'dfEvaluatePanelFunctionDropDown');
            this.FunctionDropDown.Layout.Row = 1;
            this.FunctionDropDown.Layout.Column = 2;
        end
        
        function createVectorTextField(this, parentGrid)
            this.VectorTextFieldLabel = uilabel(parentGrid, ...
                'Tag', 'dfEvaluatePanelVectorTextFieldLabel');
            this.VectorTextFieldLabel.Layout.Row = 2;
            this.VectorTextFieldLabel.Layout.Column = 1;
            
            this.VectorTextField = uieditfield(parentGrid, ...
                'ValueChangedFcn', @(~,~) this.parameterValuesChanged, ...
                'Tag', 'dfEvaluatePanelVectorTextField');
            this.VectorTextField.Layout.Row = 2;
            this.VectorTextField.Layout.Column = 2;
        end
        
        function createConfidenceBoundsCheckBox(this, parentGrid)
            this.ConfidenceBoundsCheckBox = uicheckbox(parentGrid,...
                'Text', getString(message('stats:dfittool:checkbox_showConfBnds')), ...
                'ValueChangedFcn', @(~,~) this.parameterValuesChanged, ...
                'Tag', 'dfEvaluatePanelConfidenceBoundsCheckBox');
            this.ConfidenceBoundsCheckBox.Layout.Row = 4;
            this.ConfidenceBoundsCheckBox.Layout.Column = 1;
        end
        
        function createConfidenceBoundLevelTextField(this, parentGrid)
            this.ConfidenceBoundsLevelLabel = uilabel(parentGrid, ...
                'Text', ['    ', getString(message('stats:dfittool:label_level'))], ...
                'Tag', 'dfEvaluatePanelConfidenceBoundsLevelLabel');
            this.ConfidenceBoundsLevelLabel.Layout.Row = 1;
            this.ConfidenceBoundsLevelLabel.Layout.Column = 1;
            
            this.ConfidenceBoundsTextField = uieditfield(parentGrid, ...
                'numeric', ...
                'Value', iDefaultConfidenceBoundValue(), ...
                'ValueChangedFcn', @(~,~) this.parameterValuesChanged, ...
                'Tag', 'dfEvaluatePanelConfidenceBoundsTextField');
            this.ConfidenceBoundsTextField.Limits = [0, 100];
            this.ConfidenceBoundsTextField.Layout.Row = 1;
            this.ConfidenceBoundsTextField.Layout.Column = 2;
            
            this.ConfidenceBoundsPercentLabel = uilabel(parentGrid, ...
                'Text', '%', ...
                'Tag', 'dfEvaluatePanelConfidenceBoundsPercentLabel');
            this.ConfidenceBoundsPercentLabel.Layout.Row = 1;
            this.ConfidenceBoundsPercentLabel.Layout.Column = 3;
        end
        
        function createApplyButton(this, parentGrid)
            iCreateDummyMinSizeButton(parentGrid, 1, 2);
            
            this.ApplyButton = uibutton(parentGrid,...
                'Text', getString(message('stats:dfittool:button_apply')), ...
                'ButtonPushedFcn', @(~, ~) this.applyButtonClickedCallback(), ...
                'Tag', 'dfEvaluatePanelApplyButton');
            this.ApplyButton.Layout.Row = 1;
            this.ApplyButton.Layout.Column = 2;
        end
        
        function createPlotFunctionCheckBox(this, parentGrid)
            this.PlotFunctionCheckBox = uicheckbox(parentGrid,...
                'Text', getString(message('stats:dfittool:checkbox_plotFunction')), ...
                'ValueChangedFcn', @this.plotFunctionCheckBoxValueChangedCallback, ...
                'Tag', 'dfEvaluatePanelPlotFunctionCheckBox');
            this.PlotFunctionCheckBox.Layout.Row = 8;
            this.PlotFunctionCheckBox.Layout.Column = 1;
        end
        
        function createHelpButton(this, parentGrid, column)
            iCreateDummyMinSizeButton(parentGrid, 1, column);

            this.HelpButton = uibutton(parentGrid,...
                'Text', getString(message('stats:dfittool:button_help')), ...
                'ButtonPushedFcn', @(~, ~)this.helpButtonClickedCallback(), ...
                'Tag', 'dfEvaluatePanelHelpButton');
            this.HelpButton.Layout.Row = 1;
            this.HelpButton.Layout.Column = column;
        end
        
        function createExportToWorkspaceButton(this, parentGrid, column)
            iCreateDummyMinSizeButton(parentGrid, 1, column);

            this.ExportToWorkspaceButton = uibutton(parentGrid,...
                'Text', getString(message('stats:dfittool:button_exportToWS')), ...
                'ButtonPushedFcn', @(~, ~) this.exportToWorkspaceButtonClickedCallback(), ...
                'Tag', 'dfEvaluatePanelExportToWorkspaceButton');
            this.ExportToWorkspaceButton.Layout.Row = 1;
            this.ExportToWorkspaceButton.Layout.Column = column;
        end
        
        function createCloseButton(this, parentGrid, column)
            iCreateDummyMinSizeButton(parentGrid, 1, column);

            this.CloseButton = uibutton(parentGrid, 'Text', ...
                getString(message('stats:dfittool:button_close')), ...
                'ButtonPushedFcn', @(~, ~)this.closeButtonClickedCallback(), ...
                'Tag', 'dfEvaluatePanelCloseButton');
            this.CloseButton.Layout.Row = 1;
            this.CloseButton.Layout.Column = column;
        end
       
        function createClickApplyMessageLabel(this, parentGrid)
            this.ClickApplyMessageLabel = uilabel(parentGrid, ...
                'Text', iGetMessageString({'stats:dfittool:label_evaluatePressApply'}), ...
                'Tag', 'dfEvaluatePanelClickApplyMessageLabel');
            this.ClickApplyMessageLabel.Layout.Row = 2;
            this.ClickApplyMessageLabel.Layout.Column = 2;
        end
        
        function createEvaluateResultsTable(this, parentGrid)
            this.EvaluateResultsTable = uitable(parentGrid, ...
                'RowStriping', 'off', ...
                'Tag', 'dfEvaluatePanelEvaluateResultsTable');
            this.EvaluateResultsTable.RowName = [];
            this.EvaluateResultsTable.Layout.Row = [1 3];
            this.EvaluateResultsTable.Layout.Column = [1 3];
        end
        
        % Callbacks
        function functionDropDownValueChanged(this, ~, ~)
            this.notify('FunctionDropDownValueChanged');
        end
        
        function parameterValuesChanged(this, ~, ~)
            this.notify('ParameterValuesChanged');
        end
        
        function applyButtonClickedCallback(this, ~, ~)
            this.notify('ApplyButtonClicked');
        end
        
        function exportToWorkspaceButtonClickedCallback(this, ~, ~)
            this.notify('ExportToWorkspaceButtonClicked');
        end
        
        function closeButtonClickedCallback(this, ~, ~)
            this.notify('CloseButtonClicked');
            delete(this);
        end
        
        function helpButtonClickedCallback(this, ~, ~)
            this.notify('HelpButtonClicked');
        end
        
        function plotFunctionCheckBoxValueChangedCallback(this, ~, ~)
            this.notify('PlotFunctionCheckBoxValueChanged');
        end
    end
end

% helpers
function position = iGetDefaultPanelPosition()
position = get(0,'defaultfigureposition');
position([3,4]) = [630 470];
end

function cellStr = iGetMessageString(cellStr)
cellStr = cellfun(@(x) getString(message(x)), cellStr, 'UniformOutput', false);
end

function value = iDefaultConfidenceBoundValue()
value = 95;
end

function iCreateDummyMinSizeButton(grid, row, column)
dummyText = repmat('X',1,8);
minButton = uibutton(grid, 'Text', dummyText);
minButton.Layout.Row = row;
minButton.Layout.Column = column;
minButton.Visible = 'off';
end