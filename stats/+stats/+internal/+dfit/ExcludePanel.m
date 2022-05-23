classdef ExcludePanel < handle
    % ExcludePanel   View for the Exclude dialog
    
    % Copyright 2019 The MathWorks, Inc
    
    properties(Access = private)
        % Panel   (uifigure) The main window of the dialog
        Panel
        
        % UI Components
        ExclusionRuleNameTextField
        LowerLimitDropDown
        UpperLimitDropDown
        LowerLimitTextField
        UpperLimitTextField
        SelectDataDropDown
        ExcludeGraphicallyButton
        CreateExclusionRuleButton
        ExistingExclusionRulesListBox
        CopyButton
        ViewButton
        RenameButton
        DeleteButton
        CloseButton
        HelpButton
    end
    
    properties(Dependent, SetAccess = private)
        % SelectedData   (string) name of data set selected in the drop
        % down menu
        SelectedData
        
        % SelectedExclusionRules   (cell of strings) names of exclusion
        % rules selected in the Exclusion Rules list box
        SelectedExclusionRules
        
        % ExclusionRuleName   (string) value in ExclusiongRuleNameTextField
        % component
        ExclusionRuleName
        
        % IsLowerLimitLessEqual   (boolean) true if "<=" is selected in the
        % drop down menu, false if "<" is selected
        IsLowerLimitLessEqual
        
        % IsUpperLimitGreaterEqual   (boolean) true if ">=" is selected in
        % the drop down menu, false if ">" is selected
        IsUpperLimitGreaterEqual
        
        % LowerLimitTextFieldValue   (string) value in LowerLimitTextField
        % component
        LowerLimitTextFieldValue
        
        % UpperLimitTextFieldValue   (string) value in UpperLimitTextField
        % component
        UpperLimitTextFieldValue
    end
    
    events
        LowerLimitTextFieldValueChanged
        UpperLimitTextFieldValueChanged
        SelectDataDropDownValueChanged
        ExcludeGraphicallyButtonClicked
        CreateExclusionRuleButtonClicked
        ExistingExclusionRulesListBoxValueChanged
        CopyButtonClicked
        ViewButtonClicked
        RenameButtonClicked
        DeleteButtonClicked
        CloseButtonClicked
        HelpButtonClicked
    end
    
    methods
        function this = ExcludePanel()
            this.createGUIComponents();
        end
        
        function showPanel(this)
            figure(this.Panel);
        end
        
        function hidePanel(this)
            this.Panel.Visible = 'off';
        end
        
        function delete(this)
            delete(this.Panel);
        end
        
        function resetLeftPanel(this)
            this.ExclusionRuleNameTextField.Value = "";
            this.LowerLimitDropDown.Value = this.LowerLimitDropDown.Items(1);
            this.UpperLimitDropDown.Value = this.UpperLimitDropDown.Items(1);
            this.LowerLimitTextField.Value = "";
            this.UpperLimitTextField.Value = "";
            this.SelectDataDropDown.Value = this.SelectDataDropDown.Items(1);
            this.ExcludeGraphicallyButton.Enable = false;
        end
        
        function showMessageDialog(this, messageStr, titleStr)
            uialert(this.Panel, messageStr, titleStr);
        end
        
        function updateExclusionRuleNameTextFieldValue(this, value)
            this.ExclusionRuleNameTextField.Value = value;
        end
        
        function updateLowerLimitDropDownValue(this, hasEqualSign)
            if hasEqualSign
                idx = 1;
            else
                idx = 2;
            end
            this.LowerLimitDropDown.Value = this.LowerLimitDropDown.Items(idx);
        end
        
        function updateUpperLimitDropDownValue(this, hasEqualSign)
            if hasEqualSign
                idx = 1;
            else
                idx = 2;
            end
            this.UpperLimitDropDown.Value = this.UpperLimitDropDown.Items(idx);
        end
        
        function updateLowerLimitTextFieldValue(this, value)
            this.LowerLimitTextField.Value = value;
        end
        
        function updateUpperLimitTextFieldValue(this, value)
            this.UpperLimitTextField.Value = value;
        end
        
        function updateSelectDataDropDownList(this, dataSetList)
            this.SelectDataDropDown.Items = dataSetList;
        end
        
        function updateSelectDataDropDownValue(this, value)
            if isempty(value)
                this.SelectDataDropDown.Value = this.SelectDataDropDown.Items(1);
            else
                this.SelectDataDropDown.Value = value;
            end
        end
        
        function updateExistingExclusionRulesList(this, exclusionRulesList)
            this.ExistingExclusionRulesListBox.Items = exclusionRulesList;
        end
        
        function updateExistingExclusionRulesValue(this, exclusionRulesValue)
            this.ExistingExclusionRulesListBox.Value = exclusionRulesValue;
        end
        
        function enableExcludeGraphicallyButton(this, tf)
            this.ExcludeGraphicallyButton.Enable = tf;
        end
        
        function enableCopyButton(this, tf)
            this.CopyButton.Enable = tf;
        end
        
        function enableViewButton(this, tf)
            this.ViewButton.Enable = tf;
        end
        
        function enableRenameButton(this, tf)
            this.RenameButton.Enable = tf;
        end
        
        function enableDeleteButton(this, tf)
            this.DeleteButton.Enable = tf;
        end
        
        % Getter/setter methods
        function value = get.LowerLimitTextFieldValue(this)
            value = this.LowerLimitTextField.Value;
        end
        
        function value = get.UpperLimitTextFieldValue(this)
            value = this.UpperLimitTextField.Value;
        end
        
        function value = get.SelectedData(this)
            value = this.SelectDataDropDown.Value;
        end
        
        function value = get.SelectedExclusionRules(this)
            value = this.ExistingExclusionRulesListBox.Value;
        end
        
        function value = get.ExclusionRuleName(this)
            value = this.ExclusionRuleNameTextField.Value;
        end
        
        function tf = get.IsLowerLimitLessEqual(this)
            tf = isequal(this.LowerLimitDropDown.Value, '<=');
        end
        
        function tf = get.IsUpperLimitGreaterEqual(this)
            tf = isequal(this.UpperLimitDropDown.Value, '>=');
        end
    end
    
    methods(Access = private)
        function createGUIComponents(this)
            this.Panel = uifigure( ...
                'Name', getString(message('stats:dfittool:title_exclude')), ...
                'Position', iGetDefaultPanelPosition(), ...
                'CloseRequestFcn', @(~, ~)this.closeButtonClickedCallback(), ...
                'Tag', 'dfExcludePanelUIFigure', ...
                'Visible', 'off');
            
            this.createGridLayout();
            
            stats.internal.dfit.centergui(this.Panel);
        end
        
        % Grid layout
        function createGridLayout(this)
            mainGrid = uigridlayout(this.Panel, [3, 1]);
            mainGrid.RowHeight = {'1x', 1, 'fit'};
            mainGrid.ColumnWidth = {'1x'};
            
            this.createTopGrid(mainGrid);
            
            horizontalLine = uipanel(mainGrid);
            horizontalLine.BackgroundColor = [.80 .80 .80];
            horizontalLine.Layout.Row = 2;
            horizontalLine.Layout.Column = 1;
            
            this.createBottomGrid(mainGrid);
        end
        
        function createTopGrid(this, parentGrid)
            topGrid = uigridlayout(parentGrid, [1, 3]);
            topGrid.Padding = [0 0 0 10];
            topGrid.RowHeight = {'1x'};
            topGrid.ColumnWidth = {'fit', 1, 'fit'};
            topGrid.Layout.Row = 1;
            topGrid.Layout.Column = 1;
            
            this.createTopLeftGrid(topGrid);
            
            verticalLine = uipanel(topGrid);
            verticalLine.BackgroundColor = [.80 .80 .80];
            verticalLine.Layout.Row = 1;
            verticalLine.Layout.Column = 2;
            
            this.createTopRightGrid(topGrid);
        end
        
        function createTopLeftGrid(this, parentGrid)
            topLeftGrid = uigridlayout(parentGrid, [5, 2]);
            topLeftGrid.Padding = [0 0 0 0];
            topLeftGrid.RowHeight = {'fit', 100, 70, 1, 'fit'};
            topLeftGrid.ColumnWidth = {'fit', '1x'};
            topLeftGrid.Layout.Row = 1;
            topLeftGrid.Layout.Column = 1;
            
            this.createExclusionRuleNameGrid(topLeftGrid);
            this.createExcludeSectionsGrid(topLeftGrid);
            this.createExcludeGraphicallyGrid(topLeftGrid);
            this.createCreateExclusionRuleButtonGrid(topLeftGrid);
        end
        
        function createTopRightGrid(this, parentGrid)
            topRightGrid = uigridlayout(parentGrid);
            topRightGrid.Padding = [0, 0, 0, 0];
            topRightGrid.RowHeight = {'fit', 178, 'fit'};
            topRightGrid.ColumnWidth = {'fit'};
            topRightGrid.Layout.Row = 1;
            topRightGrid.Layout.Column = 3;
            
            existingExclusionRulesLabel = uilabel(topRightGrid,...
                'Text', [getString(message('stats:dfittool:label_existingExclusionRules')) '     ']);
            existingExclusionRulesLabel.Layout.Row = 1;
            existingExclusionRulesLabel.Layout.Column = 1;
            
            this.createExistingExclusionRulesListBox(topRightGrid);
            this.createCopyViewRenameDeleteButtonsGrid(topRightGrid);
            
        end

        function createBottomGrid(this, parentGrid)
            bottomGrid = uigridlayout(parentGrid, [1, 3]);
            bottomGrid.Padding = [0, 0, 0, 0];
            bottomGrid.RowHeight = {'fit'};
            bottomGrid.ColumnWidth = {'fit', '1x', 'fit'};
            bottomGrid.Layout.Row = 3;
            bottomGrid.Layout.Column = 1;
            
            this.createHelpButton(bottomGrid, 1);
            this.createCloseButton(bottomGrid, 3);
        end
        
        function createExclusionRuleNameGrid(this, parentGrid)
            exclusionRuleNameLabel = uilabel(parentGrid, 'Text', getString(message('stats:dfittool:label_exclusionRuleName')));
            exclusionRuleNameLabel.Layout.Row = 1;
            exclusionRuleNameLabel.Layout.Column = 1;
            this.createExclusionRuleNameTextField(parentGrid);
        end
        
        function createExcludeSectionsGrid(this, parentGrid)
            p = uipanel(parentGrid,...
                'Title', getString(message('stats:dfittool:label_excludeSections')));
            p.Layout.Row = 2;
            p.Layout.Column = [1 2];
            
            excludeSectionsGrid = uigridlayout(p, [2, 3]);
            excludeSectionsGrid.RowHeight = {'fit', 'fit'};
            excludeSectionsGrid.ColumnWidth = {'fit', 'fit', '1x'};
            
            excludeSectionsLabel = uilabel(excludeSectionsGrid, 'Text', getString(message('stats:dfittool:label_llim')));
            excludeSectionsLabel.Layout.Row = 1;
            excludeSectionsLabel.Layout.Column = 1;
            this.createLowerLimitTextField(excludeSectionsGrid);
            this.createLowerLimitDropDown(excludeSectionsGrid);
            
            excludeSectionsLabel = uilabel(excludeSectionsGrid, 'Text', getString(message('stats:dfittool:label_ulim')));
            excludeSectionsLabel.Layout.Row = 2;
            excludeSectionsLabel.Layout.Column = 1;
            this.createUpperLimitTextField(excludeSectionsGrid);
            this.createUpperLimitDropDown(excludeSectionsGrid);
        end
        
        function createExcludeGraphicallyGrid(this, parentGrid)
            p = uipanel(parentGrid,...
                'Title', getString(message('stats:dfittool:label_excludeGraphically')));
            p.Layout.Row = 3;
            p.Layout.Column = [1 2];
            
            excludeGraphicallyGrid = uigridlayout(p, [1, 3]);
            excludeGraphicallyGrid.RowHeight = {'fit', 'fit'};
            excludeGraphicallyGrid.ColumnWidth = {'fit', '1x', 'fit'};
            
            excludeGraphicallyLabel = uilabel(excludeGraphicallyGrid, 'Text', getString(message('stats:dfittool:label_selectDataSet')));
            excludeGraphicallyLabel.Layout.Row = 1;
            excludeGraphicallyLabel.Layout.Column = 1;
            this.createSelectDataDropDown(excludeGraphicallyGrid);
            this.createExcludeGraphicallyButton(excludeGraphicallyGrid);
        end
        
        function createCreateExclusionRuleButtonGrid(this, parentGrid)
            createExclusionRuleButtonGrid = uigridlayout(parentGrid, [1 2]);
            createExclusionRuleButtonGrid.RowHeight = {'fit'};
            createExclusionRuleButtonGrid.ColumnWidth = {'1x', 'fit'};
            createExclusionRuleButtonGrid.Layout.Row = 5;
            createExclusionRuleButtonGrid.Layout.Column = [1 2];
            
            this.createCreateExclusionRuleButton(createExclusionRuleButtonGrid);
        end
        
        function createCopyViewRenameDeleteButtonsGrid(this, parentGrid)
            copyViewRenameDeleteGrid = uigridlayout(parentGrid, [2 2]);
            copyViewRenameDeleteGrid.Padding = [0 0 0 0];
            copyViewRenameDeleteGrid.RowHeight = {'fit', 'fit'};
            copyViewRenameDeleteGrid.ColumnWidth = {'fit', 'fit'};
            copyViewRenameDeleteGrid.Layout.Row = 3;
            copyViewRenameDeleteGrid.Layout.Column = 1;
        
            this.createCopyButton(copyViewRenameDeleteGrid);
            this.createViewButton(copyViewRenameDeleteGrid);
            this.createRenameButton(copyViewRenameDeleteGrid);
            this.createDeleteButton(copyViewRenameDeleteGrid);
        end
        
        % UI Components
        function createExclusionRuleNameTextField(this, parentGrid)
            this.ExclusionRuleNameTextField = uieditfield(parentGrid, ...
                'Tag', 'dfExcludePanelExclusionRuleNameTextField');
            this.ExclusionRuleNameTextField.Layout.Row = 1;
            this.ExclusionRuleNameTextField.Layout.Column = 2;
        end
        
        function createLowerLimitDropDown(this, parentGrid)
            this.LowerLimitDropDown = uidropdown(parentGrid, ...
                'Items', {'<=', '<'},...
                'Tag', 'dfExcludePanelLowerLimitDropDown');
            this.LowerLimitDropDown.Layout.Row = 1;
            this.LowerLimitDropDown.Layout.Column = 2;
        end
        
        function createUpperLimitDropDown(this, parentGrid)
            this.UpperLimitDropDown = uidropdown(parentGrid, ...
                'Items', {'>=', '>'},...
                'Tag', 'dfExcludePanelUpperLimitDropDown');
            this.UpperLimitDropDown.Layout.Row = 2;
            this.UpperLimitDropDown.Layout.Column = 2;
        end
           
        function createLowerLimitTextField(this, parentGrid)
            this.LowerLimitTextField = uieditfield(parentGrid, ...
                'ValueChangedFcn', @(~, ~) this.lowerLimitTextFieldValueChangedCallback, ...
                'Tag', 'dfExcludePanelLowerLimitTextField');
            this.LowerLimitTextField.Layout.Row = 1;
            this.LowerLimitTextField.Layout.Column = 3;
        end
        
        function createUpperLimitTextField(this, parentGrid)
            this.UpperLimitTextField = uieditfield(parentGrid, ...
                'ValueChangedFcn', @(~, ~) this.upperLimitTextFieldValueChangedCallback, ...
                'Tag', 'dfExcludePanelUpperLimitTextField');
            this.UpperLimitTextField.Layout.Row = 2;
            this.UpperLimitTextField.Layout.Column = 3;
        end
        
        function createSelectDataDropDown(this, parentGrid)
            this.SelectDataDropDown = uidropdown(parentGrid, ...
                'Items', {}, ...
                'ValueChangedFcn', @(~, ~) this.selectDataDropDownValueChanged(), ...
                'Tag', 'dfExcludePanelSelectDataDropDown');
            this.SelectDataDropDown.Layout.Row = 1;
            this.SelectDataDropDown.Layout.Column = 2;
        end
        
        function createExcludeGraphicallyButton(this, parentGrid)
            iCreateDummyMinSizeButton(parentGrid, 1, 3);
            
            this.ExcludeGraphicallyButton = uibutton(parentGrid, ...
                'Text', getString(message('stats:dfittool:button_exclude')), ...
                'ButtonPushedFcn', @(~,~) this.excludeGraphicallyButtonClickedCallback(), ...
                'Tag', 'dfExcludePanelExcludeGraphicallyButton');
            this.ExcludeGraphicallyButton.Layout.Row = 1;
            this.ExcludeGraphicallyButton.Layout.Column = 3;
        end
        
        function createCreateExclusionRuleButton(this, parentGrid)
            iCreateDummyMinSizeButton(parentGrid, 1, 2);

            this.CreateExclusionRuleButton = uibutton(parentGrid,...
                'Text', getString(message('stats:dfittool:button_createExclusionRule')),...
                'ButtonPushedFcn', @(~,~) this.createExclusionRuleButtonClickedCallback(), ...
                'Tag', 'dfExcludePanelCreateExclusionRuleButton');
            this.CreateExclusionRuleButton.Layout.Row = 1;
            this.CreateExclusionRuleButton.Layout.Column = 2;
        end
        
        function createExistingExclusionRulesListBox(this, parentGrid)
            this.ExistingExclusionRulesListBox = uilistbox(parentGrid, ...
                'Multiselect', 'on', ...
                'Items', {},...
                'ValueChangedFcn', @(~,~) this.existingExclusionRulesListBoxValueChangedCallback(), ...
                'Tag', 'dfExcludePanelExistingExclusionRulesListBox');
            this.ExistingExclusionRulesListBox .Layout.Row = 2;
            this.ExistingExclusionRulesListBox .Layout.Column = 1;
        end
        
        function createCopyButton(this, parentGrid)
            iCreateDummyMinSizeButton(parentGrid, 1, 1);
            
            this.CopyButton = uibutton(parentGrid, ...
                'Text', getString(message('stats:dfittool:button_copy')), ...
                'ButtonPushedFcn', @(~,~) this.copyButtonClickedCallback(), ...
                'Tag', 'dfExcludePanelCopyButton');
            this.CopyButton.Layout.Row = 1;
            this.CopyButton.Layout.Column = 1;
        end
        
        function createViewButton(this, parentGrid)
            iCreateDummyMinSizeButton(parentGrid, 1, 2);

            this.ViewButton = uibutton(parentGrid, ...
                'Text', getString(message('stats:dfittool:button_view')), ...
                'ButtonPushedFcn', @(~,~) this.viewButtonClickedCallback(), ...
                'Tag', 'dfExcludePanelViewButton');
            this.ViewButton.Layout.Row = 1;
            this.ViewButton.Layout.Column = 2;
        end
        
        function createRenameButton(this, parentGrid)
            iCreateWiderDummyMinSizeButton(parentGrid, 2, 1);

            this.RenameButton = uibutton(parentGrid, ...
                'Text', getString(message('stats:dfittool:button_rename')), ...
                'ButtonPushedFcn', @(~,~) this.renameButtonClickedCallback(), ...
                'Tag', 'dfExcludePanelRenameButton');
            this.RenameButton.Layout.Row = 2;
            this.RenameButton.Layout.Column = 1;
        end
        
        function createDeleteButton(this, parentGrid)
            iCreateWiderDummyMinSizeButton(parentGrid, 2, 2);

            this.DeleteButton = uibutton(parentGrid, ...
                'Text', getString(message('stats:dfittool:button_delete')), ...
                'ButtonPushedFcn', @(~,~) this.deleteButtonClickedCallback(), ...
                'Tag', 'dfExcludePanelDeleteButton');
            this.DeleteButton.Layout.Row = 2;
            this.DeleteButton.Layout.Column = 2;
        end
        
        function createCloseButton(this, parentGrid, column)
            iCreateWiderDummyMinSizeButton(parentGrid, 1, column);

            this.CloseButton = uibutton(parentGrid, ...
                'Text', getString(message('stats:dfittool:button_close')), ...
                'ButtonPushedFcn', @(~,~) this.closeButtonClickedCallback(), ...
                'Tag', 'dfExcludePanelCloseButton');
            this.CloseButton.Layout.Row = 1;
            this.CloseButton.Layout.Column = column;
        end
    
        function createHelpButton(this, parentGrid, column)
            iCreateWiderDummyMinSizeButton(parentGrid, 1, column);

            this.HelpButton = uibutton(parentGrid, ...
                'Text', getString(message('stats:dfittool:button_help')), ...
                'ButtonPushedFcn', @(~,~) this.helpButtonClickedCallback(), ...
                'Tag', 'dfExcludePanelHelpButton');
            this.HelpButton.Layout.Row = 1;
            this.HelpButton.Layout.Column = column;
        end
        
        % Callbacks
        function lowerLimitTextFieldValueChangedCallback(this, ~, ~)
            this.notify('LowerLimitTextFieldValueChanged');
        end
        
        function upperLimitTextFieldValueChangedCallback(this, ~, ~)
            this.notify('UpperLimitTextFieldValueChanged');
        end
        
        function selectDataDropDownValueChanged(this, ~, ~)
            this.notify('SelectDataDropDownValueChanged');
        end
        
        function excludeGraphicallyButtonClickedCallback(this, ~, ~)
            this.notify('ExcludeGraphicallyButtonClicked');
        end
        
        function createExclusionRuleButtonClickedCallback(this, ~, ~)
            this.notify('CreateExclusionRuleButtonClicked');
        end
        
        function existingExclusionRulesListBoxValueChangedCallback(this, ~, ~)
            this.notify('ExistingExclusionRulesListBoxValueChanged');
        end
        
        function copyButtonClickedCallback(this, ~, ~)
            this.notify('CopyButtonClicked');
        end
        
        function viewButtonClickedCallback(this, ~, ~)
            this.notify('ViewButtonClicked');
        end
        
        function renameButtonClickedCallback(this, ~, ~)
            this.notify('RenameButtonClicked');
        end
        
        function deleteButtonClickedCallback(this, ~, ~)
            this.notify('DeleteButtonClicked');
        end
        
        function closeButtonClickedCallback(this, ~, ~)
            this.notify('CloseButtonClicked');
        end
        
        function helpButtonClickedCallback(this, ~, ~)
            this.notify('HelpButtonClicked');
        end
    end
end

% helpers
function position = iGetDefaultPanelPosition()
position = get(0,'defaultfigureposition');
position([3,4]) = [590 350];
end

function iCreateDummyMinSizeButton(grid, row, column)
dummyText = repmat('X',1,8);
minButton = uibutton(grid, 'Text', dummyText);
minButton.Layout.Row = row;
minButton.Layout.Column = column;
minButton.Visible = 'off';
end

function iCreateWiderDummyMinSizeButton(grid, row, column)
dummyText = repmat('X',1,10);
minButton = uibutton(grid, 'Text', dummyText);
minButton.Layout.Row = row;
minButton.Layout.Column = column;
minButton.Visible = 'off';
end