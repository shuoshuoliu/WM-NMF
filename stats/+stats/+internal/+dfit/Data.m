classdef Data < handle
  % Data    Distribution Fitter Data
  
  % Copyright 2019-2020 The MathWorks, Inc
  
    properties(Access = private)
        DistributionFitting
        
        Panel;
        
        DataWorkspaceDropDown;
        FrequencyWorkspaceDropDown;
        CensoringWorkspaceDropDown;
              
        % Keep track of original selected item to properly set drop down
        % values when users select a column or row.
        DataSelectedItem;
        CensoringSelectedItem;
        FrequencySelectedItem;
        
        % Keep track of original workspace values for input to "Select
        % Column or Row" functions.
        DataSelectedWorkspaceValue;
        CensoringSelectedWorkspaceValue;
        FrequencySelectedWorkspaceValue;
        
        DataSetNameEditField;
        
        DataPreviewPanel;
        PreviewAxes;
        DataSetsTable;
        NoDataSetsLabel;
        
        DataSelectColumnOrRowButton;
        CensoringSelectColumnOrRowButton;
        FrequencySelectColumnOrRowButton;
        CreateDataSetButton;

        SelectedDataSets = {};

        RenameButton;
        ViewButton;
        DeleteButton;
        SetBinRulesButton;
        
        IsAutoName = true;
        SelectColumnOrRowDialog;
        SelectColumnOrRowTable;
        
        RenamePanel;
        
        % Listeners   (cell of listener) Listeners to various events
        Listeners = {};        
    end
    
    events
        DataSetRenamed
        DataSetAdded
		DataSetsDeleted
    end
   
    methods (Static)
        function this = getInstance
            persistent Instance;
            if ~isa(Instance, 'stats.internal.dfit.Data') || ~isvalid(Instance)
                Instance = stats.internal.dfit.Data;
            end
            this = Instance;
        end
    end
   
    methods
        function showPanel(this)
            figure(this.Panel); % Brings the uifigure forward and makes it visible
        end  
        
        function delete(this)
            delete(this.Panel);
        end
       
        function openBinWidthRulesPanel(~, dsname)
            dataSet=find(stats.internal.dfit.getdsdb,'name',dsname);
            iCreateAndOpenBinRulesPanel(dataSet);
        end
        
        function updateManageDataSetsTable(this)

            dsdb = stats.internal.dfit.getdsdb;
            
            numDataSets = stats.internal.dfit.getNumDataBaseItems(dsdb);
   
            if numDataSets > 0
                d = cell(numDataSets, 3);
                r = 1;
                
                dset = down(dsdb);
                
                while(~isempty(dset))
                    if (dset.plot == 0)
                        d{r,1} = false;
                    else
                        d{r,1} = true;
                    end
                    if (dset.showbounds == 0)
                        d{r,2} = false;
                    else
                        d{r,2} = true;
                    end
                    d{r,3} = dset.name;
                    r = r+1;
                    dset = right(dset);
                end
                
                this.NoDataSetsLabel.Visible = 'off';
                this.DataSetsTable.Data = d;
                this.DataSetsTable.Visible = 'on';
            else
                this.DataSetsTable.Visible = 'off';
                this.NoDataSetsLabel.Visible = 'on';
            end
            this.enableButtons();
        end
    end
   
    methods(Access = private)
        function this = Data()
            this.DistributionFitting = stats.internal.dfit.DistributionFitting.getInstance;
            this.RenamePanel = stats.internal.dfit.Rename(stats.internal.dfit.getdsdb);
            
            this.layoutPanel();
            this.updateManageDataSetsTable();
            
            this.addListener(this.RenamePanel, 'RenameSuccessful', @this.handleRenameSuccessful);
            this.addListener(this.DistributionFitting, 'SessionClearing', @(~, ~)this.handleSessionClearing);
            this.addListener(stats.internal.dfit.BinRulesManager.getInstance(), 'SaveBinRuleAsDefault', @(~, ~)this.updateDataPreviewAndEnableCreateDataSetButton);
        end 
        
        function handleSessionClearing(this)
            if isvalid(this.DataSetsTable)
                this.delete();
            end
        end
                
        function handleRenameSuccessful(this, ~, renameEventData) 
            newName = renameEventData.NewName;
            oldName = renameEventData.OldName;
            this.updateManageDataSetsTable();
            notify(this, 'DataSetRenamed', stats.internal.dfit.RenameEventData(this, oldName, newName));
        end
        
        function handleRenameCancelled(~, renameObj)
            delete(renameObj);
        end
        
        % Value changed functions
        function cbkDataDropDownValueChanged(this)          
            if this.IsAutoName
                dsName = '';
            end
            this.DataSelectColumnOrRowButton.Enable = 'off';
            this.CreateDataSetButton.Enable = 'off';
            
            this.DataSelectedWorkspaceValue = this.DataWorkspaceDropDown.WorkspaceValue;
            this.DataSelectedItem = this.DataWorkspaceDropDown.Value;
            
            if isempty(this.DataSelectedWorkspaceValue) % drop down shows "select"
                this.showSelectData();
            else % data selected
                if (isvector(this.DataSelectedWorkspaceValue))
                    if this.IsAutoName
                        dsName = [this.DataSelectedItem ' data'];
                    end
                 else
                    if this.IsAutoName
                        dsName = [this.DataSelectedItem '(:,1) data'];
                    end
                    this.DataWorkspaceDropDown.Value = [this.DataSelectedItem '(:,1)'];
                    this.DataSelectColumnOrRowButton.Enable = 'on';
                end
            this.updateDataPreviewAndEnableCreateDataSetButton();
            end
            if this.IsAutoName
                this.DataSetNameEditField.Value = iUniqueName(dsName);
            end
        end
        
        function cbkCensoringDropDownValueChanged(this)   
            this.CensoringSelectedWorkspaceValue = this.CensoringWorkspaceDropDown.WorkspaceValue;
            this.CensoringSelectedItem = this.CensoringWorkspaceDropDown.Value;
            this.CensoringSelectColumnOrRowButton.Enable = 'off';
            if ~isempty(this.CensoringSelectedWorkspaceValue) && ...
                 ~isvector(this.CensoringSelectedWorkspaceValue)
                    this.CensoringWorkspaceDropDown.Value = [this.CensoringSelectedItem '(:,1)'];
                    this.CensoringSelectColumnOrRowButton.Enable = 'on';
            end
            this.updateDataPreviewAndEnableCreateDataSetButton();
        end
        
        function cbkFrequencyDropDownValueChanged(this)            
            this.FrequencySelectedWorkspaceValue = this.FrequencyWorkspaceDropDown.WorkspaceValue;
            this.FrequencySelectedItem = this.FrequencyWorkspaceDropDown.Value;
            this.FrequencySelectColumnOrRowButton.Enable = 'off';
            if ~isempty(this.FrequencySelectedWorkspaceValue) && ...
                ~isvector(this.FrequencySelectedWorkspaceValue)
                    this.FrequencyWorkspaceDropDown.Value = [this.FrequencySelectedItem '(:,1)'];
                    this.FrequencySelectColumnOrRowButton.Enable = 'on';
             end
            this.updateDataPreviewAndEnableCreateDataSetButton();
        end
        
        function cbkDataSelectColumnOrRowOKButton(this, f)
            this.setDropDownValue(this.DataWorkspaceDropDown, this.DataSelectedItem, f);
            if this.IsAutoName
                this.DataSetNameEditField.Value = iUniqueName([this.DataWorkspaceDropDown.Value ' data']);
            end
            this.updateDataPreviewAndEnableCreateDataSetButton;
        end
        
        function cbkSelectColumnOrRowOKButton(this, workspaceDropDown, selectedItem, f)
            this.setDropDownValue(workspaceDropDown, selectedItem, f);
            this.updateDataPreviewAndEnableCreateDataSetButton;
        end
        
        function cbkSelectColumnOrRowCancelButton(~, f)
            close(f);
        end        

        function cbkCreateDataSet(this) 
            
            if iIsNameInvalid(this, this.DataSetNameEditField.Value)
                return;
            end
            
            [ds, err, ~, ~] = stats.internal.dfit.createdataset( ...
                iGetValueWithNone(this.DataWorkspaceDropDown), ...
                iGetValueWithNone(this.CensoringWorkspaceDropDown), ...
                iGetValueWithNone(this.FrequencyWorkspaceDropDown), ...
                this.DataSetNameEditField.Value);
            if isempty(err)
                this.updateManageDataSetsTable();
                iSelectLastRow(this, ds);
                this.enableButtons();
                iSetDefaults(this);
                this.showSelectData();
                notify(this, 'DataSetAdded');
            end        
        end

        function cbkSelectColumnOrRowTableCellSelection(this, eventData)       
            indices = eventData.Indices;
            
            if all(size(indices) == [1 2])  ||  all(indices(1,2) == indices(:,2))  % all selected cells are in the same column 
                this.SelectColumnOrRowTable.Selection = iGetColumnCellsToSelect(eventData.Source.Data, indices(1,2));
            elseif all(indices(1,1) == indices(:,1))% all selected cells are in the same row
                this.SelectColumnOrRowTable.Selection = iGetRowCellsToSelect(eventData.Source.Data, indices(1,1));
            % else error - message will be displayed when user presses the
            % OK button. 
            end
        end
        
        function cbkManageDataSetsTableCellSelection(this, eventData)           
            this.SelectedDataSets = stats.internal.dfit.getSelectedObjects(...
               stats.internal.dfit.getdsdb, this.DataSetsTable, ...
               eventData.Indices, 3);
            
            this.enableButtons();
        end
        
        function cbkManageDataSetsTableCellEdited(this, eventData)   
            % Select the row being edited
            this.DataSetsTable.Selection = eventData.Indices(1);
            selectedDataSet = stats.internal.dfit.getSelectedObjects(...
               stats.internal.dfit.getdsdb, this.DataSetsTable, ...
               eventData.Indices, 3);
            if eventData.Indices(2) == 1 % plot selected
                selectedDataSet.plot = eventData.EditData;
            else % confbounds selected
                if eventData.EditData
                    stats.internal.dfit.boundwarn(selectedDataSet);
                end
                selectedDataSet.showbounds = eventData.EditData;
            end

            this.enableButtons();
        end

        function cbkViewButton(this)
            % View is enabled only if one dataSet is selected.
            dataSet = this.SelectedDataSets;
            if isempty(dataSet.viewData)
                viewData = stats.internal.dfit.ViewData(dataSet);
                dataSet.viewData = viewData;
            end
            dataSet.viewData.showPanel();
        end
        
        function cbkSetBinRulesButton(this)
            dataSet = this.SelectedDataSets;
            iCreateAndOpenBinRulesPanel(dataSet)
        end
        
        function cbkRenameButton(this)
            this.RenamePanel.showPanel(this.SelectedDataSets.name);
        end

        function cbkDeleteButton(this)                  
            selectedNames = {};
            for i = 1:length(this.SelectedDataSets)
                selectedNames{1, end + 1} = this.SelectedDataSets(i).name;  %#ok<AGROW>
            end     
            
            [warningMessage, associatedFits, associatedExclusionRules] = this.constructWarningMessage(selectedNames);
            
            if isempty(warningMessage)
                bOKtoDelete = true;
            else
                bOKtoDelete = false;
                selection = uiconfirm(this.Panel, warningMessage, ...
                    getString(message('stats:dfittool:title_deletingDataSets')), ...
                    'Icon','warning');
                if strcmp(selection, getString(message('MATLAB:uitools:uidialogs:OK')))
                    bOKtoDelete = true;
                end
            end
             
            if bOKtoDelete 

                % Delete fits first because they have a reference to data
                % sets.
                if ~isempty(associatedFits)
                    fm = stats.internal.dfit.FitsManager.getInstance;
                    fm.deleteFits(associatedFits);
                end

                if ~isempty(associatedExclusionRules)
                    associatedExclusionRuleNames = cell(1, length(associatedExclusionRules));
                    for i=1:length(associatedExclusionRules)
                        associatedExclusionRuleNames{i} = associatedExclusionRules(i).name;
                    end
                    exclude = stats.internal.dfit.Exclude.getInstance;
                    exclude.deleteExclusionRule(associatedExclusionRuleNames);
                end
                
                % It is important to send this notification after fits and
                % exclusion rules are deleted, since we don't want
                % notification going to fits and exclusions that might have
                % already been deleted since that were using the deleted
                % data set.
                
                for i = 1:length(this.SelectedDataSets)
                    if ~isempty(this.SelectedDataSets(i).viewData)
                        delete(this.SelectedDataSets(i).viewData);
                    end
                    % Actually delete the udd object
                    delete(this.SelectedDataSets(i));
                end
                this.DataSetsTable.Selection = [];
                this.SelectedDataSets = {};
                this.updateManageDataSetsTable();
                notify(this, 'DataSetsDeleted');
            end
        end

           function cbkHelpButton(~)
            stats.internal.dfit.helpviewer('import_data', 'data');
        end
        
        function cbkCloseButton(this)
            yesStr = getString(message('stats:dfittool:button_yes'));
            noStr  = getString(message('stats:dfittool:button_no'));
            cancelStr = getString(message('stats:dfittool:button_cancel'));
            
            if this.CreateDataSetButton.Enable
                selection  = uiconfirm(this.Panel, ...
                    getString(message('stats:dfittool:warning_notImported')), ...
                    getString(message('stats:dfittool:title_notImported')), ... 
                    'Options', {yesStr, noStr, cancelStr}, ...
                    'DefaultOption', 1, 'CancelOption', 3);
                if strcmp(selection, yesStr)
                    if iIsNameInvalid(this, this.DataSetNameEditField.Value)
                        return;
                    end
                    this.cbkCreateDataSet();
                    this.Panel.Visible = 'off';
                elseif strcmp(selection, noStr)
                    iSetDefaults(this);
                    this.showSelectData();
                    this.Panel.Visible = 'off';
                % else - cancel - no action 
                end
            else
                this.Panel.Visible = 'off';
            end
        end
        
        function cbkDataSetNameChanging(this)
            % once the user starts changing the the data sets name, don't
            % don't mess with it.
            this.IsAutoName = false;
        end
        
        function cbkDataSetNameChanged(this, event)
            this.dataSetName = event.Value;
            % Pressing enter in the Data set name field should create the
            % data set (if "Create Data Set" button is enabled). However we
            % can't yet tell the difference between a focus lost and
            % pressing enter (both generate a "ValueChanged". So we are not
            % currently responding to a "ValueChanged".
            
            if strcmp(this.CreateDataSetButton.Enable, 'off') 
                uialert(this.Panel,...
                    getString(message('stats:dfstrings:sprintf_DataMustBeSpecified')), ...
                    getString(message('stats:dfittool:title_invalidSelections')));
            else
                cbkCreateDataSet(this);
            end
        end
        
        function setDropDownValue(this, workspaceDropDown, selectedItem, f)
            % If the length of the indices is 2, user
            % selected a single cell - find the column number and use that.
            % Otherwise check to see if all the columns match or all the
            % rows match and then get the values from that row or column.

            value = selectedItem;
            indices = this.SelectColumnOrRowTable.Selection;
            if isempty(indices)
                stats.internal.dfit.errordlg(getString(message('stats:dfittool:error_no_selection')), getString(message('stats:dfittool:title_dataSelector')));
                return;
            end
                
            goodSelection = true;
            if length(indices) == 2  % a single cell was selected; take the column
                value = [value '(:,' num2str(indices(2)), ')'];
            elseif all(indices(:,2) == indices(1,2)) % column selected
                value = [value '(:,' num2str(indices(1,2)), ')'];
            elseif all(indices(:,1) == indices(1,1)) % row selected
                value = [value '(' num2str(indices(1,1)), ',:)'];
            else
                stats.internal.dfit.errordlg(getString(message('stats:dfittool:error_bad_selection')), getString(message('stats:dfittool:title_dataSelector')));
                goodSelection = false;
            end
            
            if goodSelection
                workspaceDropDown.Value = value;
                close(f);
            end
        end
        
        function showSelectData(this)
            stats.internal.dfit.DFViewUtilities.createMessagePreview(this.DataPreviewPanel, ...
                getString(message('stats:dfittool:label_importMsg1')));
        end
        
        function updateDataPreviewAndEnableCreateDataSetButton(this)
            % updateDataPreviewAndEnableCreateDataSetButton updates the 
            % preview panel and sets the "createDataSet" button enable state.
            
            % If all three dropdowns are empty, just display "Select data"
            if (isempty(this.DataWorkspaceDropDown.WorkspaceValue) && ...
                    isempty(this.CensoringWorkspaceDropDown.WorkspaceValue) && ...
                    isempty(this.FrequencyWorkspaceDropDown.WorkspaceValue))
                isValid = false;
                this.showSelectData();
            else
                isValid = stats.internal.dfit.DFViewUtilities.createDataSetDistributionPreviewFromExpressions(...
                    this.DataPreviewPanel, ...
                    iGetValueWithEmpty(this.DataWorkspaceDropDown), ...
                    iGetValueWithEmpty(this.CensoringWorkspaceDropDown), ...
                    iGetValueWithEmpty(this.FrequencyWorkspaceDropDown));
            end
            if isValid
                this.CreateDataSetButton.Enable = 'on';
            else
                this.CreateDataSetButton.Enable = 'off';
            end
        end
        
        function associatedFits = findMoreFitsAssociatedWithExclusionRules(~, aF, outliers)
            % Check to see if outliers are associated with any fits not
            % already in the list. Add new ones to the lists.
            associatedFits = aF;
            moreAssociatedFits = [];
            for i = 1:length(outliers)
                moreAssociatedFits = cat(1, moreAssociatedFits, find(stats.internal.dfit.getfitdb, 'exclusionrule', outliers(i)));
            end
            for i = 1:length(moreAssociatedFits)
                if (~ismember(moreAssociatedFits(i), associatedFits))
                    associatedFits(end + 1) = moreAssociatedFits(i); %#ok<AGROW>
                end
            end
        end
        
        function [warningMessage, associatedFits, associatedExclusionRules] = constructWarningMessage(this, selectedNames)
            warningMessage = {};
            line = 1;
            % For each dataset
            for i = 1:length(selectedNames)
                ds = find(stats.internal.dfit.getdsdb, 'name', selectedNames{i});
                associatedFits = find(stats.internal.dfit.getfitdb, 'dataset', ds.name);
                associatedExclusionRules = find(stats.internal.dfit.getoutlierdb, 'dataset', ds.name);

                associatedFits = this.findMoreFitsAssociatedWithExclusionRules(associatedFits, associatedExclusionRules);
                if ~isempty(associatedFits) || ~isempty(associatedExclusionRules)
                    warningMessage{line} = getString(message('stats:dfittool:warning_deletingFit', ds.name)); %#ok<AGROW>
                    line = line + 1;
                    % List fits first
                    if ~isempty(associatedFits)
                        warningMessage{line} = this.constructWarningLine( ...
                            getString(message('stats:dfittool:warning_fit')), ...
                            getString(message('stats:dfittool:warning_fits')), ...
                            associatedFits); %#ok<AGROW>
                        line = line + 1;
                    end
                    % Then list exclusion rules
                    if ~isempty(associatedExclusionRules)
                        warningMessage{line} = this.constructWarningLine( ...
                            getString(message('stats:dfittool:label_exclusionRule')), ...
                            getString(message('stats:dfittool:label_exclusionRules')), ...
                            associatedExclusionRules); %#ok<AGROW>
                        line = line + 1;
                    end
                end
                % Add an empty line between data sets (if there are any)
                if ~isempty(warningMessage)
                    warningMessage{line} = ''; %#ok<AGROW>
                    line = line +1;
                end
            end
        end
        
        function warningLine = constructWarningLine(~, oneItemString, moreThanOneItemString, dbItems)
            if ~isempty(dbItems)
                if (length(dbItems) == 1)
                    warningLine = ['    ', oneItemString, ' ', dbItems(1).name];
                else
                    warningLine = ['    ', moreThanOneItemString, ' ', dbItems(1).name];
                    for i = 2:length(dbItems)
                        warningLine = [warningLine, ', ', dbItems(i).name]; %#ok<AGROW>
                    end
                end
            end
        end
        
        function enableButtons(this)
            % New fit should always be enabled (so is not included here)
            % Delete fit should be enabled if at least one fit is selected.
            % Copy, Edit and Save to workspace should be enabled if one and
            % only one fit is selected.
            switch numel(this.SelectedDataSets)
                case 0
                    this.RenameButton.Enable = 'off';
                    this.ViewButton.Enable = 'off';
                    this.DeleteButton.Enable = 'off';
                    this.SetBinRulesButton.Enable = 'off';
                case 1
                    this.RenameButton.Enable = 'on';
                    this.ViewButton.Enable = 'on';
                    this.DeleteButton.Enable = 'on';
                    this.SetBinRulesButton.Enable = 'on';
                otherwise
                    this.RenameButton.Enable = 'off';
                    this.ViewButton.Enable = 'off';
                    this.DeleteButton.Enable = 'on';
                    this.SetBinRulesButton.Enable = 'off';
            end
        end
        
        function addListener(this, obj, eventName, callback)
            this.Listeners{end+1} = event.listener(obj, eventName, callback);
        end
        
        function layoutPanel(this)
            this.Panel = uifigure('Tag', 'dfDataUIFigure', ...
                'Name', getString(message('stats:dfittool:title_data')), ...
                'Position', [680 350 580 675], ...
                'CloseRequestFcn', @(~, ~)this.cbkCloseButton(), ...
                'Visible', 'off');
                     
            mainG = uigridlayout(this.Panel, [2, 1]);

            mainG.RowHeight = {196, '1x'};
            topG = uigridlayout(mainG, [6, 4]);
            topG.ColumnWidth = {'fit', 'fit', 'fit', 180};
            topG.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', '1x'};
            topG.Layout.Row = 1;
            topG.Layout.Column = 1;
            l1 = uilabel(topG, 'Text', getString(message('stats:dfittool:label_importVars')));
            l1.Layout.Row = 1;
            l1.Layout.Column = [1 3];
            l2 = uilabel(topG, 'Text', getString(message('stats:dfittool:label_preview')));
            l2.Layout.Row = 1;
            l2.Layout.Column = 4;
            
            l3 = uilabel(topG, 'Text', ['  ' getString(message('stats:dfittool:label_dataset'))]);
            l3.Layout.Row = 2;
            l3.Layout.Column = 1;
            l4 = uilabel(topG, 'Text', ['  ' getString(message('stats:dfittool:label_censoring'))]);
            l4.Layout.Row = 3;
            l4.Layout.Column = 1;
            l5 = uilabel(topG, 'Text', ['  ' getString(message('stats:dfittool:label_frequency'))]);
            l5.Layout.Row = 4;
            l5.Layout.Column = 1;
            l6 = uilabel(topG, 'Text', getString(message('stats:dfittool:label_dataSetName')));
            l6.Layout.Row = 5;
            l6.Layout.Column = 1;
                        
            this.DataWorkspaceDropDown = ...
                matlab.ui.control.internal.model.WorkspaceDropDown(...
                'Parent', topG, 'Editable', 'on', ...
                'Tag', 'dfDataDataWorkspaceDropDown');
            this.DataWorkspaceDropDown.FilterVariablesFcn = @(x) isnumeric(x) & ismatrix(x);
            this.DataWorkspaceDropDown.ValueChangedFcn = @(~,~)this.cbkDataDropDownValueChanged();
            this.DataWorkspaceDropDown.Layout.Row = 2;
            this.DataWorkspaceDropDown.Layout.Column = 2;
            
            this.CensoringWorkspaceDropDown = ...
                matlab.ui.control.internal.model.WorkspaceDropDown(...
                'Parent', topG, 'Editable', 'on', ...
                'Tag', 'dfDataCensoringWorkspaceDropDown');
            this.CensoringWorkspaceDropDown.FilterVariablesFcn = @(x) (isnumeric(x) || islogical(x)) & ismatrix(x);
            this.CensoringWorkspaceDropDown.ValueChangedFcn = @(~, ~)this.cbkCensoringDropDownValueChanged();
            this.CensoringWorkspaceDropDown.Layout.Row = 3;
            this.CensoringWorkspaceDropDown.Layout.Column = 2;
            
            this.FrequencyWorkspaceDropDown = ...
                matlab.ui.control.internal.model.WorkspaceDropDown(...
                'Parent', topG, 'Editable', 'on', ...
                'Tag', 'dfDataFrequencyWorkspaceDropDown');
            this.FrequencyWorkspaceDropDown.FilterVariablesFcn = @(x) isnumeric(x) & ismatrix(x);
            this.FrequencyWorkspaceDropDown.ValueChangedFcn = @(~,~)this.cbkFrequencyDropDownValueChanged();
            this.FrequencyWorkspaceDropDown.Layout.Row = 4;
            this.FrequencyWorkspaceDropDown.Layout.Column = 2;
 
            this.DataSetNameEditField = uieditfield(topG, 'ValueChangingFcn', @(~, ~)this.cbkDataSetNameChanging(),'Tag','dfDataDataSetNameEditField');
            this.DataSetNameEditField.Layout.Row = 5;
            this.DataSetNameEditField.Layout.Column = 2;
            
            iCreateDummyMinSizeButton(topG, 2, 3)
            this.DataSelectColumnOrRowButton = uibutton(topG, 'Text', ...
                getString(message('stats:dfittool:button_selectColumnOrRow')), ...
                'Enable', 'off', ...
                'Tag', 'dfDataDataSelectColumnOrRowButton', ...
                'ButtonPushedFcn', @(~,~)this.cbkDataSelectColumnOrRow());
            this.DataSelectColumnOrRowButton.Layout.Row = 2;
            this.DataSelectColumnOrRowButton.Layout.Column = 3;
            iCreateDummyMinSizeButton(topG, 3, 3)
            this.CensoringSelectColumnOrRowButton = uibutton(topG, 'Text', ...
                getString(message('stats:dfittool:button_selectColumnOrRow')), ...
                'Enable', 'off', ...
                'Tag', 'dfDataCensoringSelectColumnOrRowButton', ...
                'ButtonPushedFcn', @(src,event)this.cbkSelectColumnOrRow(src, event, this.CensoringWorkspaceDropDown, this.CensoringSelectedWorkspaceValue, this.CensoringSelectedItem));
            this.CensoringSelectColumnOrRowButton.Layout.Row = 3;
            this.CensoringSelectColumnOrRowButton.Layout.Column = 3;
            iCreateDummyMinSizeButton(topG, 4, 3)
            this.FrequencySelectColumnOrRowButton = uibutton(topG, 'Text', ...
                getString(message('stats:dfittool:button_selectColumnOrRow')), ...
                'Enable', 'off', ...
                'Tag', 'dfDataFrequencySelectColumnOrRowButton', ...
                'ButtonPushedFcn', @(src,event)this.cbkSelectColumnOrRow(src, event, this.FrequencyWorkspaceDropDown, this.FrequencySelectedWorkspaceValue, this.FrequencySelectedItem));
            this.FrequencySelectColumnOrRowButton.Layout.Row = 4;
            this.FrequencySelectColumnOrRowButton.Layout.Column = 3;
            iCreateDummyMinSizeButton(topG, 6, 3)
            this.CreateDataSetButton = uibutton(topG, 'Text', ...
                getString(message('stats:dfittool:button_createDataSet')), ...
                'ButtonPushedFcn', @(~, ~)this.cbkCreateDataSet(), ...
                'Tag', 'dfDataCreateDataSetButton', ...
                'Enable', 'off');
            this.CreateDataSetButton.Layout.Row = 6;
            this.CreateDataSetButton.Layout.Column = 3;
            
            this.DataPreviewPanel = uipanel(topG, 'BackgroundColor', [1,1,1]);
            this.DataPreviewPanel.Layout.Row = [2 6];
            this.DataPreviewPanel.Layout.Column = 4;
            
            stats.internal.dfit.DFViewUtilities.createMessagePreview(this.DataPreviewPanel, getString(message('stats:dfittool:label_importMsg1')));
            
            bottomG = uigridlayout(mainG, [8 1]);
            bottomG.RowHeight = {'fit', 3, '1x', 80, 3, 'fit', 3, 'fit'};
            bottomG.Padding = [0 0 0 0];
            bottomG.ColumnWidth = {'1x'};
            bottomG.Layout.Row = 2;
            bottomG.Layout.Column = 1;
            
            l3 = uilabel(bottomG, 'Text', getString(message('stats:dfittool:label_manageDataSets')));
            l3.Layout.Row = 1;
            l3.Layout.Column = 1;
            
            p1 = uipanel(bottomG);
            p1.BackgroundColor = [.75 .75 .75];
            p1.Layout.Row = 2;
            p1.Layout.Column = 1;
                       
            this.NoDataSetsLabel = uilabel(bottomG, 'Text', getString(message('stats:dfittool:label_noDataSets')), 'Tag', 'dfDataPanelNoDataSetsYet');
            this.NoDataSetsLabel.Layout.Row = 3;
            this.NoDataSetsLabel.Layout.Column = 1;
            
            plotExtent = iGetTextWidth([' ' getString(message('stats:dfittool:header_plot'))]);
            boundsExtent = iGetTextWidth([' ' getString(message('stats:dfittool:header_bounds'))]);
            this.DataSetsTable = uitable(bottomG, ...
                'Tag', 'dfDataDataSetsTable', 'Visible', 'off', ...
                'CellSelectionCallback', ...
                @(~, event)this.cbkManageDataSetsTableCellSelection(event), ...
                'CellEditCallback', ...
                @(~, event)this.cbkManageDataSetsTableCellEdited(event), ...
                'ColumnWidth', {plotExtent boundsExtent 'auto'}, ...
                'ColumnName', ...
                {getString(message('stats:dfittool:header_plot')), ...
                getString(message('stats:dfittool:header_bounds')), ...
                getString(message('stats:dfittool:header_setName'))}, ...
                'SelectionType', 'row', ...
                'RowStriping', 'off', ...
                'RowName', {}, ...
                'RowStriping', 'off', ...
                'ColumnEditable', [true, true, false]);
                        
            this.DataSetsTable.Layout.Row = [3 4];
            this.DataSetsTable.Layout.Column = 1;
            
            % Empty panels to simulate horizontal lines
            p1 = uipanel(bottomG);
            p1.BackgroundColor = [.75 .75 .75];
            p1.Layout.Row = 5;
            p1.Layout.Column = 1;
            
            buttonGrid1 = uigridlayout(bottomG, [1 6]);
            buttonGrid1.Padding = [0 0 0 0];
            buttonGrid1.ColumnWidth = {'1x',  'fit',  'fit',  'fit',  'fit',  '1x'};
            buttonGrid1.Layout.Row = 6;
            buttonGrid1.Layout.Column = 1;
            
            % Empty panels to simulate horizontal lines
            p1 = uipanel(bottomG);
            p1.BackgroundColor = [.75 .75 .75];
            p1.Layout.Row = 7;
            p1.Layout.Column = 1;
            
            iCreateDummyMinSizeButton(buttonGrid1, 1, 2)
            this.ViewButton = uibutton(buttonGrid1, 'Text', ...
                getString(message('stats:dfittool:button_view')), ...
                'Enable', 'off', ...
                'Tag', 'dfDataViewButton', ...
                'ButtonPushedFcn', @(~, ~)this.cbkViewButton());
            this.ViewButton.Layout.Row = 1;
            this.ViewButton.Layout.Column = 2;
            iCreateDummyMinSizeButton(buttonGrid1, 1, 3)
            this.SetBinRulesButton = uibutton(buttonGrid1, ...
                'Enable', 'off', ...
                'Tag', 'dfDataSetBinRulesButton', ...
                'Text', getString(message('stats:dfittool:button_binWidth')), ...
                'ButtonPushedFcn', @(~, ~)this.cbkSetBinRulesButton());
            this.SetBinRulesButton.Layout.Row = 1;
            this.SetBinRulesButton.Layout.Column = 3;
            
            iCreateDummyMinSizeButton(buttonGrid1, 1, 4)
            this.RenameButton = uibutton(buttonGrid1, ...
                'Text', getString(message('stats:dfittool:button_rename')), ...
                'Enable', 'off', ...
                'Tag', 'dfDataRenameButton', ...
                'ButtonPushedFcn', @(~, ~)this.cbkRenameButton());
            this.RenameButton.Layout.Row = 1;
            this.RenameButton.Layout.Column = 4;
            
            iCreateDummyMinSizeButton(buttonGrid1, 1, 5)
            this.DeleteButton = uibutton(buttonGrid1, ...
                'Text', getString(message('stats:dfittool:button_delete')), ...
                'Enable', 'off', ...
                'Tag', 'dfDataDeleteButton', ...
                'ButtonPushedFcn', @(~, ~)this.cbkDeleteButton());
            this.DeleteButton.Layout.Row = 1;
            this.DeleteButton.Layout.Column = 5;
            
            buttonGrid2 = uigridlayout(bottomG, [1 3]);
            buttonGrid2.Padding = [0 0 0 0];
            buttonGrid2.ColumnWidth = {'fit',  '1x',  'fit'};
            buttonGrid2.Layout.Row = 8;
            buttonGrid2.Layout.Column = 1;
            
            iCreateDummyMinSizeButton(buttonGrid2, 1, 1)
            helpButton = uibutton(buttonGrid2, ...
                'Text', getString(message('stats:dfittool:button_help')), ...
                'Tag', 'dfDataHelpButton', ...
                'ButtonPushedFcn', @(~, ~)this.cbkHelpButton());
            helpButton.Layout.Row = 1;
            helpButton.Layout.Column = 1;
            
            iCreateDummyMinSizeButton(buttonGrid2, 1, 3)
            closeButton = uibutton(buttonGrid2, ...
                'Text', getString(message('stats:dfittool:button_close')), ...
                'Tag', 'dfDataCloseButton', ...
                'ButtonPushedFcn', @(~, ~)this.cbkCloseButton());
            closeButton.Layout.Row = 1;
            closeButton.Layout.Column = 3;
            
            stats.internal.dfit.centergui(this.Panel);
        end
        
        function cbkDataSelectColumnOrRow(this)
            setButtonEnablement(this.DataSelectColumnOrRowButton,false)
            cleanup = onCleanup(@()setButtonEnablement(this.DataSelectColumnOrRowButton,true));
            if isempty(this.SelectColumnOrRowDialog) || ~isgraphics(this.SelectColumnOrRowDialog)
            [this.SelectColumnOrRowDialog, okButton] = this.layoutSelectColumnOrRow(this.DataWorkspaceDropDown, this.DataSelectedWorkspaceValue, this.DataSelectedItem);
            okButton.ButtonPushedFcn = @(~, ~)this.cbkDataSelectColumnOrRowOKButton(this.SelectColumnOrRowDialog);
            end
        end
        
        function cbkSelectColumnOrRow(this, source, ~, workspaceDropDown, selectedWorkspaceValue, selectedItem)
            setButtonEnablement(source,false)
            cleanup = onCleanup(@()setButtonEnablement(source,true));
            if isempty(this.SelectColumnOrRowDialog) || ~isgraphics(this.SelectColumnOrRowDialog)
            this.SelectColumnOrRowDialog = this.layoutSelectColumnOrRow(workspaceDropDown, selectedWorkspaceValue, selectedItem);
            end
        end
        
        function [f, okButton] = layoutSelectColumnOrRow(this, workspaceDropDown, selectedWorkspaceValue, selectedItem)
            f = uifigure('Name', getString(message('stats:dfittool:title_dataSelector')), 'WindowStyle', 'modal');
            
            mainG = uigridlayout(f, [2, 1]);
            mainG.RowHeight = {'1x', 'fit'};
            mainG.ColumnWidth = {'1x'};
            
            this.SelectColumnOrRowTable = uitable(mainG, 'Data', selectedWorkspaceValue, 'CellSelectionCallback', ...
                @(~, event)this.cbkSelectColumnOrRowTableCellSelection(event));
            this.SelectColumnOrRowTable.Layout.Row = 1;
            this.SelectColumnOrRowTable.Layout.Column = 1;
            
            buttonGrid = uigridlayout(mainG, [1, 3]);
            buttonGrid.RowHeight = {'fit'};
            buttonGrid.ColumnWidth = {'1x', 'fit', 'fit'};
            buttonGrid.Layout.Row = 2;
            buttonGrid.Layout.Column = 1;
            
            iCreateDummyMinSizeButton(buttonGrid, 1, 2)
            okButton = uibutton(buttonGrid, ...
                'Text', getString(message('stats:dfittool:button_OK')), ...
                'Tag', 'dfDataSelectColumnOrRowOKButton', ...
                'ButtonPushedFcn', @(~, ~)this.cbkSelectColumnOrRowOKButton(workspaceDropDown, selectedItem, f));
            okButton.Layout.Row = 1;
            okButton.Layout.Column = 2;
            iCreateDummyMinSizeButton(buttonGrid, 1, 3)
            cancelButton = uibutton(buttonGrid, ...
                'Text', getString(message('stats:dfittool:button_cancel')), ...
                'Tag', 'dfDataSelectColumnOrRowCancelButton', ...
                'ButtonPushedFcn', @(~, ~)this.cbkSelectColumnOrRowCancelButton(f));
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 3;
        end
    end
end

function dsName = iUniqueName(dsNameBase)
    copy = 2;
    notUnique = true;
    dsName = dsNameBase;
    while (notUnique)
        ds = find(stats.internal.dfit.getdsdb,'name',dsName);
        if isempty(ds)
            notUnique = false;
        else
            dsName = [dsNameBase ' (' num2str(copy) ')'];
            copy = copy + 1;
        end
    end
end
        
function value = iGetValueWithEmpty(workspaceDropDown)
    if isempty(workspaceDropDown.WorkspaceValue)
        value = [];
    else
        value = workspaceDropDown.Value;
    end
end

function value = iGetValueWithNone(workspaceDropDown)
    if isempty(workspaceDropDown.WorkspaceValue)
        % This string is intentionally not translated
        value = '(none)';
    else
        value = workspaceDropDown.Value;
    end
end

function cellsToSelect = iGetColumnCellsToSelect(d, selectedColumn)
    s = size(d);
    numRows = s(1);
    cellsToSelect = [(1:numRows)', selectedColumn*(ones(numRows,1))];
end

function cellsToSelect = iGetRowCellsToSelect(d, selectedRow)
    s = size(d);
    numColumns = s(2);
    cellsToSelect = [(selectedRow*(ones(numColumns,1))), (1:numColumns)'];
end


function iSetDefaults(this)
    this.DataWorkspaceDropDown.Value = getString(message('MATLAB:ui:defaults:select'));
    this.CensoringWorkspaceDropDown.Value = getString(message('MATLAB:ui:defaults:select'));
    this.FrequencyWorkspaceDropDown.Value = getString(message('MATLAB:ui:defaults:select'));
    this.DataSelectedItem = '';
    this.CensoringSelectedItem = '';
    this.FrequencySelectedItem = '';  
    this.DataSelectedWorkspaceValue = [];
    this.CensoringSelectedWorkspaceValue = [];
    this.FrequencySelectedWorkspaceValue = [];
    this.DataSelectColumnOrRowButton.Enable = 'off';
    this.CensoringSelectColumnOrRowButton.Enable = 'off';
    this.FrequencySelectColumnOrRowButton.Enable = 'off';
    this.CreateDataSetButton.Enable = 'off';
    this.DataSetNameEditField.Value = '';
    this.IsAutoName = true;
end

function iCreateDummyMinSizeButton(grid, row, column)
    dummyText = repmat('X',1,8);
    minButton = uibutton(grid, 'Text', dummyText);
    minButton.Layout.Row = row;
    minButton.Layout.Column = column;
    minButton.Visible = 'off';
end

function iCreateAndOpenBinRulesPanel(dataSet)
    if isempty(dataSet.binRulesPanel)
        binRulesPanel = stats.internal.dfit.DataSetBinRulesPanel(dataSet);
        dataSet.binRulesPanel = binRulesPanel;
    end
    dataSet.binRulesPanel.showPanel();
end

function width = iGetTextWidth(str)
    f = figure('Visible', 'off');
    txtCtrl = uicontrol('Parent', f, 'Style', 'Text', 'FontSize', 12, 'String', str);
    ext = txtCtrl.Extent;
    width = max(40, ext(3));
    delete(f);
end

function iSelectLastRow(this, ds)
    this.SelectedDataSets = ds;
    sizeData = size(this.DataSetsTable.Data);
    this.DataSetsTable.Selection = sizeData(1);
end

function tf = iIsNameInvalid(this, name)
% iIsNameInvalid alerts if names are empty, just spaces or already exists.
    tf = false;
    
    dataSetTitle = getString(message('stats:dfittool:label_dataSetUpperCase'));
    dataSetMessage = getString(message('stats:dfittool:label_dataSetLowerCase'));
    dataSetMessage2 = getString(message('stats:dfittool:label_dataSetMixedCase'));

    if strcmp(name, "")        
        uialert(this.Panel, ...
            getString(message('stats:dfittool:error_noName', dataSetMessage)), ...
            getString(message('stats:dfittool:title_noObjectTypeName', dataSetTitle))); 
        tf = true;
        return;
    end
    
    % just spaces?
    trimmedName = strtrim(name);
    if strcmp(trimmedName, "")               
        uialert(this.Panel, ...
            getString(message('stats:dfittool:error_noOnlySpacesAllowed')), ...
            getString(message('stats:dfittool:title_invalidObjectTypeName', dataSetTitle)));  
        tf = true;
        return;
    end
    
    this.DataSetNameEditField.Value = trimmedName;
              
    names = stats.internal.dfit.getDataBaseNames(stats.internal.dfit.getdsdb);

    % Does name already exist?
    if ismember(trimmedName, names)
        uialert(this.Panel, ...
            getString(message('stats:dfittool:error_duplicateName', dataSetMessage2, trimmedName)), ...
            getString(message('stats:dfittool:title_duplicateName', dataSetTitle)));
        tf = true;
        return;
    end
end

function setButtonEnablement(buttonObj,tf)
buttonObj.Enable = tf;
drawnow();
end
