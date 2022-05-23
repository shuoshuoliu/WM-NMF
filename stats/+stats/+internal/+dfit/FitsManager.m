classdef FitsManager < handle
   %FITSMANAGER  Dialog to manage fits
   
   %   Copyright 2019-2020 The MathWorks, Inc. 
   
   properties (GetAccess = 'private', SetAccess = 'private')
             
      DistributionFitting;      
      Panel;
      
      FitTable;
      SelectedFits = {};
      
      DisabledStyle;
      
      NewFitButton;
      CopyButton
      EditButton;
      SaveToWorkspaceButton;
      DeleteButton;
      
      SessionClearingListener
      DataSetRenamedListener;
      ExclusionRuleRenamedListener;
   end
   
   events
       FitAdded
       FitChanged
       FitsDeleted
   end
   
   methods
       
       function this = showPanel(this)
           figure(this.Panel); % Brings the uifigure forward and sets it visible
       end
       
       function delete(this)
           delete(this.Panel);
       end
       
       function fitAdded(this)
           this.updateFitsTable();
           notify(this, 'FitAdded');
       end
       
       function fitChanged(this, fitObj, fitFrame, oldName, newName)
           eventData = stats.internal.dfit.FittingEventData(fitObj, fitFrame, oldName, newName);
           notify(this, 'FitChanged', eventData);
           this.updateFitsTable();
       end
              
       function deleteFits(this, fits) % fits are the udd objects.
           for i = 1:length(fits)
               iUpdateSelectedFitsList(this, fits(i));
               delete(fits(i).fitframe);
               delete(fits(i));
           end
           this.updateFitsTable();
           notify(this, 'FitsDeleted');
       end
       
       function updateFitsTable(this)
           fitdb = stats.internal.dfit.getfitdb;
           numFits = stats.internal.dfit.getNumDataBaseItems(fitdb);
           
           removeStyle(this.FitTable);
           if isempty(this.SelectedFits)
                selection = [];
           else
                selection = zeros(1, length(this.SelectedFits));
                selectionIndex = 1;
           end
           
           if numFits > 0
               fits = cell(numFits, 5);
               r = 1;
               f = down(fitdb);

               while(~isempty(f))
                   if (f.plot == 0)
                       fits{r, 1} = false;
                   else
                       fits{r, 1} = true;
                   end
                   if (f.showbounds == 0)
                       fits{r, 2} = false;
                   else
                       fits{r, 2} = true;
                   end
                   fits{r, 3} = f.name;
                   fits{r, 4} = f.dataset;
                   fits{r, 5} = iGetDistributionDisplayName(this, f);
                   % if "bad" fit, uncheck plot and conf bounds (but do not
                   % change those values in the database)
                   if ~(f.isgood)
                       fits{r, 1} = false;
                       fits{r, 2} = false;
                       addStyle(this.FitTable, this.DisabledStyle, 'cell', [r, 1; r, 2]);
                   end
                   if iIsSelectedFits(this, f.name)
                       selection(selectionIndex) = r; %#ok<AGROW>
                       selectionIndex = selectionIndex + 1;
                   end
                   r = r+1;
                   f = right(f);
               end
           else
               fits = {};
           end
           
           if ishandle(this.FitTable)
               this.FitTable.Data = fits;
               this.FitTable.Selection = selection;
           end
           this.enableButtons();
       end
   end
   
  methods (Static)
       function this = getInstance
           persistent Instance;
           if ~isa(Instance, 'stats.internal.dfit.FitsManager') || ~isvalid(Instance)
               Instance = stats.internal.dfit.FitsManager;
               Instance.Panel = layoutPanel(Instance);
               Instance.updateFitsTable();
           end
           this = Instance;
       end
   end
   
   methods(Access = private)
       function this = FitsManager()
           this.DistributionFitting = stats.internal.dfit.DistributionFitting.getInstance;
           this.SessionClearingListener = event.listener(this.DistributionFitting, 'SessionClearing', @(~, ~)this.sessionClearing);
           this.DataSetRenamedListener = event.listener(stats.internal.dfit.Data.getInstance, 'DataSetRenamed', @(~, evt) this.handleDataSetRenamed(evt));
           this.ExclusionRuleRenamedListener = event.listener(stats.internal.dfit.Exclude.getInstance, 'ExclusionRuleRenamed', @(~, evt) this.handleExclusionRuleRenamed(evt));         
           this.DisabledStyle = uistyle('BackgroundColor', [0.9400 0.9400 0.9400]);
       end
       
       function sessionClearing(this)
           if isvalid(this)
               delete(this);
           end
       end
       
       function handleDataSetRenamed(this, evt)
           fit = down(stats.internal.dfit.getfitdb);
           while(~isempty(fit))
               if strcmp(fit.dataset, evt.OldName)
                   fit.dataset = evt.NewName;
               end
               fit = right(fit);
           end
           this.updateFitsTable();
       end
       
       function handleExclusionRuleRenamed(~, evt)
           fit = down(stats.internal.dfit.getfitdb);
           while(~isempty(fit))
               if strcmp(fit.exclusionrulename, evt.OldName)
                   fit.exclusionrulename = evt.NewName;
               end
               fit = right(fit);
           end
       end

       function cbEditButton(this)
           if isempty(this.SelectedFits.fitframe) % probably from load session
               stats.internal.dfit.FitEditor(this.SelectedFits);
           else
               showPanel(this.SelectedFits.fitframe);
           end
       end
       
       function cbCopyButton(this)
           % There can only be 1 selected fit for the Copy Button to be
           % enabled.
           newUddFit = copyfit(this.SelectedFits);
           stats.internal.dfit.FitEditor(newUddFit);
       end
       
       function cbNewFitButton(~)
           stats.internal.dfit.FitEditor();
       end
       
       function cbSaveToWorkspaceButton(this)
            stats.internal.dfit.save2ws(this.SelectedFits.name);
       end

       function panel = layoutPanel(this)
           panel = uifigure('Tag', 'dfFitsManagerUIFigure', ...
               'Name', getString(message('stats:dfittool:title_fitManager')), ...
               'Position', [680 650 560 380], ...
               'CloseRequestFcn', @(~, ~)this.cbCloseButton(), ...
               'Visible', 'off');
           
           stats.internal.dfit.centergui(panel);
           
           mainGrid = uigridlayout(panel, [5, 1]);
           mainGrid.RowHeight = {'fit', '1x', 'fit', 2, 'fit'};
           mainGrid.ColumnWidth = {'1x'};
           
           l = uilabel(mainGrid, 'Text', getString(message('stats:dfittool:label_tableOfFits')));
           l.Layout.Row = 1;
           l.Layout.Column = 1;
           
           buttonGrid1 = uigridlayout(mainGrid, [1, 7]);
           buttonGrid1.Padding = [0 0 0 0];
           buttonGrid1.RowHeight = {'fit'};
           buttonGrid1.ColumnWidth = {'1x', 'fit', 'fit', 'fit', 'fit', 'fit', '1x'};
           buttonGrid1.Layout.Row = 3;
           buttonGrid1.Layout.Column = 1;
           
           % Empty panels to simulate horizontal lines
           p = uipanel(mainGrid);
           p.BackgroundColor = [.75 .75 .75];
           p.Layout.Row = 4;
           p.Layout.Column = 1;
           
           buttonGrid2 = uigridlayout(mainGrid, [1, 3]);
           buttonGrid2.Padding = [0 0 0 0];
           buttonGrid2.RowHeight = {'fit'};
           buttonGrid2.ColumnWidth = {'fit', '1x', 'fit'};
           buttonGrid2.Layout.Row = 5;
           buttonGrid2.Layout.Column = 1;
           
           % Add the table
           this.FitTable = uitable(mainGrid, ...
               'Tag', 'dfFitsManagerTable', ...
               'CellSelectionCallback', ...
               @(~, event)this.cbManageFitsTableCellSelection(event), ...
               'CellEditCallback', ...
               @(~, event)this.cbManageFitsTableCellEdited(event), ...
               'CellDoubleClickedFcn', ...
               @(~, ~)this.cbEditButton(), ...
               'ColumnName', ...
               { getString(message('stats:dfittool:header_plot')), ...
               getString(message('stats:dfittool:header_bounds')), ...
               getString(message('stats:dfittool:header_name')), ...
               getString(message('stats:dfittool:header_setName')), ...
               getString(message('stats:dfittool:header_distribution'))}, ...
               'SelectionType', 'row', ...
               'RowName', {}, ...
               'RowStriping', 'off', ...
               'ColumnEditable', [true, true, false, false, false]);
           this.FitTable.Layout.Row = 2;
           this.FitTable.Layout.Column = 1;
           
           % Add the fits manager buttons
           iCreateDummyMinSizeButton(buttonGrid1, 1, 2);
           this.NewFitButton = uibutton(buttonGrid1, 'Text', ...
               getString(message('stats:dfittool:button_new')), ...
               'Tag', 'dfFitsManagerNewFitButton', ...
               'ButtonPushedFcn', @(~, ~)this.cbNewFitButton());
           this.NewFitButton.Layout.Row = 1;
           this.NewFitButton.Layout.Column = 2;
           
           iCreateDummyMinSizeButton(buttonGrid1, 1, 3);
           this.CopyButton = uibutton(buttonGrid1, 'Text', ...
               getString(message('stats:dfittool:button_copy')), ...
               'Tag', 'dfFitsManagerCopyButton', ...
               'ButtonPushedFcn', @(~, ~)this.cbCopyButton(), ...
               'Enable', 'off');
           this.CopyButton.Layout.Row = 1;
           this.CopyButton.Layout.Column = 3;
           
           iCreateDummyMinSizeButton(buttonGrid1, 1, 4);
           this.EditButton = uibutton(buttonGrid1, 'Text', ...
               getString(message('stats:dfittool:button_edit')), ...
               'Tag', 'dfFitsManagerEditButton', ...
               'ButtonPushedFcn', @(~, ~)this.cbEditButton(), ...
               'Enable', 'off');
           this.EditButton.Layout.Row = 1;
           this.EditButton.Layout.Column = 4;
           
           iCreateDummyMinSizeButton(buttonGrid1, 1, 5);
           this.SaveToWorkspaceButton = uibutton(buttonGrid1, 'Text', ...
               getString(message('stats:dfittool:button_save')), ...
               'Tag', 'dfFitsManagerSaveToWorkspaceButton', ...
               'ButtonPushedFcn', @(~, ~)this.cbSaveToWorkspaceButton(), ...
               'Enable', 'off');
           this.SaveToWorkspaceButton.Layout.Row = 1;
           this.SaveToWorkspaceButton.Layout.Column = 5;
           
           iCreateDummyMinSizeButton(buttonGrid1, 1, 6);
           this.DeleteButton = uibutton(buttonGrid1, ...
               'Text', getString(message('stats:dfittool:button_delete')), ...
               'Tag', 'dfFitsManagerDeleteButton', ...
               'ButtonPushedFcn', @(~, ~)this.cbDeleteButton(), ...
               'Enable', 'off');
           this.DeleteButton.Layout.Row = 1;
           this.DeleteButton.Layout.Column = 6;
           
           iCreateDummyMinSizeButton(buttonGrid2, 1, 1);
           b = uibutton(buttonGrid2, ...
               'Text', getString(message('stats:dfittool:button_help')), ...
               'Tag', 'dfFitsManagerHelpButton', ...
               'ButtonPushedFcn', @(~, ~)this.cbHelpButton());
           b.Layout.Row = 1;
           b.Layout.Column = 1;
           
           iCreateDummyMinSizeButton(buttonGrid2, 1, 3);
           b = uibutton(buttonGrid2, ...
               'Text', getString(message('stats:dfittool:button_close')), ...
               'Tag', 'dfFitsManagerCloseButton', ...
               'ButtonPushedFcn', @(~, ~)this.cbCloseButton());
           b.Layout.Row = 1;
           b.Layout.Column = 3;
       end
       
       function cbManageFitsTableCellSelection(this, eventData)
           this.SelectedFits = stats.internal.dfit.getSelectedObjects(...
               stats.internal.dfit.getfitdb, this.FitTable, ...
               eventData.Indices, 3);

           this.enableButtons();
       end
       
       function cbManageFitsTableCellEdited(this, eventData)  
            % Select the row being edited
            this.FitTable.Selection = eventData.Indices(1);
            selectedFit = stats.internal.dfit.getSelectedObjects(...
               stats.internal.dfit.getfitdb, this.FitTable, ...
               eventData.Indices, 3);
           
            if eventData.Indices(2) == 1 % plot selected
                if selectedFit.isgood
                    selectedFit.plot = eventData.EditData;
                else                   
                    this.updateFitsTable();
                    uialert(this.Panel, ...
                        getString(message('stats:dfstrings:badFit_CannotPlot')), ...
                        getString(message('stats:dfstrings:dlg_DistributionFittingMessage')),...
                        'Icon', 'info');
                end      
            else % confbounds selected
                if selectedFit.isgood
                    if eventData.EditData
                        stats.internal.dfit.boundwarn(selectedFit);
                    end
                    selectedFit.showbounds = eventData.EditData;
                else
                    this.updateFitsTable();
                    uialert(this.Panel, ...
                        getString(message('stats:dfstrings:badFit_CannotDisplayConfidenceBounds')), ...
                        getString(message('stats:dfstrings:dlg_DistributionFittingMessage')),...
                        'Icon', 'info');
                end
            end
        end

       function enableButtons(this)
           % New fit should always be enabled (so is not included here)
           % Delete fit should be enabled if at least one fit is selected.
           % Copy, Edit should be enabled if one and only one fit is selected.
           % Save to Workspace should be enabled if one and only one fit is
           % selected AND that fit's 'isgood' property is true.
           switch numel(this.SelectedFits)
               case 0
                   this.CopyButton.Enable = 'off';
                   this.EditButton.Enable = 'off';
                   this.SaveToWorkspaceButton.Enable = 'off';
                   this.DeleteButton.Enable = 'off';
               case 1
                   this.CopyButton.Enable = 'on';
                   this.EditButton.Enable = 'on';
                   if (this.SelectedFits(1).isgood)
                       this.SaveToWorkspaceButton.Enable = 'on';
                   else
                       this.SaveToWorkspaceButton.Enable = 'off';
                   end
                   this.DeleteButton.Enable = 'on';
               otherwise
                   this.CopyButton.Enable = 'off';
                   this.EditButton.Enable = 'off';
                   this.SaveToWorkspaceButton.Enable = 'off';
                   this.DeleteButton.Enable = 'on';
           end
       end
       
       function cbDeleteButton(this)
           this.deleteFits(this.SelectedFits);
       end
       
       function cbCloseButton(this)
           this.Panel.Visible = 'off';
       end
       
       function cbHelpButton(~)
           stats.internal.dfit.helpviewer('manage_fits', 'manage fits');
       end
   end
end

function displayName = iGetDistributionDisplayName(this, f)
    if strcmp(f.fittype, 'param')
        distributionCodeName = f.distname; 
    else
        distributionCodeName = 'nonparametric';
    end
    displayName = this.DistributionFitting.getFitTypeDisplayName(distributionCodeName);
end

function iUpdateSelectedFitsList(this, fit)
    for i = 1:length(this.SelectedFits)
        if strcmp(this.SelectedFits(i).name, fit.name)
            this.SelectedFits(i) = [];
            break;
        end
    end
end

function tf = iIsSelectedFits(this, name)
    tf = false;
    for i = 1:length(this.SelectedFits)
        if strcmp(this.SelectedFits(i).name, name)
            tf = true;
            return;
        end
    end
end

function iCreateDummyMinSizeButton(grid, row, column)
    dummyText = repmat('X',1,8);
    minButton = uibutton(grid, 'Text', dummyText);
    minButton.Layout.Row = row;
    minButton.Layout.Column = column;
    minButton.Visible = 'off';
end