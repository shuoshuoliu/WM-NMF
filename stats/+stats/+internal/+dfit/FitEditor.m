classdef FitEditor < handle
    % FitEditor   Dialog to edit fits
    
    %   Copyright 2019-2020 The MathWorks, Inc.

   properties(Access = private)
      DistributionFitting;   
      Panel;
      FitNameEditField;
      FitName = '';
      DataDropDown;
      DistributionDropDown;
      ExclusionRuleDropDown;
      ResultsTextArea;
      ParamPanels = {};
      CurrentParamPanel;
      CurrentParamPanelObj;
      Fit = [];
      ApplyButton;
      SaveToWorkspaceButton;
      
      % MainGrid is the parent of the parameter panels. 
      % Current parameter panel is 'Visible' 'on'.
      % All other should be 'Visible', 'off'.
      MainGrid;
      ParameterPanelGridRow = 2;
      ParameterPanelGridColumn = 1; 
      IsFitChanged = false;
      DataSetAddedListener;
      DataSetRenamedListener;
      DataSetDeletedListener;
      ExclusionRulesDeletedListener;
      ExclusionRuleAddedListener;
      ExclusionRuleRenamedListener;
      
      SessionClearingListener
   end
                   
   methods
       
	  function this = FitEditor(varargin)
          this.DistributionFitting = stats.internal.dfit.DistributionFitting.getInstance;
          if (nargin == 0)  %new fit
              this.FitName = stats.internal.dfit.getfitname();
              this.setupPanel();
          else % From copy fit or load session
              copiedFit = varargin{1};
              this.Fit = copiedFit;
              this.FitName = copiedFit.name;
              copiedFit.fitframe = this;
              this.setupPanel();
              this.updateCopiedPanel(copiedFit);
          end
          this.init();
       end

       function this = showPanel(this)
           figure(this.Panel); % Brings the uifigure forward and makes it visible.
       end
       
       function setPanelName(this, name)
           this.Panel.Name = name;
       end
       
       function fit = getFit(this)
           fit = this.Fit;
       end
       
       function fitName = getFitName(this)
           fitName = this.FitName;
       end
       
       function dataName = getDataName(this)
           dataName = this.DataDropDown.Value;
       end
      
       function exclusionRule = getExclusionRule(this)
           if strcmp(this.ExclusionRuleDropDown.Value, iNone())
               exclusionRule = '';
           else
               exclusionRule = this.ExclusionRuleDropDown.Value;
           end
       end  
       
       function setResults(this, resultsText)
           this.ResultsTextArea.Value = resultsText;
       end
                   
       function clearResults(this)
           if strcmp(this.ApplyButton.Enable, 'on')
               this.ResultsTextArea.Value = getString(message('stats:dfittool:label_fitPressApply'));
           else
               this.ResultsTextArea.Value = '';
           end
           this.SaveToWorkspaceButton.Enable = 'off';
       end
       
       function delete(this)
           delete(this.Panel);
       end

	   function setFitChanged(this)
	       this.IsFitChanged = true;
	   end
   end
   
   methods(Access = private)
       
       function updateCopiedPanel(this, newUDDObj)
           
           % "Copy fit" copies information from the uddobj as it was when 
           % last "applied". Thus, it does not necessarily match current 
           % GUI values.
           
           this.Panel.Name = getString(message('stats:dfittool:title_editFit'));

           this.DataDropDown.Value = newUDDObj.dataset;
           
           if strcmp(newUDDObj.fittype, 'param')
               distribution = newUDDObj.distname;
           else
               distribution = 'nonparametric';
           end
           
           this.DistributionDropDown.Value = this.DistributionFitting.getFitTypeDisplayName(distribution);
             
           if ~isempty( newUDDObj.exclusionrulename)
                this.ExclusionRuleDropDown.Value = newUDDObj.exclusionrulename;
           end

           pp = this.findParameterPanelPanelFromName(this.DistributionDropDown.Value);
           ppObj = this.findParameterPanelClassObjFromName(this.DistributionDropDown.Value);
           
           this.updateParamPanel(pp, ppObj);
           ppObj.updateFixedValues(newUDDObj);

           this.ResultsTextArea.Value = newUDDObj.resultstext;
           fm = stats.internal.dfit.FitsManager.getInstance;
            
           % enable saveToWorkspace
           if newUDDObj.isgood
               this.SaveToWorkspaceButton.Enable = 'on';
           else
               this.SaveToWorkspaceButton.Enable = 'off';
           end
           
           % The only time the apply button is disabled is when there are
           % no data sets. This happens only when pressing the "New Fit"
           % button when no data sets have been created. Users can change
           % data sets, but can't set the "Data" field to empty. Thus a fit
           % that can be copied must have had a dataset. So we should
           % enable the "Apply" button.

           this.ApplyButton.Enable = 'on';
           
           fm.fitAdded();
       end     
              
       function setupPanel(this)
           % The order here is important. Layout the fit panel first in 
           % order to get the "parent grid" for the parameter panels.
           
           % In order to avoid going through all the fittypes twice (once
           % to populate the distribution drop down and once to setup the
           % parameter panels), create the parameter panels first, and
           % then use the information therein to populate the distribution 
           % drop down.  
           
           this.Panel = this.layoutPanel();
           this.createParameterPanels();
           
           % List the compatible distributions
           this.listCompatibleDistributions();
           this.DistributionDropDown.Value = getString(message('stats:dfittool:NameNormal'));
                      
           % Set "Normal" as the default parameter panel.
           pp = findParameterPanelPanelFromName(this, getString(message('stats:dfittool:NameNormal')));
           this.CurrentParamPanel = pp;
           pp.Visible = 'on';
           ppObj = findParameterPanelClassObjFromName(this, getString(message('stats:dfittool:NameNormal')));
           this.CurrentParamPanelObj = ppObj;         
           this.clearResults();
       end

       function this = init(this)
           this.DataSetAddedListener = event.listener(stats.internal.dfit.Data.getInstance, 'DataSetAdded', @(~, ~) this.handleDataSetAdded);
           this.DataSetDeletedListener = event.listener(stats.internal.dfit.Data.getInstance, 'DataSetsDeleted', @(~, ~) this.handleDataSetDeleted);
           this.DataSetRenamedListener = event.listener(stats.internal.dfit.Data.getInstance, 'DataSetRenamed', @(~, evt) this.handleDataSetRenamed(evt));
           
           this.ExclusionRuleAddedListener = event.listener(stats.internal.dfit.Exclude.getInstance, 'ExclusionRuleAdded', @(~, ~) this.handleExclusionRuleAdded);
           this.ExclusionRulesDeletedListener = event.listener(stats.internal.dfit.Exclude.getInstance, 'ExclusionRulesDeleted', @(~, ~) this.handleExclusionRulesDeleted);
           this.ExclusionRuleRenamedListener = event.listener(stats.internal.dfit.Exclude.getInstance, 'ExclusionRuleRenamed', @(~, evt) this.handleExclusionRuleRenamed(evt));
           
           this.SessionClearingListener = event.listener(this.DistributionFitting, 'SessionClearing', @(~, ~)this.sessionClearing);
       end
       
       function sessionClearing(this)
           this.delete();
       end
                      
       function handleExclusionRuleAdded(this, ~)
           % Just update the list - the Value will be stay the same.  
           this.updateExclusionRuleDropDown();
       end
       
       function updateExclusionRuleDropDown(this)
           this.ExclusionRuleDropDown.Items = [{iNone()}, stats.internal.dfit.getDataBaseNames(stats.internal.dfit.getoutlierdb)];
       end
       
       function handleExclusionRuleRenamed(this, evt)
           % First update the dropdown; make sure to preserve the selected
           % exclusion rule;
           selectedExclusionRenamed = false;
           if strcmp(this.ExclusionRuleDropDown.Value, evt.OldName)
               selectedExclusionRenamed = true;
           end
           
           % update the dropdown
           this.updateExclusionRuleDropDown();
           
           % reset the value if the selected dataset was renamed
           if selectedExclusionRenamed
               this.ExclusionRuleDropDown.Value = evt.NewName;
           end
       end
       
       function handleExclusionRulesDeleted(this, ~)
           this.updateExclusionRuleDropDown();
       end
       
       function handleDataSetRenamed(this, renameEvent)           
           % First update the dropdown; make sure to preserve the selected
           % data set;
           selectedDataSetRenamed = false;
           if strcmp(this.DataDropDown.Value, renameEvent.OldName)
               selectedDataSetRenamed = true;
           end
           
           % update the dropdown
           this.DataDropDown.Items = stats.internal.dfit.getDataBaseNames(stats.internal.dfit.getdsdb);
           
           % reset the value if the selected dataset was renamed
           if selectedDataSetRenamed
               this.DataDropDown.Value = renameEvent.NewName;
           end
       end
       
       function handleDataSetAdded(this)
           % Just update the list - the Value will be stay the same.  
           this.DataDropDown.Items = stats.internal.dfit.getDataBaseNames(stats.internal.dfit.getdsdb);
           
           % If we were going from no data sets to one data set, enable
           % apply
           dsdb = stats.internal.dfit.getdsdb;
           if (stats.internal.dfit.getNumDataBaseItems(dsdb) == 1)
               this.ApplyButton.Enable = 'on';
               this.clearResults();  % should go from no message to "press apply" message
           end
       end
       
       function handleDataSetDeleted(this)
           dsNames = stats.internal.dfit.getDataBaseNames(stats.internal.dfit.getdsdb);
           if isempty(dsNames)
               this.DataDropDown.Items = {};
               this.ApplyButton.Enable = 'off';
               this.clearResults();
           else
               this.DataDropDown.Items = dsNames;
           end
       end
              
       function panel = layoutPanel(this)
           panel = uifigure('Tag', 'dfFitEditorUIFigure', ...
               'Name', getString(message('stats:dfittool:title_newFit')), ...
               'Position', [700 350 430 720], ...
               'CloseRequestFcn', @(~, ~)this.cbCloseButton(), 'Visible', 'off');

           this.MainGrid = uigridlayout(panel, [6, 1]);
           this.MainGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', '1x', 'fit'};
           this.MainGrid.ColumnWidth = {'1x'};
           
           fitInfoGrid = uigridlayout(this.MainGrid, [4, 2]);
           fitInfoGrid.RowHeight = {'1x'};
           fitInfoGrid.ColumnWidth = {'fit', '1x'};
           fitInfoGrid.Layout.Row = 1;
           fitInfoGrid.Layout.Column = 1;
           
           applyGrid = uigridlayout(this.MainGrid, [1, 2]);
           applyGrid.RowHeight = {'fit'};
           applyGrid.ColumnWidth = {'1x', 'fit'};
           applyGrid.Layout.Row = 3;
           applyGrid.Layout.Column = 1;
           bottomButtonGrid = uigridlayout(this.MainGrid, [1, 6]);
           bottomButtonGrid.RowHeight = {'fit'};
           bottomButtonGrid.ColumnWidth = {'fit', '1x' 'fit',  'fit',  '1x', 'fit'};
           bottomButtonGrid.Layout.Row = 6;
           bottomButtonGrid.Layout.Column = 1;
           
           % Populate the fitInfoGrid - labels (first column)
           l = uilabel(fitInfoGrid, 'Text', getString(message('stats:dfittool:label_fitName')));
           l.Layout.Row = 1;
           l.Layout.Column = 1;
           l = uilabel(fitInfoGrid, 'Text', getString(message('stats:dfittool:label_data')));
           l.Layout.Row = 2;
           l.Layout.Column = 1;
           l = uilabel(fitInfoGrid, 'Text', getString(message('stats:dfittool:label_distribution')));
           l.Layout.Row = 3;
           l.Layout.Column = 1;
           l = uilabel(fitInfoGrid, 'Text', getString(message('stats:dfittool:label_exclusionRule')));
           l.Layout.Row = 4;
           l.Layout.Column = 1;

           this.FitNameEditField = uieditfield(fitInfoGrid, 'Tag', 'dfFitEditorFitNameEditField', 'Value', this.FitName);
           this.FitNameEditField.Layout.Row = 1;
           this.FitNameEditField.Layout.Column = 2;
           
           dataSets = stats.internal.dfit.getDataBaseNames(stats.internal.dfit.getdsdb);
           
           this.DataDropDown = uidropdown(fitInfoGrid, 'Tag', 'dfFitEditorDataDropDown', ...
               'Items', dataSets, 'ValueChangedFcn', @(~, ~)this.cbDataOrExclusionRulesChanged());
           this.DataDropDown.Layout.Row = 2;
           this.DataDropDown.Layout.Column = 2;
            
           this.DistributionDropDown = uidropdown(fitInfoGrid, ...
                'Tag', 'dfFitEditorDistributionDropDown', ...
                'ValueChangedFcn', @(src, ~)this.cbDistributionChanged(src));
           this.DistributionDropDown.Layout.Row = 3;
           this.DistributionDropDown.Layout.Column = 2;
           
           this.ExclusionRuleDropDown = uidropdown(fitInfoGrid, ...
               'Tag', 'dfFitEditorExclusionRuleDropDown', 'Items', {iNone()}, ...
               'ValueChangedFcn', @(~, ~)this.cbDataOrExclusionRulesChanged());
           this.updateExclusionRuleDropDown();
           this.ExclusionRuleDropDown.Layout.Row = 4;
           this.ExclusionRuleDropDown.Layout.Column = 2;
                                
           % Add the Apply button
           iCreateDummyMinSizeButton(applyGrid, 1, 2);
           this.ApplyButton = uibutton(applyGrid, ...
               'Tag', 'dfFitEditorApplyButton', ...
               'Enable', 'off', ...
               'Text', getString(message('stats:dfittool:button_apply')), ...
               'ButtonPushedFcn', @(~, ~)this.cbApplyButton());
           this.ApplyButton.Layout.Row = 1;
           this.ApplyButton.Layout.Column = 2;
           
           % Enable the apply button if there are any data sets
           dsdb = stats.internal.dfit.getdsdb;
           if (stats.internal.dfit.getNumDataBaseItems(dsdb) > 0)
               this.ApplyButton.Enable = 'on';
               this.clearResults(); 
           end
           
           % Add the Results label
           l = uilabel(this.MainGrid, 'Text', getString(message('stats:dfittool:title_results')));
           l.Layout.Row = 4;
           l.Layout.Column = 1;
           
           % Add the Results text area
           this.ResultsTextArea = uitextarea(this.MainGrid, ...
               'Tag', 'dfFitEditorResultsTextArea');
           this.ResultsTextArea.Layout.Row = 5;
           this.ResultsTextArea.Layout.Column = 1;
           this.ResultsTextArea.Value = getString(message('stats:dfittool:label_fitPressApply'));
           
           % Populate the bottom button grid
           iCreateDummyMinSizeButton(bottomButtonGrid, 1, 1);
           b = uibutton(bottomButtonGrid, ...
               'Text', getString(message('stats:dfittool:button_help')), ...
               'Tag', 'dfFitEditorHelpButton', ...    
               'ButtonPushedFcn', @(~, ~)this.cbHelpButton());
           b.Layout.Row = 1;
           b.Layout.Column = 1;
           
           iCreateDummyMinSizeButton(bottomButtonGrid, 1, 3);
           this.SaveToWorkspaceButton = uibutton(bottomButtonGrid, ...
               'Tag', 'dfFitEditorSaveToWorkspaceButton', ...  
               'Text', getString(message('stats:dfittool:button_save')),...
               'Enable', 'off', ... 
               'ButtonPushedFcn', @(~, ~)this.cbSaveToWorkspaceButton());
           this.SaveToWorkspaceButton.Layout.Row = 1;
           this.SaveToWorkspaceButton.Layout.Column = 3;
           
           iCreateDummyMinSizeButton(bottomButtonGrid, 1, 4);
           b = uibutton(bottomButtonGrid, ...
               'Tag', 'dfFitEditorManageFitsButton', ...  
               'Text', getString(message('stats:dfittool:button_manageFits')), ...
               'ButtonPushedFcn', @(~, ~)this.cbManageFits());
           b.Layout.Row = 1;
           b.Layout.Column = 4;
           
           iCreateDummyMinSizeButton(bottomButtonGrid, 1, 6);
           b = uibutton(bottomButtonGrid, ...
               'Tag', 'dfFitEditorCloseButton', ... 
               'Text', getString(message('stats:dfittool:button_close')), ...
               'ButtonPushedFcn', @(~, ~)this.cbCloseButton());
           b.Layout.Row = 1;
           b.Layout.Column = 6;
                      
           stats.internal.dfit.centergui(panel); 
           panel.Visible = 'on';
       end     
       
       function cbSaveToWorkspaceButton(this, ~)
            stats.internal.dfit.save2ws(this.FitNameEditField.Value())
       end
       
       function cbManageFits(~, ~)
           fm = stats.internal.dfit.FitsManager.getInstance;
           fm.showPanel;
       end
                     
       function cbDistributionChanged(this, src)
           this.clearResults();
           this.SaveToWorkspaceButton.Enable = 'off';
           pp = findParameterPanelPanelFromName(this, src.Value);
           ppObj = findParameterPanelClassObjFromName(this, src.Value);
           this.updateParamPanel(pp, ppObj);
           this.IsFitChanged = true;
       end
       
       function cbDataOrExclusionRulesChanged(this)
           this.clearResults();
           this.SaveToWorkspaceButton.Enable = 'off';
           this.listCompatibleDistributions();
           this.IsFitChanged = true;
       end
       
       function pp = findParameterPanelPanelFromName(this, name)
           % findParameterPanelPanelFromName returns the actual panel for a 
           % given distribution
           pp = [];
           for i = 1:length(this.ParamPanels)
                if strcmp(name, this.ParamPanels{i}.getDisplayName())
                    pp = this.ParamPanels{i}.getParameterPanel();
                    break;
                end
           end

           % pp should never be empty.
           if isempty(pp)
                assert(false, getString(message('stats:dfittool:assert_findParameterPanelPanelFromName')));
           end
       end
       
       function ppObj = findParameterPanelClassObjFromName(this, name)
           % findParamPanelClassObjFromName returns either the 
           % ParametricParameterPanel or the NonParametricParameterPanel class object 
           ppObj = [];
           for i = 1:length(this.ParamPanels)
                if strcmp(name, this.ParamPanels{i}.getDisplayName())
                    ppObj = this.ParamPanels{i};
                    break;
                end
           end

           % ppObj should never be empty.
           if isempty(ppObj)
               assert(false, getString(message('stats:dfittool:assert_findParameterPanelClassObjFromName')));
           end
       end
                      
       function updateParamPanel(this, pp, ppObj) 
           % Make the current panel invisible.
           this.CurrentParamPanel.Visible = 'off';
           % Make the new parameter panel visible and set it as the 
           % current panel.
           pp.Visible = 'on';
           this.CurrentParamPanel = pp;
           this.CurrentParamPanelObj = ppObj;
       end
       
       function cbApplyButton(this)           
           if iIsNameInvalid(this)
               return;
           end
           
           if this.hasInvalidValues()
               return;
           end
           
           oldName = this.FitName;
           trimmedName = strtrim(this.FitNameEditField.Value);
           this.FitNameEditField.Value = trimmedName;
           fm = stats.internal.dfit.FitsManager.getInstance;
           
           newFit = isNewFit(this.Panel.Name);
           
           % If only the name has changed, don't bother with anything
           % else.
           if ~(this.IsFitChanged) && ~newFit
               % no changes except perhaps a name change
               if strcmp(this.FitName, trimmedName)
                   return;
               else
                   this.FitName = trimmedName;
                   this.Fit.name = trimmedName;
                   fm.fitChanged(this.Fit, this.Panel, oldName, trimmedName);
                   return;
               end
           end
           
           this.FitName = trimmedName;
           
           fitType = getSelectedFitType(this.DistributionFitting, this.DistributionDropDown.Value);
           
           ppObj = this.CurrentParamPanelObj;
           
           fitArgs = ppObj.getFitArgs(fitType, this);
           
           this.Fit = distributionFitter(fitArgs{:});
           
           this.SaveToWorkspaceButton.Enable = 'off';
           
           if ~isempty(this.Fit)  % fit added or changed
               if newFit
                   fm.fitAdded();
               else
                   fm.fitChanged(this.Fit, this.Panel, this.FitName, this.FitName);
               end
               if this.Fit.isgood == true
                   this.SaveToWorkspaceButton.Enable = 'on';
               end
           end
           
           % reset fit changed flag
           this.IsFitChanged = false;
       end         
       
       function tf = hasInvalidValues(this)
           ppObj = this.CurrentParamPanelObj;
           tf = ppObj.isBadValue();
       end
       
       function cbCloseButton(this)
           if isNewFit(this.Panel.Name)
               this.delete();
           else
               this.Panel.Visible = 'off';
           end
       end
       
       function cbHelpButton(~)
           stats.internal.dfit.helpviewer('new_fit', 'fits');       
       end
          
       function createParameterPanels(this)
           fitTypes = this.DistributionFitting.getFitTypes();
           this.ParamPanels = cell(1, length(fitTypes));
           for i = 1:length(fitTypes)
               ft = fitTypes{i};        
               if strcmp(ft.getDisplayName(), getString(message('stats:dfittool:NameNonparametric')))
                   this.ParamPanels{i} = stats.internal.dfit.NonParametricParameterPanel(this, this.MainGrid, this.ParameterPanelGridRow, this.ParameterPanelGridColumn, ft);
               else
                   this.ParamPanels{i} = stats.internal.dfit.ParametricParameterPanel(this, this.MainGrid,  this.ParameterPanelGridRow, this.ParameterPanelGridColumn, ft);
               end
           end
       end  

       function listCompatibleDistributions(this)
           if isempty(this.DataDropDown.Value)  %No data set specified
               isint = true;
               iscens = false;
               datahi = -Inf;
               datalo = Inf;
           elseif isequal(this.ExclusionRuleDropDown.Value, iNone())  %Data, but no exclusion rule
               db = find(stats.internal.dfit.getdsdb,'name', this.DataDropDown.Value);
               isint = db.isinteger;
               iscens = db.iscensored;
               datalo = db.datalo;
               datahi = db.datahi;
           else % exclusion rule specified
               [isint,iscens,datalo,datahi] = this.getIncludedDatastats();
           end
           oldDistribution = this.DistributionDropDown.Value;
           newDistribution = getString(message('stats:dfittool:NameNormal'));
           j = 1;
           for i = 1:length(this.ParamPanels)
               ft = this.ParamPanels{i}.getFitType();
               % If its compatible, add it, otherwise continue
               if (ft.getIntegerOnly() && ~isint)
                   continue;
               end
               if (iscens && ~ft.getCensoredOK())
                   continue;
               end
               ll = ft.getLowerLimit();
               if (ft.getLowerLimitOK())
                   if (datalo < ll)
                       continue;
                   end
               else
                   if (datalo <= ll)
                       continue;
                   end
               end
               ul = ft.getUpperLimit();
               if (ft.getUpperLimitOK())
                   if (datahi > ul)
                       continue;
                   end
               else
                   if (datahi >= ul)
                       continue;
                   end
               end
               % If we get here, this distribution is OK. Add it to the
               % list, check to see if this was the distribution selected
               % before, if so, save it ...
               distributionDropDownItems{j} = this.ParamPanels{i}.getDisplayName(); %#ok<AGROW>
               j = j+ 1;
               if isequal(this.ParamPanels{i}.getDisplayName(), oldDistribution)
                   newDistribution = oldDistribution;
               end
           end
           this.DistributionDropDown.Items = distributionDropDownItems;          
           this.DistributionDropDown.Value = newDistribution;
       end
       
       function [isint,iscens,datalo,datahi] = getIncludedDatastats(this)
           % GETINCLUDEDDATASTATS Get stats on included data, to determine available fits          
           
           exclname = this.ExclusionRuleDropDown.Value;
           dsname = this.DataDropDown.Value; 
           % Return defaults in case the unexpected happens
           isint = true;
           iscens = false;
           datalo = Inf;
           datahi = -Inf;
           
           % Get exclusion rule by name
           if ~isempty(exclname)
               excl = find(stats.internal.dfit.getoutlierdb, 'name', exclname);
           else
               excl = [];
           end
           
           % Get data set
           if ~isempty(dsname)
               ds = find(stats.internal.dfit.getdsdb,'name',dsname);
           else
               ds = [];
           end
           if ~isempty(ds)
               if isempty(excl)
                   % If there's no exclusion rule, return information from data set
                   isint  = ds.isinteger;
                   iscens = ds.iscensored;
                   datalo = ds.datalo;
                   datahi = ds.datahi;
               else
                   % Compute information from included data
                   [y,cens] = getincludeddata(ds,excl);
                   if ~isempty(cens)
                       y(cens==1) = [];
                       iscens = any(cens(:));
                   end
                   y = y(:);
                   if ~isempty(y)
                       isint = all(y == round(y));
                       datalo = min(y);
                       datahi = max(y);
                   end
               end
           end
       end
   end 
end

% helpers
function str = iNone()
str = strcat('(',getString(message('stats:dfittool:label_none')),')');
end

function tf = iIsNameInvalid(this)
    tf = false;
    
    fitMessage = getString(message('stats:dfittool:label_fit'));
    fitTitle = getString(message('stats:dfittool:label_fitUpperCase'));
    
    fitName = this.FitNameEditField.Value;

    % empty name?
    if strcmp(fitName, "")
        uialert(this.Panel, ...
            getString(message('stats:dfittool:error_noName', fitMessage)), ...
            getString(message('stats:dfittool:title_noObjectTypeName', fitTitle)));
        tf = true;
        return;
    end

    % just spaces?
    trimmedName = strtrim(fitName);
    if strcmp(trimmedName, "")
        uialert(this.Panel, ...
            getString(message('stats:dfittool:error_noOnlySpacesAllowed')), ...
            getString(message('stats:dfittool:title_invalidObjectTypeName', fitTitle)));
        tf = true;
        return;
    end

    % fit name the same?
    f = down(stats.internal.dfit.getfitdb);
    while(~isempty(f))
        if ~isequal(f, this.Fit) % Don't count the current item
            if isequal(f.name, trimmedName)
                uialert(this.Panel, ...
                    getString(message('stats:dfittool:error_duplicateFitName',  ...
                    trimmedName)), ...
                    getString(message('stats:dfittool:title_duplicateName', fitTitle)));
                tf = true;
                break;
            end
        end
        f = right(f);
    end
end

function TF = isNewFit(name)
    if isequal(name, getString(message('stats:dfittool:title_newFit')))
        TF = true;
    else
        TF = false;
    end
end

function iCreateDummyMinSizeButton(grid, row, column)
    dummyText = repmat('X',1,8);
    minButton = uibutton(grid, 'Text', dummyText);
    minButton.Layout.Row = row;
    minButton.Layout.Column = column;
    minButton.Visible = 'off';
end