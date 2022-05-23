classdef Exclude < handle
    % Exclude   Dialog to create rules for excluding specified data values
    
    %   Copyright 2019-2020 The MathWorks, Inc.
    
    properties(Access = private)
        % DistributionFitting   (stats.internal.dift.DistributionFitting)
        DistributionFitting
        
        % ExclusionRuleDatabase   (stats.outlierdb)
        ExclusionRuleDatabase
        
        % Data   (stats.internal.dfit.Data)
        Data
        
        % ExcludePanel   (stats.internal.dfit.ExcludePanel)
        ExcludePanel
        
        % RenamePanel   (stats.internal.dfit.Rename)
        RenamePanel
        
        % ExclusionRuleConstraint
        % (stats.internal.dfit.ExclusionRuleConstraint)
        ExclusionRuleConstraint
        
        % SelectedData   (string) name of selected data set used to create
        % exclusion rule
        SelectedData
        
        % SelectedExclusionRules   (cell of strings) name(s) of exclusion
        % selected rules object 
        SelectedExclusionRules
        
        % Listeners   (cell of listener) Listeners to various events
        Listeners = {};
    end
    
    events
        ExclusionRulesDeleted
        ExclusionRuleAdded
        ExclusionRuleRenamed
    end
    
    methods(Static)
        function obj = getInstance()
            persistent Instance;
            if ~isa(Instance, 'stats.internal.dfit.Exclude') || ~isvalid(Instance)
                Instance = stats.internal.dfit.Exclude();
            end
            obj = Instance;
        end
    end
    
    methods
        function showPanel(this)
            this.ExcludePanel.showPanel();
        end
        
        function deleteExclusionRule(this, exclusionRuleNames)
            if iIsOKToDeleteExclusionRules(exclusionRuleNames)
                cellfun(@(x) iDeleteExclusionRule(x), exclusionRuleNames);
                
                newExclusionRulesList = iExistingExclusionRulesList();
                this.ExcludePanel.updateExistingExclusionRulesList(newExclusionRulesList);
                this.SelectedExclusionRules = intersect(this.SelectedExclusionRules, newExclusionRulesList);
                this.ExcludePanel.updateExistingExclusionRulesValue(this.SelectedExclusionRules);
                this.updateExistingExclusionRulesButtonStates();
                
                this.notify('ExclusionRulesDeleted');
            end 
        end
        
        function updateBoundsFields(this, lowerLimitStrValue, upperLimitStrValue)
            this.ExcludePanel.updateLowerLimitTextFieldValue(lowerLimitStrValue);
            this.ExcludePanel.updateUpperLimitTextFieldValue(upperLimitStrValue);
        end
        
        function delete(this)
            stats.internal.dfit.delgraphexclude();
            delete(this.ExcludePanel);
            delete(this.RenamePanel);
        end
    end
    
    methods(Access = private)
        function this = Exclude()
            this.DistributionFitting = stats.internal.dfit.DistributionFitting.getInstance();
            this.Data = stats.internal.dfit.Data.getInstance();
            this.ExclusionRuleDatabase = iOutlierDatabase();
            this.ExcludePanel = stats.internal.dfit.ExcludePanel();
            this.RenamePanel = stats.internal.dfit.Rename(this.ExclusionRuleDatabase);
            this.ExclusionRuleConstraint = stats.internal.dfit.ExclusionRuleConstraint();
            
            this.initializeDataSetList();
            this.initializeExistingExclusionRulesList();
            
            this.addListener(this.DistributionFitting, 'SessionClearing', @(~,~) this.delete);
            this.addListener(this.Data, 'DataSetAdded', @this.dataSetAddedOrDeletedCallback);
            this.addListener(this.Data, 'DataSetsDeleted', @this.dataSetAddedOrDeletedCallback);
            this.addListener(this.Data, 'DataSetRenamed', @this.dataSetRenamedCallback);
            this.addListener(this.ExcludePanel, 'LowerLimitTextFieldValueChanged', @this.lowerLimitTextFieldValueChangedCallback);
            this.addListener(this.ExcludePanel, 'UpperLimitTextFieldValueChanged', @this.upperLimitTextFieldValueChangedCallback);
            this.addListener(this.ExcludePanel, 'SelectDataDropDownValueChanged', @this.selectDataDropDownValueChangedCallback);
            this.addListener(this.ExcludePanel, 'ExcludeGraphicallyButtonClicked', @this.excludeGraphicallyButtonClickedCallback);
            this.addListener(this.ExcludePanel, 'CreateExclusionRuleButtonClicked', @this.createExclusionRuleButtonClickedCallback);
            this.addListener(this.ExcludePanel, 'ExistingExclusionRulesListBoxValueChanged', @this.existingExclusionRulesListBoxValueChangedCallback);
            this.addListener(this.ExcludePanel, 'CopyButtonClicked', @this.copyButtonClickedCallback);
            this.addListener(this.ExcludePanel, 'ViewButtonClicked', @this.viewButtonClickedCallback);
            this.addListener(this.ExcludePanel, 'RenameButtonClicked', @this.renameButtonClickedCallback);
            this.addListener(this.ExcludePanel, 'DeleteButtonClicked', @this.deleteButtonClickedCallback);
            this.addListener(this.ExcludePanel, 'CloseButtonClicked', @this.closeButtonClickedCallback);
            this.addListener(this.ExcludePanel, 'HelpButtonClicked', @this.helpButtonClickedCallback);
            this.addListener(this.RenamePanel, 'RenameSuccessful', @this.renameSuccessfulCallback);
        end
        
        function initializeDataSetList(this)
            this.SelectedData = iNone();
            this.ExcludePanel.updateSelectDataDropDownList(iDataSetList());
            this.ExcludePanel.updateSelectDataDropDownValue(this.SelectedData);
            this.ExcludePanel.enableExcludeGraphicallyButton(false);
        end
        
        function initializeExistingExclusionRulesList(this)
            this.ExcludePanel.updateExistingExclusionRulesList(iExistingExclusionRulesList());
            this.updateExistingExclusionRulesButtonStates();
        end
        
        function updateExistingExclusionRulesButtonStates(this)
            switch numel(this.SelectedExclusionRules)
                case 0
                    this.ExcludePanel.enableCopyButton(false);
                    this.ExcludePanel.enableViewButton(false);
                    this.ExcludePanel.enableRenameButton(false);
                    this.ExcludePanel.enableDeleteButton(false);
                case 1
                    this.ExcludePanel.enableCopyButton(true);
                    this.ExcludePanel.enableViewButton(true);
                    this.ExcludePanel.enableRenameButton(true);
                    this.ExcludePanel.enableDeleteButton(true);
                otherwise
                    this.ExcludePanel.enableCopyButton(false);
                    this.ExcludePanel.enableViewButton(false);
                    this.ExcludePanel.enableRenameButton(false);
                    this.ExcludePanel.enableDeleteButton(true);
            end
        end
        
        function resetPanel(this)
            this.SelectedData = iNone();
            this.ExcludePanel.resetLeftPanel();
        end
        
        function createExclusionRule(this, exclusionRuleName, selectedData, lowerLimit, upperLimit, isLowerLimitLessEqual, isUpperLimitGreaterEqual)
            stats.outlier(exclusionRuleName, selectedData, lowerLimit, upperLimit, isLowerLimitLessEqual, isUpperLimitGreaterEqual);
            
            this.ExcludePanel.updateExistingExclusionRulesList(iExistingExclusionRulesList());
            this.updateExistingExclusionRulesButtonStates();
            this.notify('ExclusionRuleAdded');
        end
        
        function [newname, dataset, yl, yh, yle, yhe] = copyExclusionRule(this, exclusionRuleName)
            exruleObj = find(this.ExclusionRuleDatabase, 'name', exclusionRuleName);
            
            %create copy name
            COPY = sprintf(' %s ', iGetMessageString('stats:dfstrings:assignment_Copy'));
            index = strfind(exclusionRuleName, COPY);
            if isempty(index)
                sourcename = sprintf('%s%s', exclusionRuleName, COPY);
            else
                sourcename = sprintf('%s%s', exclusionRuleName(1:index-1), COPY);
            end
            cn = 1;
            newname = sprintf('%s%d', sourcename, cn);
            %loop until unique name is found
            while true
                if isempty(find(this.ExclusionRuleDatabase, 'name', newname))
                    break;
                else
                    cn = cn+1;
                    newname = sprintf('%s%d', sourcename, cn);
                end
            end
            
            %make sure dataset still exists
            if isempty(find(iDatasetDatabase(), 'name', exruleObj.dataset))
                dataset = iNone(); % not to be translated
            else
                dataset = exruleObj.dataset;
            end
            
            yl = exruleObj.YLow;
            yh = exruleObj.YHigh;
            yle = exruleObj.YLowLessEqual;
            yhe = exruleObj.YHighGreaterEqual;
        end
                
        function addListener(this, obj, eventName, callback)
            this.Listeners{end+1} = event.listener(obj, eventName, callback);
        end
                
        % Callbacks
        function dataSetAddedOrDeletedCallback(this, ~, ~)
            newDataSetList = iDataSetList();
            if iIsSelectedDataDeleted(this.SelectedData, newDataSetList)
                this.SelectedData = iNone();
                stats.internal.dfit.delgraphexclude();
            end
            this.ExcludePanel.updateSelectDataDropDownList(newDataSetList);
            this.ExcludePanel.updateSelectDataDropDownValue(this.SelectedData);
            this.ExcludePanel.enableExcludeGraphicallyButton(~iIsNone(this.SelectedData));
        end
        
        function dataSetRenamedCallback(this, ~, renameEventData)
            oldname = renameEventData.OldName;
            newName = renameEventData.NewName;
            
            affectedExclusionRules = find(this.ExclusionRuleDatabase, 'dataset', oldname);
            for i = 1:numel(affectedExclusionRules)
                 affectedExclusionRules(i).dataset = newName;
            end

            if strcmp(this.SelectedData, oldname)
                this.SelectedData = newName;
            end
            this.ExcludePanel.updateSelectDataDropDownList(iDataSetList());
            this.ExcludePanel.updateSelectDataDropDownValue(this.SelectedData);
        end
      
        function renameSuccessfulCallback(this, ~, renameEventData)
            newName = renameEventData.NewName;
            this.initializeExistingExclusionRulesList();
            this.SelectedExclusionRules = {newName};
            this.ExcludePanel.updateExistingExclusionRulesValue(newName);
            this.notify('ExclusionRuleRenamed', renameEventData);
        end
        
        function lowerLimitTextFieldValueChangedCallback(this, ~, ~)
            textFieldValue = strtrim(this.ExcludePanel.LowerLimitTextFieldValue);
            isLowerLimitLessEqual = this.ExcludePanel.IsLowerLimitLessEqual;
            this.ExcludePanel.updateLowerLimitTextFieldValue(textFieldValue);
            [~, errMessage, messageTitle] = this.ExclusionRuleConstraint.getLowerLimitValueAndErrorMessage(textFieldValue, isLowerLimitLessEqual);
            if ~isempty(errMessage)
                this.ExcludePanel.showMessageDialog(errMessage, messageTitle);
                return;
            end
            
            if iHasExcludeGraphicallyWindowOpened()
                this.excludeGraphicallyButtonClickedCallback();
            end
        end
        
        function upperLimitTextFieldValueChangedCallback(this, ~, ~)
            textFieldValue = strtrim(this.ExcludePanel.UpperLimitTextFieldValue);
            isUpperLimitGreaterEqual = this.ExcludePanel.IsUpperLimitGreaterEqual;
            this.ExcludePanel.updateUpperLimitTextFieldValue(textFieldValue);
            [~, errMessage, messageTitle] = this.ExclusionRuleConstraint.getUpperLimitValueAndErrorMessage(textFieldValue, isUpperLimitGreaterEqual);
            if ~isempty(errMessage)
                this.ExcludePanel.showMessageDialog(errMessage, messageTitle);
                return;
            end
            
            if iHasExcludeGraphicallyWindowOpened()
                this.excludeGraphicallyButtonClickedCallback();
            end
        end
       
        function selectDataDropDownValueChangedCallback(this, ~, ~)
            stats.internal.dfit.delgraphexclude();
            this.SelectedData = this.ExcludePanel.SelectedData;
            this.ExcludePanel.enableExcludeGraphicallyButton(~iIsNone(this.SelectedData));
        end
        
        function excludeGraphicallyButtonClickedCallback(this, ~, ~)
            lowerLimitTextFieldValue = this.ExcludePanel.LowerLimitTextFieldValue;
            upperLimitTextFieldValue = this.ExcludePanel.UpperLimitTextFieldValue;
            lowerLimitValue = this.ExclusionRuleConstraint.getLowerLimitValueAndErrorMessage(lowerLimitTextFieldValue, true);
            upperLimitValue = this.ExclusionRuleConstraint.getUpperLimitValueAndErrorMessage(upperLimitTextFieldValue, true);
            
            stats.internal.dfit.graphexclude(this.SelectedData, lowerLimitValue, upperLimitValue);
        end
        
        function createExclusionRuleButtonClickedCallback(this, ~, ~)
            exclusionRuleName = this.ExcludePanel.ExclusionRuleName;
            selectedData = this.SelectedData;
            lowerLimit = this.ExcludePanel.LowerLimitTextFieldValue;
            upperLimit = this.ExcludePanel.UpperLimitTextFieldValue;
            isLowerLimitLessEqual = this.ExcludePanel.IsLowerLimitLessEqual;
            isUpperLimitGreaterEqual = this.ExcludePanel.IsUpperLimitGreaterEqual;
            
            trimmedExclusionRuleName = strtrim(exclusionRuleName);
            this.ExcludePanel.updateExclusionRuleNameTextFieldValue(trimmedExclusionRuleName);
            
            [errMessage, messageTitle] = this.ExclusionRuleConstraint.getErrorMessage(...
                exclusionRuleName, lowerLimit, upperLimit,...
                isLowerLimitLessEqual, isUpperLimitGreaterEqual);
            
            if ~isempty(errMessage)
                this.ExcludePanel.showMessageDialog(errMessage, messageTitle);
            else
                this.createExclusionRule(trimmedExclusionRuleName, selectedData, lowerLimit, upperLimit, isLowerLimitLessEqual, isUpperLimitGreaterEqual);
                this.resetPanel();
                stats.internal.dfit.delgraphexclude();
            end
        end
        
        function existingExclusionRulesListBoxValueChangedCallback(this, ~, ~)
            this.SelectedExclusionRules = this.ExcludePanel.SelectedExclusionRules;
            this.updateExistingExclusionRulesButtonStates();
        end
        
        function copyButtonClickedCallback(this, ~, ~)
            [newRuleName, datasetName, lowerLimitValue, upperLimitValue,...
                isLowerLimitLessEqual, isUpperLimitGreaterEqual] = this.copyExclusionRule(this.SelectedExclusionRules{1});
            
            this.ExcludePanel.updateExclusionRuleNameTextFieldValue(newRuleName);
            this.ExcludePanel.updateLowerLimitDropDownValue(isLowerLimitLessEqual);
            this.ExcludePanel.updateUpperLimitDropDownValue(isUpperLimitGreaterEqual);
            this.ExcludePanel.updateLowerLimitTextFieldValue(lowerLimitValue);
            this.ExcludePanel.updateUpperLimitTextFieldValue(upperLimitValue);
            this.ExcludePanel.updateSelectDataDropDownValue(datasetName);
            this.ExcludePanel.updateExistingExclusionRulesValue({});
            
            this.selectDataDropDownValueChangedCallback();
            this.existingExclusionRulesListBoxValueChangedCallback();
        end
        
        function viewButtonClickedCallback(this, ~, ~)
            exclusionRule = find(this.ExclusionRuleDatabase, 'name', this.SelectedExclusionRules{1});
            if isempty(exclusionRule.viewExclusionRule)
                viewExclusionRuleObj = stats.internal.dfit.ViewExclusionRule(exclusionRule);
                exclusionRule.viewExclusionRule = viewExclusionRuleObj;
            end
            exclusionRule.viewExclusionRule.showPanel();
        end
        
        function renameButtonClickedCallback(this, ~, ~)
            oldName = this.ExcludePanel.SelectedExclusionRules{1};
            this.RenamePanel.showPanel(oldName);
        end

        function deleteButtonClickedCallback(this, ~, ~)
            exclusionRulesToDelete = this.SelectedExclusionRules;
            this.deleteExclusionRule(exclusionRulesToDelete);
        end
        
        function closeButtonClickedCallback(this, ~, ~)
            this.RenamePanel.hidePanel();
            this.ExcludePanel.hidePanel();
        end
        
        function helpButtonClickedCallback(~, ~, ~)
            stats.internal.dfit.helpviewer('exclude_data', 'exclude');
        end
    end
end

% helpers
function outputStr = iGetMessageString(inputStr)
outputStr = getString(message(inputStr));
end

function str = iNone()
str = strcat('(',getString(message('stats:dfittool:label_none')),')');
end

function tf = iIsNone(str)
tf = isequal(str, iNone()) || isempty(str);
end

function tf = iIsSelectedDataDeleted(selectedData, dataSetList)
tf = ~ismember(selectedData, dataSetList);
end

function databaseObj = iDatasetDatabase()
databaseObj = stats.internal.dfit.getdsdb;
end

function databaseObj = iOutlierDatabase()
databaseObj = stats.internal.dfit.getoutlierdb;
end

function names = iDataSetList()
names = [iNone(),...
    stats.internal.dfit.getDataBaseNames(iDatasetDatabase())];
end

function names = iExistingExclusionRulesList()
names = stats.internal.dfit.getDataBaseNames(iOutlierDatabase());
end

function iDeleteExclusionRule(exclusionRuleName)
exclusionRuleObj = find(iOutlierDatabase(), 'name', exclusionRuleName);
if ~isempty(exclusionRuleObj.viewExclusionRule)
    delete(exclusionRuleObj.viewExclusionRule);
end
delete(exclusionRuleObj);
end

function tf = iHasExcludeGraphicallyWindowOpened()
t = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');
c = get(0,'Child');
f = findobj(c,'flat','Type','figure','Tag','dfexcludegraph');
set(0,'ShowHiddenHandles',t);

tf = ~isempty(f);
end

function tf = iIsOKToDeleteExclusionRules(names)
fitdb = stats.internal.dfit.getfitdb;
fit = down(fitdb);
tf = true;
fitsToDelete = {};
associatedFits = {};

if ~isempty(fit)
    msg = '';
    for i=1:length(names)
        fitCnt = 0;
        fitnames = '';
        m='';
        fit = down(fitdb);
        while(~isempty(fit))
            if strcmp(names{i}, fit.exclusionrulename)
                fitsToDelete{1, end + 1} = fit.name; %#ok<*AGROW>
                associatedFits{1, end + 1} = fit;
                if fitCnt > 0
                    fitnames = [fitnames, ', '];
                end
                fitCnt = fitCnt + 1;
                fitnames = [fitnames, fit.name];
            end
            fit = right(fit);
        end
        if fitCnt == 1
            m = getString(message('stats:dfstrings:sprintf_FitWillBeDeleted', names{i}, fitnames));
        elseif fitCnt > 1
            m = getString(message('stats:dfstrings:sprintf_FitsWillBeDeleted', names{i}, fitnames));
        end
        msg = [msg, m]; 
    end
    if ~isempty(msg)
        button = questdlg(msg, getString(message('stats:dfstrings:dlg_DeletingExclusionRules')),...
                         getString(message('stats:dfstrings:button_Ok')), ...
                         getString(message('stats:dfstrings:button_Cancel')), ...
                         getString(message('stats:dfstrings:button_Ok')));
        if ~strcmp(button, getString(message('stats:dfstrings:button_Ok')))
            tf = false;
        end
    end 
end

if tf
    if ~isempty(associatedFits)
        fm = stats.internal.dfit.FitsManager.getInstance;
        fm.deleteFits(cell2mat(associatedFits));
    end
end
end