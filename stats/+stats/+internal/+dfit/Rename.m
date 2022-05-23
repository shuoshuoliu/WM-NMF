classdef Rename < handle
    % Rename   Dialog to rename objects in Distribution Fitter database
    
    %   Copyright 2019-2020 The MathWorks, Inc.
    
    properties(Access = private)
        % Panel   (uifigure) The main window of the Rename dialog
        Panel
        
        % UI Components
        NameTextField
        OKButton
        CancelButton
        
        % DatabaseObj   (UDD object) database object to access
        DatabaseObj
        
        % Initialname   (string) initial string to show in the 
        % NameTextField when the dialog is made visible 
        InitialName = '';
        
        % SessionClearingListener   (listener) Listener to SessionClearing event in
        % DistributionFitting object
        SessionClearingListener
    end
    
    events
        RenameSuccessful
    end
    
    methods
        function this = Rename(databaseObj)
            this.DatabaseObj = databaseObj;
            
            this.SessionClearingListener = event.listener(stats.internal.dfit.DistributionFitting.getInstance(), 'SessionClearing', @(~,~) this.delete);
        end
        
        function showPanel(this, initialName)
            if isempty(this.Panel) || ~isvalid(this.Panel)
                this.createGUIComponents(iGetWindowTitle(this.DatabaseObj));
            end
            this.InitialName = initialName;
            this.NameTextField.Value = initialName;
            figure(this.Panel);
        end
        
        function hidePanel(this)
            delete(this.Panel);
        end
        
        function delete(this)
            delete(this.Panel);
        end
    end
    
    methods(Access = private) 
        function createGUIComponents(this, title)
            this.Panel = uifigure(...
                'Name', title, ...
                'Position', iGetDefaultPanelPosition(), ...
                'Visible', 'off', ...
                'CloseRequestFcn', @(~,~) this.hidePanel(), ...
                'WindowStyle', 'modal', ...
                'Tag', 'dfRenameUIFigure');
            
            grid = uigridlayout(this.Panel, [2, 3]);
            grid.ColumnWidth = {'1x',  'fit',  'fit'};
            grid.RowHeight = {'fit', 'fit'};
            
            this.createNameTextField(grid);
            this.createOKButton(grid);
            this.createCancelButton(grid);
        end
        
        function renameObject(this, oldName, newName)
            if strcmp(newName, "")
                messageStr = getString(message('stats:dfittool:error_noName', iGetObjectType(this.DatabaseObj)));
                titleStr = getString(message('stats:dfittool:title_noObjectTypeName', iGetTitleObjectType(this.DatabaseObj)));
                this.showMessageDialog(messageStr, titleStr);
                return
            end
            
            trimmedNewName = strtrim(newName);
            if strcmp(trimmedNewName, "")
                messageStr = iGetMessageString('stats:dfittool:error_noOnlySpacesAllowed');
                titleStr = getString(message('stats:dfittool:title_invalidObjectTypeName', iGetTitleObjectType(this.DatabaseObj)));
                this.showMessageDialog(messageStr, titleStr);
                return
            end
            
            if strcmp(trimmedNewName, oldName)
                this.hidePanel();
                return
            end
            
            existingNames = stats.internal.dfit.getDataBaseNames(this.DatabaseObj);
            
            if ismember(trimmedNewName, existingNames)
                messageStr = getString(message('stats:dfittool:error_duplicateName', iGetMixedObjectType(this.DatabaseObj), trimmedNewName));
                titleStr = getString(message('stats:dfittool:title_duplicateName', iGetTitleObjectType(this.DatabaseObj)));
                this.showMessageDialog(messageStr, titleStr);
                return
            end
            
            if ~strcmp(oldName, '')
                objToBeRenamed = find(this.DatabaseObj, 'name', oldName);
                objToBeRenamed.name = trimmedNewName;
                renameEventData = iCreateRenameEventData(this.DatabaseObj, oldName, newName);
                this.notify('RenameSuccessful', renameEventData);
                this.hidePanel();
            end
        end
        
        function showMessageDialog(this, messageStr, titleStr)
            figure(this.Panel);
            stats.internal.dfit.errordlg(messageStr, titleStr);
        end
        
        % UI components
        function createNameTextField(this, parentGrid)
            this.NameTextField = uieditfield(parentGrid,'Tag','dfRenameEditField'); %, ...
            %'ValueChangedFcn', @(~,~) this.oKButtonClickedCallback());
            this.NameTextField.Layout.Row = 1;
            this.NameTextField.Layout.Column = [1 3];
        end
        
        function createOKButton(this, parentGrid)
            iCreateDummyMinSizeButton(parentGrid, 2, 2);
            
            this.OKButton = uibutton(parentGrid, ...
                'Text', iGetMessageString('stats:dfittool:button_OK'), ...
                'ButtonPushedFcn', @(~, ~)this.oKButtonClickedCallback(),...
                'Tag','dfRenameOKButton');
            this.OKButton.Layout.Row = 2;
            this.OKButton.Layout.Column = 2;
        end
        
        function createCancelButton(this, parentGrid)
            iCreateDummyMinSizeButton(parentGrid, 2, 3);

            this.CancelButton = uibutton(parentGrid, ...
                'Text', iGetMessageString('stats:dfittool:button_cancel'), ...
                'ButtonPushedFcn', @(~, ~)this.cancelButtonClickedCallback(),...
                'Tag','dfRenameCancelButton');
            this.CancelButton.Layout.Row = 2;
            this.CancelButton.Layout.Column = 3;
        end
        
        % Callbacks
        function oKButtonClickedCallback(this, ~, ~)
            oldName = this.InitialName;
            newName = this.NameTextField.Value;
            this.renameObject(oldName, newName);
        end
        
        function cancelButtonClickedCallback(this, ~, ~)
            this.hidePanel();
        end
    end
end

% helper function
function position = iGetDefaultPanelPosition()
width = 250;
height = 75;

callingPanelPosition = get(gcbf, 'Position');
x = callingPanelPosition(1) + (callingPanelPosition(3) - width)/2;
y = callingPanelPosition(2) + (callingPanelPosition(4) - height)/2;

position = [x, y, width, height];
end

function iCreateDummyMinSizeButton(grid, row, column)
dummyText = repmat('X',1,8);
minButton = uibutton(grid, 'Text', dummyText);
minButton.Layout.Row = row;
minButton.Layout.Column = column;
minButton.Visible = 'off';
end

function outputStr = iGetMessageString(inputStr)
outputStr = getString(message(inputStr));
end

function titleStr = iGetWindowTitle(databaseObj)
switch class(databaseObj)
    case 'stats.dsdb'
        titleStr = iGetMessageString('stats:dfittool:title_renameDataSet');
    case 'stats.outlierdb'
        titleStr = iGetMessageString('stats:dfittool:title_renameOutlier');
    otherwise
        titleStr = '';
end
end

function objType = iGetObjectType(databaseObj)
switch class(databaseObj)
    case 'stats.dsdb'
        objType = iGetMessageString('stats:dfittool:label_dataSetLowerCase');
    case 'stats.outlierdb'
        objType = iGetMessageString('stats:dfittool:label_exclusionRuleLowerCase');
    otherwise
        objType = '';
end
end

function objType = iGetMixedObjectType(databaseObj)
switch class(databaseObj)
    case 'stats.dsdb'
        objType = iGetMessageString('stats:dfittool:label_dataSetMixedCase');
    case 'stats.outlierdb'
        objType = iGetMessageString('stats:dfittool:label_exclusionRuleMixedCase');
    otherwise
        objType = '';
end
end

function objType = iGetTitleObjectType(databaseObj)
switch class(databaseObj)
    case 'stats.dsdb'
        objType = iGetMessageString('stats:dfittool:label_dataSetUpperCase');
    case 'stats.outlierdb'
        objType = iGetMessageString('stats:dfittool:label_exclusionRuleUpperCase');
    otherwise
        objType = '';
end
end

function eventData = iCreateRenameEventData(dataBaseObj, oldName, newName)
eventData = stats.internal.dfit.RenameEventData(class(dataBaseObj), oldName, newName);
end


    