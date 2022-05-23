classdef MainPanel < handle
    % MainPanel   Main dialog for Distribution Fitter App
    
    %   Copyright 2019-2020 The MathWorks, Inc.
    
    properties
        % DistributionFitterUIFigure   (uifigure) The main window of the
        % distribution fitter app
        DistributionFitterUIFigure
    end
    
    properties(Access = private)
        % Grid layout
        MainGrid           
        DisplayTypeAndDistributionGrid
        ButtonsGrid                   
        MainAxesGrid
        
        % Menu Components
        FileMenu
        ImportDataMenu     
        ClearSessionMenu               
        LoadSessionMenu                
        SaveSessionMenu             
        GenerateCodeMenu             
        DefineCustomDistributionsMenu  
        ImportCustomDistributionsMenu 
        PrinttoFigureMenu    
        CloseDistributionFitterMenu    
        ViewMenu                       
        LegendMenu           
        GridMenu                       
        ConfidenceLevelMenu            
        Confidence90Menu              
        Confidence95Menu               
        Confidence99Menu               
        ConfidenceOtherMenu           
        ClearPlotMenu                 
        ToolsMenu                   
        AxesLimitControlMenu          
        DefaultAxesLimitMenu        
        SetDefaultBinRulesMenu      
        HelpMenu          
        StatisticsandMachineLearningToolboxHelpMenu 
        DistributionFitterHelpMenu 
        ExamplesMenu               
        
        % Display and distribution drop down components
        DisplayTypeLabel               
        DisplayTypeDropDown         
        DistributionLabel              
        DistributionDropDown        
        
        % Buttons to other panels
        DataButton                    
        NewFitsButton                 
        ManageFitsButton              
        EvaluateButton               
        ExcludeButton                
        
        % Main axes
        MainAxes                   
        
        % Axes limit control components
        YUpperLimitGrid  
        YLowerLimitGrid    
        XUpperLimitGrid    
        XLowerLimitGrid      
        YUpperLimitLabel 
        YLowerLimitLabel  
        XLowerLimitLabel
        XUpperLimitLabel  
        YUpperLimitSpinner
        YLowerLimitSpinner    
        XLowerLimitSpinner 
        XUpperLimitSpinner 
    end
    
    methods
        function this = MainPanel()
            this.createGUIComponents();
        end

        function delete(this)            
            %delete method closes the figure without verification
            this.DistributionFitterUIFigure.CloseRequestFcn = [];
            
            % Clear current session
            stats.internal.dfit.session('clear');
            
            % Delete any dfittool-related figures
            hEvaluate = dfgetset('evaluateFigure');
            if ~isempty(hEvaluate) && ishghandle(hEvaluate)
                delete(hEvaluate);
            end
            delete(this.DistributionFitterUIFigure);
            stats.internal.dfit.delgraphexclude;
            
            distributionFitting = stats.internal.dfit.DistributionFitting.getInstance();
            delete(distributionFitting);
        end
    end
    
    methods(Access = private)
        function createGUIComponents(this)
            this.DistributionFitterUIFigure = uifigure('Visible', 'off', ...
                'color',get(0,'defaultuicontrolbackgroundcolor'),...
                'Tag', 'Distribution Fitter Figure', ...
                'Name', getString(message('stats:dfstrings:dlg_DistributionFittingTool')), ...
                'numbertitle','off',...
                'Position', iGetInitialFigurePosition(), ...
                'CloseRequestFcn', @(~,~) this.closeFigureCallback(), ...
                'DeleteFcn', @(~,~) this.delete(), ...
                'WindowButtonDownFcn', @stats.internal.dfit.tips);
            
            dfgetset('dffig', this.DistributionFitterUIFigure);
            
            % Enable some zoom functionality
            h = zoom(this.DistributionFitterUIFigure);
            h.UseLegacyExplorationModes = 'on'; 

            this.createMenuSection();
            
            this.MainGrid = uigridlayout(this.DistributionFitterUIFigure);
            this.MainGrid.ColumnWidth = {'1x'};
            this.MainGrid.RowHeight = {35, 'fit', '1x'};
            this.MainGrid.Padding = [0 10 0 0];
            
            this.createDisplayTypeAndDistributionDropDownsSection();
            this.createButtonsSection();
            this.createMainAxesSection();
            
            feature('EnableUIComponentsInUIFigure', 1);
            stats.internal.dfit.docontext('create', this.DistributionFitterUIFigure);
            
            this.DistributionFitterUIFigure.Visible = 'on';
        end
        
        % GUI creation
        function createMenuSection(this)
            % Create File Menu
            this.FileMenu = uimenu(this.DistributionFitterUIFigure, ...
                'Text', getString(message('stats:dfstrings:label_File')), 'Tag', 'filemenu');
            this.ImportDataMenu = uimenu(this.FileMenu, ...
                'Text', getString(message('stats:dfstrings:label_ImportData')), ...
                'MenuSelectedFcn', @(~,~) this.dataButtonClickedCallback(), ...
                'Tag','dfitMenuImportData');
            this.ClearSessionMenu = uimenu(this.FileMenu, ...
                'Text', getString(message('stats:dfstrings:label_ClearSession')), ...
                'Callback', @(~,~) this.clearSessionMenuClickedCallback(), ...
                'Tag','dfitMenuImportClearSession', ...
                'Separator', 'on');
            this.LoadSessionMenu = uimenu(this.FileMenu, ...
                'Text', getString(message('stats:dfstrings:label_LoadSession')), ...
                'Callback', @(~,~) this.loadSessionMenuClickedCallback(), ...
                'Tag','dfitMenuLoadSession');
            this.SaveSessionMenu = uimenu(this.FileMenu, ...
                'Text', getString(message('stats:dfstrings:label_SaveSession')), ...
                'Callback', @(~,~) this.saveSessionMenuClickedCallback(), ...
                'Tag','dfitMenuSaveSession');
            this.GenerateCodeMenu = uimenu(this.FileMenu, ...
                'Text', getString(message('stats:dfstrings:label_GenerateCode')), ...
                'Callback', @(~,~) this.generateCodeMenuClickedCallback(), ...
                'Tag','dfitMenuGenCode');
            this.DefineCustomDistributionsMenu = uimenu(this.FileMenu, ...
                'Text', getString(message('stats:dfstrings:label_DefineCustomDistributions')),...
                'Callback',{@stats.internal.dfit.customdist, 'define'}', ...
                'Tag','dfitMenuDefineCustom', ...
                'Separator', 'on');
            this.ImportCustomDistributionsMenu = uimenu(this.FileMenu, ...
                'Text', getString(message('stats:dfstrings:label_ImportCustomDistributions')), ...
                'Callback',{@stats.internal.dfit.customdist, 'import'}, ...
                'Tag','dfitMenuImportCustom');
            this.PrinttoFigureMenu = uimenu(this.FileMenu,...
                'Text', getString(message('stats:dfstrings:label_PrintToFigure')) ,...
                'Callback', @(~,~) this.printToFigureMenuClickedCallback(), ...
                'Tag','dfitMenuPrint2Fig', ...
                'Separator', 'on');
            this.CloseDistributionFitterMenu = uimenu(this.FileMenu, ...
                'Text', getString(message('stats:dfstrings:label_CloseDistributionFitting')), ...
                'MenuSelectedFcn', @(~,~) this.closeFigureCallback(), ...
                'Tag','dfitMenuClose', ...
                'Separator', 'on');
            
            % Create View Menu
            this.ViewMenu = uimenu(this.DistributionFitterUIFigure, ...
                'Text', getString(message('stats:dfstrings:label_View')), ...
                'Tag', 'viewmenu');
            this.LegendMenu = uimenu(this.ViewMenu, ...
                'Text', getString(message('stats:dfstrings:label_Legend')), 'Separator','off',...
                'Callback', @(~,~) this.toggleLegendMenuCallback(), ...
                'Checked','on',...
                'Tag','showlegend');
            this.GridMenu = uimenu(this.ViewMenu, ...
                'Text', getString(message('stats:dfstrings:label_Grid')),...
                'Callback', @(~,~) this.toggleGridMenuCallback(), ...
                'Checked','off', ...
                'Tag','showgrid');
            this.ConfidenceLevelMenu = uimenu(this.ViewMenu, ...
                'Text', getString(message('stats:dfstrings:label_ConfidenceLevel')),...
                'Separator','on');
            this.Confidence90Menu = uimenu(this.ConfidenceLevelMenu, ...
                'Text', '9&0%',  ...
                'Callback', @(~,~) this.setConfidenceBoundMenuCallback(0.90), ...
                'Tag','conflev90');
            this.Confidence95Menu = uimenu(this.ConfidenceLevelMenu, ...
                'Text', '9&5%', ...
                'Callback', @(~,~) this.setConfidenceBoundMenuCallback(0.95), ...
                'Checked', 'on',...
                'Tag','conflev95');
            this.Confidence99Menu = uimenu(this.ConfidenceLevelMenu, ...
                'Text', '9&9%',  ...
                'Callback', @(~,~) this.setConfidenceBoundMenuCallback(0.99), ...
                'Tag','conflev99');
            this.ConfidenceOtherMenu = uimenu(this.ConfidenceLevelMenu, ...
                'Text', getString(message('stats:dfstrings:label_Other')),  ...
                'Callback', @(~,~) this.setConfidenceBoundMenuCallback([]), ...
                'Tag','conflevOther');
            this.ClearPlotMenu = uimenu(this.ViewMenu, ...
                'Text', getString(message('stats:dfstrings:label_ClearPlot')),...
                'Callback', @(~,~) this.clearPlotMenuClickedCallback(), ...
                'Tag','dfitMenuClearPlot');

            % Create Tools Menu
            this.ToolsMenu = uimenu(this.DistributionFitterUIFigure, ...
                'Text', getString(message('stats:dfstrings:label_Tools')), ...
                'Tag', 'toolsmenu');
            this.AxesLimitControlMenu = uimenu(this.ToolsMenu, ...
                'Text', getString(message('stats:dfstrings:label_AxesLimitControl')), ...
                'Callback', @(~,~) this.toggleAxesLimitControlCallback(),...
                'Tag', 'showaxlimctrl', ...
                'Checked', 'off', ...
                'Separator', 'on');
            this.DefaultAxesLimitMenu = uimenu(this.ToolsMenu, ...
                'Text', getString(message('stats:dfstrings:label_DefaultAxesLimits')), ...
                'Callback', @(~,~) this.defaultAxesMenuClickedCallback(), ...
                'Tag', 'dfitMenuDefaultAxes');
            this.SetDefaultBinRulesMenu = uimenu(this.ToolsMenu, ...
                'Text', getString(message('stats:dfstrings:label_SetDefaultBinRules')),...
                'Callback', @(~,~) this.defaultBinWidthRulesClickedCallback(), ...
                'Tag', 'setbinrules', ...
                'Separator','on');

            % Create Help Menu
            this.HelpMenu = uimenu(this.DistributionFitterUIFigure, ...
                'Text', getString(message('stats:dfstrings:label_Help')), ...
                'Tag', 'helpmenu');
            this.StatisticsandMachineLearningToolboxHelpMenu = uimenu(this.HelpMenu, ....
                'Text', getString(message('stats:dfstrings:label_StatisticsToolboxHelp')),...
                'Callback', 'doc stats',...
                'Tag','dfitMenuHelpTbx');
            this.DistributionFitterHelpMenu = uimenu(this.HelpMenu, ...
                'Text',  getString(message('stats:dfstrings:label_DistributionFittingToolHelp')), ...
                'Callback', @(varargin) stats.internal.dfit.helpviewer('distribution_fitting', 'dfittool'), ...
                'Tag', 'dfitMenuHelpDfit');
            this.ExamplesMenu = uimenu(this.HelpMenu, ...
                'Text', getString(message('stats:dfstrings:label_Demos')),...
                'Callback', 'demo toolbox stat', ...
                'Tag', 'dfitMenuDemos', ...
                'Separator', 'on');
            
            dfgetset('showlegend','on');
            dfgetset('showgrid','off');
            dfgetset('conflev',0.95);
            dfgetset('showaxlimctrl','off');
        end
        
        function createDisplayTypeAndDistributionDropDownsSection(this)
            p = uipanel(this.MainGrid);
            p.Layout.Row = 1;
            p.Layout.Column = 1;
            
            this.DisplayTypeAndDistributionGrid = uigridlayout(p);
            this.DisplayTypeAndDistributionGrid.ColumnWidth = {'fit', 'fit', '0.1x', 'fit', 'fit', '1x'};
            this.DisplayTypeAndDistributionGrid.RowHeight = {'fit'};
            this.DisplayTypeAndDistributionGrid.Padding = [10 10 10 6];
            
            this.DisplayTypeLabel = uilabel(this.DisplayTypeAndDistributionGrid, ...
                'Tag','displaytext',...
                'Text',getString(message('stats:dfstrings:uicontrol_DisplayType')), ...
                'HorizontalAlignment','left','FontWeight','bold');
            this.DisplayTypeLabel.Layout.Row = 1;
            this.DisplayTypeLabel.Layout.Column = 1;
            
            this.DisplayTypeDropDown = uidropdown(this.DisplayTypeAndDistributionGrid, ...
                'Tag','displaylist',...
                'Items',iGetDisplayTypeChoices(), ...
                'BackgroundColor', iWhite(), ...
                'ItemsData', [1 2 3 4 5 6], ...
                'ValueChangedFcn', @(~,~) this.displayTypeAndDistributionDropDownValueChangedCallback());
            this.DisplayTypeDropDown.Layout.Row = 1;
            this.DisplayTypeDropDown.Layout.Column = 2;
            setappdata(this.DisplayTypeDropDown, 'codenames',...
                {'pdf' 'cdf' 'icdf' 'probplot' 'survivor' 'cumhazard'});
            
            this.DistributionLabel = uilabel(this.DisplayTypeAndDistributionGrid, ...
                'Tag','typetext',...
                'Text',getString(message('stats:dfstrings:uicontrol_Distribution')), ...
                'HorizontalAlignment','left', ...
                'FontWeight','bold', ...
                'Enable', 'off');
            this.DistributionLabel.Layout.Row = 1;
            this.DistributionLabel.Layout.Column = 4;
            
            alldist = stats.internal.dfit.getdistributions();
            dnames = {alldist.name};
            dcodes = {alldist.code};
            islocscale = [alldist.islocscale];
            choices = dnames(islocscale);
            normalname = getString(message('stats:dfittool:NameNormal'));
            default = find(strcmpi(normalname,choices));
            if length(default)~=1
                default = 1;
            end
            itemsData = 1:length(choices);

            this.DistributionDropDown = uidropdown(this.DisplayTypeAndDistributionGrid, ...
                'Tag','typelist',...
                'Items',choices, ...
                'BackgroundColor', iWhite(), ...
                'ItemsData', itemsData, ...
                'Value', default, 'Enable','off', ...
                'ValueChangedFcn', @(~,~) this.displayTypeAndDistributionDropDownValueChangedCallback());
            this.DistributionDropDown.Layout.Row = 1;
            this.DistributionDropDown.Layout.Column = 5;   
            
            setappdata(this.DistributionDropDown, 'allfullnames', choices);
            setappdata(this.DistributionDropDown, 'allcodenames', dcodes(islocscale));
            setappdata(this.DistributionDropDown, 'alldistspec', alldist(islocscale));
            setappdata(this.DistributionDropDown, 'okcodenames', dcodes(islocscale));
            
            h0 = uibutton(this.DistributionFitterUIFigure, 'Visible', 'off');
            hVec = [h0, this.DisplayTypeLabel, this.DisplayTypeDropDown, this.DistributionLabel, this.DistributionDropDown];
            setappdata(this.DistributionFitterUIFigure, 'selectioncontrols', hVec);
        end
        
        function createButtonsSection(this)
            this.ButtonsGrid = uigridlayout(this.MainGrid);
            this.ButtonsGrid.ColumnWidth = {'1x', 'fit', 'fit', 'fit', 'fit', 'fit', '1x'};
            this.ButtonsGrid.RowHeight = {'1x'};
            this.ButtonsGrid.Padding = [0 0 0 0];
            this.ButtonsGrid.Layout.Row = 2;
            this.ButtonsGrid.Layout.Column = 1;
            
            iCreateDummyMinSizeButton(this.ButtonsGrid, 1, 2);
            iCreateDummyMinSizeButton(this.ButtonsGrid, 1, 3);
            iCreateDummyMinSizeButton(this.ButtonsGrid, 1, 4);
            iCreateDummyMinSizeButton(this.ButtonsGrid, 1, 5);
            iCreateDummyMinSizeButton(this.ButtonsGrid, 1, 6);
            
            this.DataButton = uibutton(this.ButtonsGrid, 'push',...
                'Text', getString(message('stats:dfstrings:button_Data')),...
                'Tooltip', getString(message('stats:dfstrings:button_ImportViewEtc')), ...
                'ButtonPushedFcn', @(~,~) this.dataButtonClickedCallback, ...
                'Tag', 'dfMainDataButton');
            this.DataButton.Layout.Row = 1;
            this.DataButton.Layout.Column = 2;
            
            this.NewFitsButton = uibutton(this.ButtonsGrid, 'push',...
                'Text', getString(message('stats:dfstrings:button_NewFit')),...
                'Tooltip', getString(message('stats:dfstrings:button_AddAFittedDistribution')), ...
                'ButtonPushedFcn', @(~,~) this.newFitButtonClickedCallback, ...
                'Tag', 'dfMainNewFitButton');
            this.NewFitsButton.Layout.Row = 1;
            this.NewFitsButton.Layout.Column = 3;
            
            this.ManageFitsButton = uibutton(this.ButtonsGrid, 'push',...
                'Text', getString(message('stats:dfstrings:button_ManageFits')),...
                'Tooltip', getString(message('stats:dfstrings:button_EditViewEtc')), ...
                'ButtonPushedFcn', @(~,~) this.manageFitButtonClickedCallback, ...
                'Tag', 'dfMainManageFitButton');
            this.ManageFitsButton.Layout.Row = 1;
            this.ManageFitsButton.Layout.Column = 4;
            
            this.EvaluateButton = uibutton(this.ButtonsGrid, 'push',...
                'Text', getString(message('stats:dfstrings:button_Evaluate')),...
                'Tooltip', getString(message('stats:dfstrings:button_EvaluateFitsEtc')), ...
                'ButtonPushedFcn', @(~,~) this.evaluateButtonClickedCallback, ...
                'Tag', 'dfMainEvaluateButton');
            this.EvaluateButton.Layout.Row = 1;
            this.EvaluateButton.Layout.Column = 5;
            
            this.ExcludeButton = uibutton(this.ButtonsGrid, 'push',...
                'Text', getString(message('stats:dfstrings:button_Exclude')),...
                'Tooltip', getString(message('stats:dfstrings:button_DefineRulesEtc')), ...
                'ButtonPushedFcn', @(~,~) this.excludeButtonClickedCallback, ...
                'Tag', 'dfMainExcludeButton');
            this.ExcludeButton.Layout.Row = 1;
            this.ExcludeButton.Layout.Column = 6;
            
            hButtons = [this.DataButton, this.NewFitsButton, this.ManageFitsButton, this.ExcludeButton];
            setappdata(this.DistributionFitterUIFigure,'buttoncontrols', hButtons);
        end
        
        function createMainAxesSection(this)
            this.MainAxesGrid = uigridlayout(this.MainGrid);
            this.MainAxesGrid.ColumnWidth = {'fit', 'fit', '1x', 'fit', 50};
            this.MainAxesGrid.RowHeight = {'fit', '1x', 'fit', 'fit'};
            this.MainAxesGrid.ColumnSpacing = 0;
            this.MainAxesGrid.RowSpacing = 0;
            this.MainAxesGrid.Padding = [0 0 0 0];
            this.MainAxesGrid.Layout.Row = 3;
            this.MainAxesGrid.Layout.Column = 1;
            
            % Create a panel to hold Axes object
            p = uipanel(this.MainAxesGrid, 'BorderType', 'none','Tag','dfMainAxesPanel');
            p.Layout.Row = [1 3];
            p.Layout.Column = [2 4];
            
            % Create Axes
            this.MainAxes = axes(p, ...
                'Tag', 'dfMainAxes', ...
                'units', 'normalized', ...
                'FontSize', 14, ...
                'Box', 'on', ...
                'XLimMode','manual', ...
                'YLimMode','manual', ...
                'ZLimMode','manual',...
                'CLimMode','manual', ...
                'AlimMode','manual');
                        
            % Update curves to cover current x limits
            addlistener(this.MainAxes, 'XLim', 'PostSet', @(src, evt) iLocalUpdateCurves(this.MainAxes));
            addlistener(this.MainAxes, 'XLim', 'PostSet', @this.updateSpinnerValues);
            addlistener(this.MainAxes, 'YLim', 'PostSet', @this.updateSpinnerValues);

            this.createAxesLimitControl();            
        end
        
        function createAxesLimitControl(this)
            % Create grids
            this.YUpperLimitGrid = uigridlayout(this.MainAxesGrid);
            this.YUpperLimitGrid.ColumnWidth = {'1x'};
            this.YUpperLimitGrid.ColumnSpacing = 0;
            this.YUpperLimitGrid.RowSpacing = 1;
            this.YUpperLimitGrid.Padding = [10 1 1 1];
            this.YUpperLimitGrid.Layout.Row = 1;
            this.YUpperLimitGrid.Layout.Column = 1;
            
            this.YLowerLimitGrid = uigridlayout(this.MainAxesGrid);
            this.YLowerLimitGrid.ColumnWidth = {'1x'};
            this.YLowerLimitGrid.ColumnSpacing = 1;
            this.YLowerLimitGrid.RowSpacing = 1;
            this.YLowerLimitGrid.Padding = [10 1 1 1];
            this.YLowerLimitGrid.Layout.Row = 3;
            this.YLowerLimitGrid.Layout.Column = 1;
            
            this.XLowerLimitGrid = uigridlayout(this.MainAxesGrid);
            this.XLowerLimitGrid.ColumnWidth = {'1x'};
            this.XLowerLimitGrid.ColumnSpacing = 1;
            this.XLowerLimitGrid.RowSpacing = 1;
            this.XLowerLimitGrid.Padding = [10 1 1 1];
            this.XLowerLimitGrid.Layout.Row = 4;
            this.XLowerLimitGrid.Layout.Column = 2;
            
            this.XUpperLimitGrid = uigridlayout(this.MainAxesGrid);
            this.XUpperLimitGrid.ColumnWidth = {'1x'};
            this.XUpperLimitGrid.ColumnSpacing = 1;
            this.XUpperLimitGrid.RowSpacing = 1;
            this.XUpperLimitGrid.Padding = [10 1 1 1];
            this.XUpperLimitGrid.Layout.Row = 4;
            this.XUpperLimitGrid.Layout.Column = 4;
            
            % Create labels
            this.YUpperLimitLabel = uilabel(this.YUpperLimitGrid, ...
                'Text', getString(message('stats:dfstrings:xlate_YUpperLimit')), ...
                'Tag', 'YUpperLimitLabel');
            this.YUpperLimitLabel.Layout.Row = 1;
            this.YUpperLimitLabel.Layout.Column = 1;
            
            this.YLowerLimitLabel = uilabel(this.YLowerLimitGrid, ...
                'Text', getString(message('stats:dfstrings:xlate_YLowerLimit')), ...
                'Tag', 'YLowerLimitLabel');
            this.YLowerLimitLabel.Layout.Row = 1;
            this.YLowerLimitLabel.Layout.Column = 1;
            
            this.XLowerLimitLabel = uilabel(this.XLowerLimitGrid, ...
                'Text', getString(message('stats:dfstrings:xlate_XLowerLimit')), ...
                'Tag', 'XLowerLimitLabel');
            this.XLowerLimitLabel.Layout.Row = 2;
            this.XLowerLimitLabel.Layout.Column = 1;
            
            this.XUpperLimitLabel = uilabel(this.XUpperLimitGrid, ...
                'Text', getString(message('stats:dfstrings:xlate_XUpperLimit')), ...
                'Tag', 'XUpperLimitLabel');
            this.XUpperLimitLabel.Layout.Row = 2;
            this.XUpperLimitLabel.Layout.Column = 1;
            
            % Create spinners
            this.YLowerLimitSpinner = uispinner(this.YLowerLimitGrid, ...
                'HorizontalAlignment', 'center', ...
                'LowerLimitInclusive', 'off', ...
                'UpperLimitInclusive', 'off', ...
                'ValueChangingFcn', @this.spinnerValueChangingCallback, ...
                'Tag', 'YLowerLimitSpinner');
            this.YLowerLimitSpinner.Layout.Row = 2;
            this.YLowerLimitSpinner.Layout.Column = 1;
            
            this.YUpperLimitSpinner = uispinner(this.YUpperLimitGrid, ...
                'HorizontalAlignment', 'center', ...
                'LowerLimitInclusive', 'off', ...
                'UpperLimitInclusive', 'off', ...
                'ValueChangingFcn', @this.spinnerValueChangingCallback, ...
                'Tag', 'YUpperLimitSpinner');
            this.YUpperLimitSpinner.Layout.Row = 2;
            this.YUpperLimitSpinner.Layout.Column = 1;
            
            this.XLowerLimitSpinner = uispinner(this.XLowerLimitGrid, ...
                'HorizontalAlignment', 'center', ...
                'LowerLimitInclusive', 'off', ...
                'UpperLimitInclusive', 'off', ...
                'ValueChangingFcn', @this.spinnerValueChangingCallback, ...
                'Tag', 'XLowerLimitSpinner');
            this.XLowerLimitSpinner.Layout.Row = 1;
            this.XLowerLimitSpinner.Layout.Column = 1;
            
            this.XUpperLimitSpinner = uispinner(this.XUpperLimitGrid, ...
                'HorizontalAlignment', 'center', ...
                'LowerLimitInclusive', 'off', ...
                'UpperLimitInclusive', 'off', ...
                'ValueChangingFcn', @this.spinnerValueChangingCallback, ...
                'Tag', 'XUpperLimitSpinner');
            this.XUpperLimitSpinner.Layout.Row = 1;
            this.XUpperLimitSpinner.Layout.Column = 1;
            
            this.showAxesLimitSpinners(false);
        end
        
        % Menu callbacks
        function clearSessionMenuClickedCallback(this, ~, ~)
            delete(findall(gcbf,'Tag','dfstarthint'));
            isOK = stats.internal.dfit.asksavesession(this.DistributionFitterUIFigure);
            if isOK
                stats.internal.dfit.session('clear');
            end
        end
        
        function saveSessionMenuClickedCallback(~, ~, ~)
            stats.internal.dfit.session('save');
        end
        
        function loadSessionMenuClickedCallback(this, ~, ~)
            delete(findall(gcbf,'Tag','dfstarthint'));
            isOK = stats.internal.dfit.asksavesession(this.DistributionFitterUIFigure);
            if isOK
                stats.internal.dfit.session('load');
            end
        end
        
        function generateCodeMenuClickedCallback(~, ~, ~)
            stats.internal.dfit.fig2m();
        end
        
        function defaultBinWidthRulesClickedCallback(~, ~, ~)
            defaultBinRulesPanel = stats.internal.dfit.DefaultBinRulesPanel.getInstance();
            defaultBinRulesPanel.showPanel();
        end
        
        function printToFigureMenuClickedCallback(this, ~, ~)
            stats.internal.dfit.dupfigure(this.DistributionFitterUIFigure);
        end
        
        function toggleLegendMenuCallback(this, ~, ~)
            stats.internal.dfit.togglelegend(this.DistributionFitterUIFigure);
        end
        
        function toggleGridMenuCallback(this, ~, ~)
            stats.internal.dfit.togglegrid(this.DistributionFitterUIFigure);
        end
        
        function setConfidenceBoundMenuCallback(this, value, ~, ~)
            if (stats.internal.dfit.setconflev(this.DistributionFitterUIFigure, value))
                dfgetset('dirty',true);   % session has changed since last save
            end
        end
        
        function clearPlotMenuClickedCallback(this, ~, ~)
            delete(findall(this.DistributionFitterUIFigure,'Tag','dfstarthint'));
            stats.internal.dfit.cbkclear();
        end
        
        function toggleAxesLimitControlCallback(this, ~, ~)
            % Get new state
            onoff = stats.internal.dfit.on2off(get(this.AxesLimitControlMenu, 'Checked'));
            dfgetset('showaxlimctrl', onoff);
            
            % Add or remove controls
            this.showAxesLimitSpinners(onoff);
            
            % Change menu state
            this.AxesLimitControlMenu.Checked = onoff;
        end
        
        function defaultAxesMenuClickedCallback(~, ~, ~)
            stats.internal.dfit.updatexlim([],true,true);
            stats.internal.dfit.updateylim(true);
        end
        
        function closeFigureCallback(this, ~, ~)
            isOK = stats.internal.dfit.asksavesession();
            if isOK
                delete(this);
            end
        end
        
        % Drop-down callbacks
        function displayTypeAndDistributionDropDownValueChangedCallback(this, ~, ~)
            % Get the requested function and distribution types
            drawnow('expose'); 
            stats.internal.dfit.setplottype(this.DistributionFitterUIFigure);
            
            % Update legend position
            stats.internal.dfit.updatelegend(this.DistributionFitterUIFigure, true);
            
            this.showAxesLimitSpinners('off');
            set(this.AxesLimitControlMenu, 'Checked', 'off');
        end
                
        % Axis limit control
        function showAxesLimitSpinners(this, onoff)
            this.YUpperLimitLabel.Visible = onoff;
            this.YLowerLimitLabel.Visible = onoff;
            this.XLowerLimitLabel.Visible = onoff;
            this.XUpperLimitLabel.Visible = onoff;
            this.YUpperLimitSpinner.Visible = onoff;
            this.YLowerLimitSpinner.Visible = onoff;
            this.XLowerLimitSpinner.Visible = onoff;
            this.XUpperLimitSpinner.Visible = onoff;
            
            % Resize MainAxes grid layout
            if isequal(onoff, 'on')
                this.MainAxesGrid.ColumnWidth = {'fit', 'fit', '1x', 'fit', 25};
                this.MainAxesGrid.RowHeight = {'fit', '1x', 'fit', 'fit'};
                this.updateSpinnerValues();
            else
                this.MainAxesGrid.ColumnWidth = {0, 'fit', '1x', 'fit',0};
                this.MainAxesGrid.RowHeight = {0, '1x', 'fit', 0};
            end
        end
        
        function updateSpinnerValues(this, ~, ~)
            this.updateSpinnerLimits();
            iLocalUpdateText(this.MainAxes, 'x', this.XLowerLimitSpinner, this.XUpperLimitSpinner);
            iLocalUpdateText(this.MainAxes, 'y', this.YLowerLimitSpinner, this.YUpperLimitSpinner);
            this.updateSpinnerStepSizes();
        end 
        
        function spinnerValueChangingCallback(this, ~, eventData)
            curVal = eventData.Value;
            if ~isnumeric(curVal)
                curVal = str2num(curVal);
            end
                        
            switch eventData.Source.Tag
                case 'XLowerLimitSpinner'
                    xLim = this.MainAxes.XLim;
                    xLim(1) = curVal;
                    stats.internal.dfit.updatexlim(xLim);
                case 'XUpperLimitSpinner'
                    xLim = this.MainAxes.XLim;
                    xLim(2) = curVal;
                    stats.internal.dfit.updatexlim(xLim);
                case 'YLowerLimitSpinner'
                    yLim = this.MainAxes.YLim;
                    yLim(1) = iMapYSpinnerValueToAxesValue(this.DistributionFitterUIFigure, this.MainAxes, curVal);
                    this.MainAxes.YLim = yLim;
                case 'YUpperLimitSpinner'
                    yLim = this.MainAxes.YLim;
                    yLim(2) = iMapYSpinnerValueToAxesValue(this.DistributionFitterUIFigure, this.MainAxes, curVal);
                    this.MainAxes.YLim = yLim;
            end
        end
        
        function updateSpinnerLimits(this)
            this.XLowerLimitSpinner.Limits = [-Inf, this.MainAxes.XLim(2)];
            this.XUpperLimitSpinner.Limits = [this.MainAxes.XLim(1), Inf];
            
            invcdffun = getappdata(this.MainAxes, 'InverseCdfFunction');
            if isempty(invcdffun)
                this.YLowerLimitSpinner.LowerLimitInclusive = 'off';
                this.YLowerLimitSpinner.UpperLimitInclusive = 'off';
                this.YUpperLimitSpinner.LowerLimitInclusive = 'off';
                this.YUpperLimitSpinner.UpperLimitInclusive = 'off';
                
                this.YLowerLimitSpinner.Limits = [-Inf, this.MainAxes.YLim(2)];
                this.YUpperLimitSpinner.Limits = [this.MainAxes.YLim(1), Inf];
            else
                % For probability plot:
                this.YLowerLimitSpinner.LowerLimitInclusive = 'on';
                this.YLowerLimitSpinner.UpperLimitInclusive = 'off';
                this.YUpperLimitSpinner.LowerLimitInclusive = 'off';
                this.YUpperLimitSpinner.UpperLimitInclusive = 'on';
                
                ySpinnerLimit = iMapYLimValuesToLowerAndUpperSpinnerLimits(this.MainAxes);
                
                this.YLowerLimitSpinner.Limits = [0, ySpinnerLimit(2)];
                this.YUpperLimitSpinner.Limits = [ySpinnerLimit(1), 1];
            end
        end
        
        function updateSpinnerStepSizes(this)
            this.XLowerLimitSpinner.Step = this.MainAxes.XTick(2) - this.MainAxes.XTick(1);
            this.XUpperLimitSpinner.Step = this.MainAxes.XTick(end) - this.MainAxes.XTick(end-1);
            
            if isequal(this.MainAxes.YTickMode, 'manual')
                this.YLowerLimitSpinner.Step = iComputeYSpinnerStepSize(this.MainAxes, 'lower', this.YLowerLimitSpinner.Value);
                this.YUpperLimitSpinner.Step = iComputeYSpinnerStepSize(this.MainAxes, 'upper', this.YUpperLimitSpinner.Value);
            else
                this.YLowerLimitSpinner.Step = this.MainAxes.YTick(2) - this.MainAxes.YTick(1);
                this.YUpperLimitSpinner.Step = this.MainAxes.YTick(end) - this.MainAxes.YTick(end-1);
            end
        end
        
        % Button callbacks
        function dataButtonClickedCallback(~, ~, ~)
            delete(findall(gcbf,'Tag','dfstarthint'));
            d = stats.internal.dfit.Data.getInstance;
            d.showPanel();
        end
        
        function newFitButtonClickedCallback(~, ~, ~)
            stats.internal.dfit.FitEditor();
        end
        
        function manageFitButtonClickedCallback(~, ~, ~)
            fitsManager = stats.internal.dfit.FitsManager.getInstance;
            fitsManager.showPanel();
        end
        
        function evaluateButtonClickedCallback(~, ~, ~)
            evaluate = stats.internal.dfit.Evaluate.getInstance;
            evaluate.showPanel();
        end
        
        function excludeButtonClickedCallback(~, ~, ~)
            exclude = stats.internal.dfit.Exclude.getInstance;
            exclude.showPanel();
        end
    end
end

% Helpers
function position = iGetInitialFigurePosition()
tempFigure=uifigure('visible','off','units','pixels');
dfp=get(tempFigure,'position');
dfop=get(tempFigure,'outerposition');
diffp = dfop - dfp;
xmargin = diffp(3);
ymargin = diffp(4);
close(tempFigure)
oldu = get(0,'units');
set(0,'units','pixels');
screenSize=get(0,'screensize');
screenWidth=screenSize(3); 
screenHeight=screenSize(4);
set(0,'units',oldu');

% Get the desired width and height
width=dfp(3)*1.2 + xmargin;
height=dfp(4)*1.3 + ymargin;
if width > screenWidth
  width = screenWidth-10-xmargin;
end
if height > screenHeight
  height = screenHeight-10-ymargin;
end

% Calculate the position on the screen
leftEdge=min((screenWidth/3)+10+xmargin/2, screenWidth-width-10-2*xmargin);
bottomEdge=(screenHeight-height)/2;
width = max(1,width);
height = max(1,height);

position = [leftEdge, bottomEdge, width, height];
end

function choices = iGetDisplayTypeChoices()
choices = {getString(message('stats:dfstrings:dropdown_Density')) ...
           getString(message('stats:dfstrings:dropdown_CDF')) ...
           getString(message('stats:dfstrings:dropdown_Quantile')) ...
           getString(message('stats:dfstrings:dropdown_ProbabilityPlot'))...
           getString(message('stats:dfstrings:dropdown_SurvivorFunction')) ...
           getString(message('stats:dfstrings:dropdown_CumulativeHazard'))};
end

function rgbArray = iWhite()
rgbArray = [1,1,1];
end

function iLocalUpdateCurves(ax)
xlim = get(ax,'xlim');
fitdb = stats.internal.dfit.getfitdb;
ft = down(fitdb);
while(~isempty(ft))                    % loop over all fits
    fitxlim = get(ft,'xlim');
    if isempty(fitxlim)
        % fit has no x limits yet, so update it
        doit = true;
        newlim = xlim;
    else
        xrange = abs(diff(xlim));
        frange = abs(diff(fitxlim));
        if frange>2*xrange || xrange>2*frange
            % drastic change in scale, re-compute limits
            newlim = xlim;
            doit = true;
        elseif fitxlim(1)>xlim(1) || fitxlim(2)<xlim(2)
            % minor change in scale, expand limits
            newlim = [min(xlim(1),fitxlim(1)), max(xlim(2),fitxlim(2))];
            doit = true;
        else
            % old scale is still okay, so avoid performance hit
            doit = false;
        end
    end
    if doit
        try
            updateplot(ft,newlim);
        catch
            % likely to be dealt with elsewhere
        end
    end
    ft = right(ft);
end
end

function iCreateDummyMinSizeButton(grid, row, column)
dummyText = repmat('X',1,8);
minButton = uibutton(grid, 'Text', dummyText);
minButton.Layout.Row = row;
minButton.Layout.Column = column;
minButton.Visible = 'off';
end

% Helper related to axes limit control
function iLocalUpdateText(ax, xory, lowerSpinner, upperSpinner)
% Get information about the axis in question
tickloc = [];
if xory=='x'
    lim = get(ax,'XLim');
    if isequal(get(ax,'XTickLabelMode'),'manual')
        tickloc = get(ax,'XTick');
        ticklabel = get(ax,'XTickLabel');
    end
    cdffun = '';
    logscale = isequal(get(ax,'XScale'),'log');
else
    lim = get(ax,'YLim');
    if isequal(get(ax,'YTickLabelMode'),'manual')
        tickloc = get(ax,'YTick');
        ticklabel = get(ax,'YTickLabel');
    end
    cdffun = getappdata(ax,'CdfFunction');
    if ~isempty(cdffun)
        params = getappdata(ax,'DistributionParameters');
    end
    logscale = isequal(get(ax,'XScale'),'log');
end

% Loop over the low and high values, and update text box with limits
hvec = [lowerSpinner upperSpinner];
small = max(abs(lim))*sqrt(eps);
for jLim=1:2
    newval = lim(jLim);
    if ~isempty(tickloc) && any(abs(tickloc-newval) < small)
        % If this tick has a label already made, use it
        j = find(abs(tickloc-newval) < small);
        j = j(1);
        newval = deblank(ticklabel(j,:));
    else
        % Otherwise make a new label, calling cdf function if necessary
        if ~isempty(cdffun)
            newval = feval(cdffun,newval,params{:}); %#ok<FVAL>
        end
        if abs(newval)<small && ~logscale
            % Just to get 0 instead of ~eps when appropriate
            newval = 0;
        end
    end
    if ischar(newval)
        newval = str2num(newval); 
    end
    hvec(jLim).Value = newval;
end
end

function axesValue = iMapYSpinnerValueToAxesValue(fHandle, axHandle, spinnerValue)
axesValue = spinnerValue;
invcdffun = getappdata(axHandle,'InverseCdfFunction');
if ~isempty(invcdffun)
    params = getappdata(axHandle,'DistributionParameters');
    axesValue = feval(invcdffun, spinnerValue, params{:});
end
if ~isfinite(axesValue)
    uialert(fHandle, ...
        getString(message('stats:dfstrings:dlg_BadAxisLimitValue')),...`
        getString(message('stats:dfstrings:dlg_DistributionFittingTool')), ...
        'Icon', 'warning');
    
    axesValue = spinnerValue;
end
end

function spinnerLimits = iMapYLimValuesToLowerAndUpperSpinnerLimits(axHandle)
lim = get(axHandle,'YLim');
if isequal(get(axHandle,'YTickLabelMode'),'manual')
    tickloc = get(axHandle,'YTick');
    ticklabel = get(axHandle,'YTickLabel');
else
    tickloc = [];
    ticklabel = [];
end
cdffun = getappdata(axHandle,'CdfFunction');
if ~isempty(cdffun)
    params = getappdata(axHandle,'DistributionParameters');
end
logscale = isequal(get(axHandle,'XScale'),'log');

small = max(abs(lim))*sqrt(eps);
for jLim=1:2
    newval = lim(jLim);
    if ~isempty(tickloc) && any(abs(tickloc-newval) < small)
        % If this tick has a label already made, use it
        j = find(abs(tickloc-newval) < small);
        j = j(1);
        newval = deblank(ticklabel(j,:));
    else
        % Otherwise make a new label, calling cdf function if necessary
        if ~isempty(cdffun)
            newval = feval(cdffun,newval,params{:});
        end
        if abs(newval)<small && ~logscale
            % Just to get 0 instead of ~eps when appropriate
            newval = 0;
        end
    end
    if ischar(newval)
        newval = str2num(newval); 
    end
    spinnerLimits(jLim) = newval; %#ok<AGROW>
end
end

function stepSize = iComputeYSpinnerStepSize(axHandle, lowerOrUpperSpinner, spinnerValue)
tickLabelValues = str2num(axHandle.YTickLabel); %#ok<*ST2NM>
[~, idx] = min(abs(spinnerValue - tickLabelValues));

stepSize = 1;
if isequal(lowerOrUpperSpinner, 'lower')
    if idx ~= 1
        stepSize = abs(spinnerValue - tickLabelValues(idx-1));
    end
else
    if idx ~= numel(tickLabelValues)
        stepSize = abs(tickLabelValues(idx+1) - spinnerValue);
    end
end
if ~isfinite(stepSize)
    stepSize = 1;
end
end

