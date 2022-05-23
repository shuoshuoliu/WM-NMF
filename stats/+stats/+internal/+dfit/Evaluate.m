classdef Evaluate < handle
    % Evaluate   Dialog to evaluate fitted distribution at any data points
    % chosen
    
    %   Copyright 2019 The MathWorks, Inc.
    
    properties(Access = private)
        % DistributionFitting   (stats.internal.dfit.DistributionFitting)
        DistributionFitting
        
        % FitsManager   (stats.internal.dfit.FitsManager)
        FitsManager
        
        % EvaluatePanel   (stats.internal.dfit.EvaluatePanel)
        EvaluatePanel
        
        % SelectedFunction   (stats.internal.dfit.EvaluateFunction) Selected
        % proability function to evaluate on the fit(s)
        SelectedFunction
        
        % EvaluateResultsTable   (table) Store evaluated data
        EvaluateResultsTable
        
        % SessionClearingListener (listener) Listener to SessionClearingListener
        % event from the DistributionFitting object
        SessionClearingListener
        
        % Listeners   (cell of listener) Listeners to various events
        Listeners = {};
    end
    
    methods(Static)
        function this = getInstance()
            persistent Instance;
            if ~isa(Instance, 'stats.internal.dfit.Evaluate') || ~isvalid(Instance)
                Instance = stats.internal.dfit.Evaluate();
            end
            this = Instance;
        end
    end
    
    methods
        function showPanel(this)
            if ~isvalid(this.EvaluatePanel)
                this.createAndInitializeEvaluatePanel();
            end
            this.EvaluatePanel.showPanel();
        end
        
        function updatePlotCheckBoxValue(this, tf)
            this.EvaluatePanel.updatePlotCheckBox(tf);
        end
        
        function delete(this)
            delete(this.EvaluatePanel);
            % delete the evaluate plot (if there is one)
            stats.internal.dfit.evaluateplot(false);
        end
    end
    
    methods(Access = private)
        function this = Evaluate()
            this.DistributionFitting = stats.internal.dfit.DistributionFitting.getInstance();
            this.FitsManager = stats.internal.dfit.FitsManager.getInstance();
            this.createAndInitializeEvaluatePanel();
            
            this.SessionClearingListener = event.listener(this.DistributionFitting, 'SessionClearing', @(~,~) this.delete);
        end
        
        function createAndInitializeEvaluatePanel(this)
            fitList = iGetGoodFitNames();
            functionType = iGetDistributionFitterCurrentPlotType();
            vectorValue = iGetDistributionFitterHorizontalAxesLimit();
            this.SelectedFunction = iEvaluateFunctionFactory(functionType);
            this.EvaluateResultsTable = iEmptyTable();
            vectorLabel = this.SelectedFunction.getAppropriateVectorTextFieldLabel();
            isConfidenceBoundEnabled = this.SelectedFunction.canHaveConfidenceBoundOptionsEnabled();
            
            this.EvaluatePanel = stats.internal.dfit.EvaluatePanel(fitList, functionType, vectorLabel,...
                vectorValue, isConfidenceBoundEnabled, this.EvaluateResultsTable);
            
            this.addListener(this.FitsManager, 'FitAdded', @this.fitsAddedOrDeleted);
            this.addListener(this.FitsManager, 'FitsDeleted', @this.fitsAddedOrDeleted);
            this.addListener(this.FitsManager, 'FitChanged', @this.fitChanged);
            this.addListener(this.EvaluatePanel, 'ApplyButtonClicked', @this.applyButtonCallback);
            this.addListener(this.EvaluatePanel, 'ExportToWorkspaceButtonClicked', @this.exportToWorkspaceButtonCallback);
            this.addListener(this.EvaluatePanel, 'CloseButtonClicked', @this.closeButtonClickedCallback);
            this.addListener(this.EvaluatePanel, 'HelpButtonClicked', @this.helpButtonCallback);
            this.addListener(this.EvaluatePanel, 'FunctionDropDownValueChanged', @this.functionDropDownValueChangedCallback);
            this.addListener(this.EvaluatePanel, 'ParameterValuesChanged', @this.parameterValuesChangedCallback);
            this.addListener(this.EvaluatePanel, 'PlotFunctionCheckBoxValueChanged', @this.plotFunctionCheckBoxValueChangedCallback);
        end
       
        % Callbacks
        function fitsAddedOrDeleted(this, ~, ~)
            newFitList = iGetGoodFitNames();
            currentSelectedFits = this.EvaluatePanel.SelectedFits;
            
            newSelectedFits = intersect(currentSelectedFits, newFitList);
            if isempty(newSelectedFits) && ~isempty(newFitList)
                newSelectedFits = newFitList(1);
            end
            
            this.EvaluatePanel.updateFitListBox(newFitList, newSelectedFits);            
        end
        
        function fitChanged(this, ~, fittingEventData)
            newFitList = iGetGoodFitNames();
            currentSelectedFits = this.EvaluatePanel.SelectedFits;
            
            oldFitName = fittingEventData.OldName;
            newFitName = fittingEventData.NewName;
            idx = find(ismember(currentSelectedFits, oldFitName));
            if ~isempty(idx)
                currentSelectedFits{idx} = newFitName;
            end
            
            newSelectedFits = intersect(currentSelectedFits, newFitList);
            this.EvaluatePanel.updateFitListBox(newFitList, newSelectedFits);
        end
        
        function applyButtonCallback(this, ~, ~) %#ok<*INUSD>
            errmsg = iCheckPossibleErrorsInPanelBeforeEvaluating(this.EvaluatePanel);
            
            if ~isempty(errmsg)
                this.EvaluatePanel.showMessageDialog(iGetMessageString(errmsg));
                this.EvaluateResultsTable = iEmptyTable();
            else
                selectedFits = this.EvaluatePanel.SelectedFits;
                vectorValue = this.EvaluatePanel.VectorValue;
                showBounds = this.EvaluatePanel.IsComputeConfidenceBounds;
                levelValue = num2str(this.EvaluatePanel.LevelValue);
                plotFun = this.EvaluatePanel.IsPlotFunction;
                
                [errmsg, this.EvaluateResultsTable] = this.SelectedFunction.evaluateFits(...
                    selectedFits, vectorValue, showBounds, levelValue, plotFun);
                
                if ~isempty(errmsg)
                    this.EvaluatePanel.showMessageDialog(errmsg);
                end
            end
            
            this.EvaluatePanel.updateEvaluateResultsTable(this.EvaluateResultsTable);
        end
        
        function exportToWorkspaceButtonCallback(this, ~, ~)
            results = dfgetset('evaluateResults');
            export2wsdlg({getString(message('stats:dfstrings:cellstr_SaveToVariable'))}, ...
                {'evaluateresults'}, {results});
        end
        
        function helpButtonCallback(this, ~, ~)
             stats.internal.dfit.helpviewer('evaluate_fits', 'evaluate');
        end
        
        function functionDropDownValueChangedCallback(this, ~, ~)
            this.parameterValuesChangedCallback();
            
            this.SelectedFunction = iEvaluateFunctionFactory(this.EvaluatePanel.SelectedFunction);

            label = this.SelectedFunction.getAppropriateVectorTextFieldLabel();
            isConfidenceBoundEnabled = this.SelectedFunction.canHaveConfidenceBoundOptionsEnabled();
            
            this.EvaluatePanel.updateVectorTextFieldLabel(label);
            this.EvaluatePanel.updateConfidenceBoundComponentStates(isConfidenceBoundEnabled);
        end
        
        function parameterValuesChangedCallback(this, ~, ~)
            this.EvaluateResultsTable = iEmptyTable();
            this.EvaluatePanel.updateEvaluateResultsTable(this.EvaluateResultsTable);
            stats.internal.dfit.evaluateplot(false);
        end
        
        function plotFunctionCheckBoxValueChangedCallback(this, ~, ~)
            if ~isempty(this.EvaluateResultsTable)
                stats.internal.dfit.evaluateplot(this.EvaluatePanel.IsPlotFunction);
            end
        end
        
        function closeButtonClickedCallback(this, ~, ~)
            stats.internal.dfit.evaluateplot(false);
            this.deleteListeners();
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
function obj = iEvaluateFunctionFactory(functionName)
switch functionName
    case 'pdf'
        obj = stats.internal.dfit.PDFEvaluateFunction();
    case 'cdf'
        obj = stats.internal.dfit.CDFEvaluateFunction();
    case 'icdf'
        obj = stats.internal.dfit.InverseCDFEvaluateFunction();
    case 'survivor'
        obj = stats.internal.dfit.SurvivorEvaluateFunction();
    case 'cumhazard'
        obj = stats.internal.dfit.CumulativeHazardEvaluateFunction();
    case 'hazrate'
        obj = stats.internal.dfit.HazardRateEvaluateFunction();
    otherwise
        % shouldn't reach here, but if it does set default as PDF
        obj = stats.internal.dfit.PDFEvaluateFunction();
end
end

function tbl = iEmptyTable()
tbl = table([]);
end

function errmsg = iCheckPossibleErrorsInPanelBeforeEvaluating(evaluatePanel)
fitList = evaluatePanel.FitList;
errmsg = '';
if isempty(fitList)
    errmsg = 'stats:dfittool:error_no_fit';
elseif isempty(evaluatePanel.SelectedFits)
    errmsg = 'stats:dfittool:error_no_fit_selected';
elseif isempty(evaluatePanel.VectorValue)
    errmsg = 'stats:dfittool:error_no_expression';
end
end

function outputStr = iGetMessageString(inputStr)
outputStr = getString(message(inputStr));
end

function xstr = iGetDistributionFitterHorizontalAxesLimit()
xlims = dfgetset('xminmax');

if isempty(xlims)
    xstr = '';
else
    % Start out by choosing a rounding that will give the smaller (in
    % magnitude) of min(x) and max(x) a single sig digit, but at most two sig
    % digits in the larger.  The smaller may round to zero.  If the min and max
    % round to the same thing, use more digits until we get rounded numbers
    % that differ.
    xmag = max(min(abs(xlims)), max(abs(xlims))/10);
    rounder = 10^floor(log10(xmag));
    while true
        xmin = round(xlims(1)./rounder) * rounder;
        xmax = round(xlims(2)./rounder) * rounder;
        if xmin < xmax, break; end
        rounder = rounder/10;
    end
    
    % Create 10 steps, where the step size will have three more significant
    % digits than the endpoints.
    stepRounder = rounder ./ 1000;
    step = floor((xmax-xmin)./(10*stepRounder)) * stepRounder;
    
    % Figure out how many digits we need display in order to distinguish the
    % endpoints.  That's the number of digits to the left of the decimal, plus
    % however many we've rounded to on the right.  Allow for at least four so
    % that we'll get, e.g., "4000", and not "4e+03".
    xminDigits = max(round(log10(max(abs(xmin),1))-log10(rounder)),4);
    xmaxDigits = max(round(log10(max(abs(xmax),1))-log10(rounder)),4);
    xstr = sprintf('%0.*g:%g:%0.*g', xminDigits, xmin, step, xmaxDigits, xmax);
end
end

function functionName = iGetDistributionFitterCurrentPlotType()
functionName = dfgetset('ftype');
if strcmp(functionName,'probplot')
    functionName = 'cdf'; % can't evaluate a prob plot, use CDF instead
% else {'pdf' 'cdf' 'survivor' 'icdf' 'cumhazard' 'hazrate'}
    % otherwise use the current setting
end

if isempty(functionName) % shouldn't reach here, but set it as PDF for default
    functionName = 'pdf';
end
end

function databaseObj = iFitDatabase()
databaseObj = stats.internal.dfit.getfitdb;
end

function names = iGetGoodFitNames()
goodFits = find(iFitDatabase(), 'isgood', 1);
names = arrayfun(@(x) x.name, goodFits, 'UniformOutput', false);
end
