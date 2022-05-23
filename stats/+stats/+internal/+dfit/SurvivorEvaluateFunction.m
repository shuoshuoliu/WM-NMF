classdef SurvivorEvaluateFunction < stats.internal.dfit.EvaluateFunction
    % SurvivorEvaluateFunction   Function that computes a survivor function
    % on the fits
    
    %   Copyright 2019 The MathWorks, Inc.
    
    methods
        function[errmsg, evaluateTable] = evaluateFits(this, selectedFits, vectorValue, wantBounds, confLevel, plotFun) %#ok<*INUSL>
            [errmsg, x, values] = stats.internal.dfit.evaluatefits(selectedFits, vectorValue, 'survivor', wantBounds, confLevel, plotFun);
            if ~isempty(errmsg)
                evaluateTable = iEmptyTable();
            else
                evaluateTable = iCreateEvaluatedTable(selectedFits, x, values, wantBounds);
            end
        end
        
        function label = getAppropriateVectorTextFieldLabel(this) %#ok<*MANU>
            label = iGetMessageString('stats:dfittool:label_atX');
        end
        
        function tf = canHaveConfidenceBoundOptionsEnabled(this)
            tf = true;
        end
    end
end

% helper functions
function tbl = iCreateEvaluatedTable(fitNames, x, values, showBounds)
tbl = array2table([x, values]);
if numel(fitNames) == 1
    tbl.Properties.VariableNames = iCreateColumnVariableNamesForOneFit(showBounds);
else
    tbl.Properties.VariableNames = iCreateColumnvariableNamesForMultipleFits(fitNames, showBounds);
end
end

function varNames = iCreateColumnVariableNamesForOneFit(showBounds)
xLabel = iGetMessageString('stats:dfittool:evalHeading_X');
SXLabel = iGetMessageString('stats:dfittool:evalHeading_survivor');
LBLabel = 'LB';
UBLabel = 'UB';
if ~showBounds
    varNames = {xLabel, SXLabel};
else
    varNames = {xLabel, SXLabel, LBLabel, UBLabel};
end
end

function varNames = iCreateColumnvariableNamesForMultipleFits(fitNames, showBounds)
xLabel = iGetMessageString('stats:dfittool:evalHeading_X');
SXLabel = iGetMessageString('stats:dfittool:evalHeading_survivor');
LBLabel = 'LB';
UBLabel = 'UB';
if ~showBounds
    varNames = {strcat(fitNames,{': '}, SXLabel)};
else
    varNames = cellfun(@(x) [strcat(x,{': '},SXLabel), strcat(x,{': '},LBLabel), strcat(x,{': '},UBLabel)],...
        fitNames,'Uniform',false);
end
varNames = [xLabel, varNames{:}];
end

function tbl = iEmptyTable()
tbl = table();
end

function outputString = iGetMessageString(inputString)
outputString = getString(message(inputString));
end