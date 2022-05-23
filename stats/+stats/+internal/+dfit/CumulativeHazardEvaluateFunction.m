classdef CumulativeHazardEvaluateFunction < stats.internal.dfit.EvaluateFunction
    % CumulativeHazardEvaluateFunction   Function that computes a
    % cumulative hazard function on the fits
    
    %   Copyright 2019 The MathWorks, Inc.
    
    methods
        function[errmsg, evaluateTable] = evaluateFits(this, selectedFits, vectorValue, wantBounds, confLevel, plotFun) %#ok<*INUSL>
            [errmsg, x, values] = stats.internal.dfit.evaluatefits(selectedFits, vectorValue, 'cumhazard', wantBounds, confLevel, plotFun);
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
HXLabel = iGetMessageString('stats:dfittool:evalHeading_cumHaz');
LBLabel = 'LB';
UBLabel = 'UB';
if ~showBounds
    varNames = {xLabel, HXLabel};
else
    varNames = {xLabel, HXLabel, LBLabel, UBLabel};
end
end

function varNames = iCreateColumnvariableNamesForMultipleFits(fitNames, showBounds)
xLabel = iGetMessageString('stats:dfittool:evalHeading_X');
HXLabel = iGetMessageString('stats:dfittool:evalHeading_cumHaz');
LBLabel = 'LB';
UBLabel = 'UB';
if ~showBounds
    varNames = {strcat(fitNames,{': '}, HXLabel)};
else
    varNames = cellfun(@(x) [strcat(x,{': '},HXLabel), strcat(x,{': '},LBLabel), strcat(x,{': '},UBLabel)],...
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