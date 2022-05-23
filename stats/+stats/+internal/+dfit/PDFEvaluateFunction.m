classdef PDFEvaluateFunction < stats.internal.dfit.EvaluateFunction
    % PDFEValuateFunction   Function that computes probability density
    % function (PDF) on the fits
    
    %   Copyright 2019 The MathWorks, Inc.
    
    methods
        function[errmsg, evaluateTable] = evaluateFits(this, selectedFits, vectorValue, wantBounds, confLevel, plotFun) %#ok<*INUSL>
            [errmsg, x, values] = stats.internal.dfit.evaluatefits(selectedFits, vectorValue, 'pdf', false, confLevel, plotFun);
            if ~isempty(errmsg)
                evaluateTable = iEmptyTable();
            else
                evaluateTable = iCreateEvaluatedTable(selectedFits, x, values);
            end
        end
        
        function label = getAppropriateVectorTextFieldLabel(this) %#ok<*MANU>
            label = iGetMessageString('stats:dfittool:label_atX');
        end
        
        function tf = canHaveConfidenceBoundOptionsEnabled(this)
            tf = false;
        end
    end
end

% helper functions
function tbl = iCreateEvaluatedTable(fitNames, x, values)
tbl = array2table([x, values]);
xLabel = iGetMessageString('stats:dfittool:evalHeading_X');
fXLabel = iGetMessageString('stats:dfittool:evalHeading_pdf');
if numel(fitNames) == 1
    tbl.Properties.VariableNames = {xLabel, fXLabel};
else
    tbl.Properties.VariableNames = [xLabel, strcat(fitNames,{': '},fXLabel)];
end
end

function tbl = iEmptyTable()
tbl = table();
end

function outputString = iGetMessageString(inputString)
outputString = getString(message(inputString));
end