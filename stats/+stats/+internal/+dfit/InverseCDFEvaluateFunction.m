classdef InverseCDFEvaluateFunction < stats.internal.dfit.EvaluateFunction
    % InverseCDFEvaluateFunction   Function that computes a quantile
    % (inverse CDF) function on the fits
    
    %   Copyright 2019 The MathWorks, Inc.
    
    methods
        function[errmsg, evaluateTable] = evaluateFits(this, selectedFits, vectorValue, wantBounds, confLevel, plotFun) %#ok<*INUSL>
            [errmsg, p, values] = stats.internal.dfit.evaluatefits(selectedFits, vectorValue, 'icdf', wantBounds, confLevel, plotFun);
            if ~isempty(errmsg)
                evaluateTable = iEmptyTable();
            else
                evaluateTable = iCreateEvaluatedTable(selectedFits, p, values, wantBounds);
            end
        end
        
        function label = getAppropriateVectorTextFieldLabel(this) %#ok<*MANU>
            label = iGetMessageString('stats:dfittool:label_atP');
        end
        
        function tf = canHaveConfidenceBoundOptionsEnabled(this)
            tf = true;
        end
    end
end

% helper functions
function tbl = iCreateEvaluatedTable(fitNames, p, values, showBounds)
tbl = array2table([p, values]);
if numel(fitNames) == 1
    tbl.Properties.VariableNames = iCreateColumnVariableNamesForOneFit(showBounds);
else
    tbl.Properties.VariableNames = iCreateColumnvariableNamesForMultipleFits(fitNames, showBounds);
end
end

function varNames = iCreateColumnVariableNamesForOneFit(showBounds)
pLabel = iGetMessageString('stats:dfittool:evalHeading_P');
FinvXLabel = iGetMessageString('stats:dfittool:evalHeading_invCDF');
LBLabel = 'LB';
UBLabel = 'UB';
if ~showBounds
    varNames = {pLabel, FinvXLabel};
else
    varNames = {pLabel, FinvXLabel, LBLabel, UBLabel};
end
end

function varNames = iCreateColumnvariableNamesForMultipleFits(fitNames, showBounds)
pLabel = iGetMessageString('stats:dfittool:evalHeading_P');
FinvXLabel = iGetMessageString('stats:dfittool:evalHeading_invCDF');
LBLabel = 'LB';
UBLabel = 'UB';
if ~showBounds
    varNames = {strcat(fitNames, {': '}, FinvXLabel)};
else
    varNames = cellfun(@(x) [strcat(x,{': '},FinvXLabel), strcat(x,{': '},LBLabel), strcat(x,{': '},UBLabel)],...
        fitNames,'Uniform',false);
end
varNames = [pLabel, varNames{:}];
end

function tbl = iEmptyTable()
tbl = table();
end

function outputString = iGetMessageString(inputString)
outputString = getString(message(inputString));
end