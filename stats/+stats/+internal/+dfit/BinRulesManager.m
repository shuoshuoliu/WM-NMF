classdef BinRulesManager < handle
    % BinRulesManager    manages default and data set bin rules
    
    %   Copyright 2019-2020 The MathWorks, Inc.
    
    events
        ApplyBinRuleToAllDataSets
        SaveBinRuleAsDefault
    end
    
    methods(Static)
        function this = getInstance()
            persistent Instance;
            if ~isa(Instance, 'stats.internal.dfit.BinRulesManager') || ~isvalid(Instance)
                Instance = stats.internal.dfit.BinRulesManager();
            end
            this = Instance;
        end
    end
    
    methods(Access = private)
        function this = BinRulesManager()
        end
    end
    
    methods
        function applyBinRuleToDataSet(~, dataSetObj)
            if ~isempty(dataSetObj.viewData)
                dataSetObj.viewData.updateDistributionPreviewAxes();
            end
        end
        
        function applyBinRuleToAllDataSets(this)
            this.notify('ApplyBinRuleToAllDataSets');
        end
        
        function saveBinRuleAsDefault(this)
            this.notify('SaveBinRuleAsDefault');
        end
        
        function [isValid, errMessage] = checkBinRulesValidity(~, selectedBinRule, selectedBinWidthRule,...
                numberOfBinsEditFieldValue, binWidthEditFieldValue, binBoundaryAtEditFieldValue)
            isValid = true;
            errMessage = {};
            switch selectedBinRule
                case 'NumberOfBins'
                    numberOfBinsEditFieldValue = strtrim(numberOfBinsEditFieldValue);
                    if isempty(numberOfBinsEditFieldValue)
                        errMessage{end+1} = iGetMessageString('stats:dfittool:error_noNumBins');
                    else
                        numberOfBinsValue = str2double(numberOfBinsEditFieldValue);
                        if isnan(numberOfBinsValue) || ~isequal(numel(numberOfBinsValue),1) ||...
                                ~iIsInteger(numberOfBinsValue) || numberOfBinsValue < 1 ||...
                                numberOfBinsValue > 1000
                            errMessage{end+1} = iGetMessageString('stats:dfittool:error_badNumBins');
                        end
                    end
                case 'BinWidth'
                    binWidthEditFieldValue = strtrim(binWidthEditFieldValue);
                    binBoundaryAtEditFieldValue = strtrim(binBoundaryAtEditFieldValue);
                    if isempty(binWidthEditFieldValue)
                        errMessage{end+1} = iGetMessageString('stats:dfittool:error_noBinWidth');
                    else
                        binWidthValue = str2double(binWidthEditFieldValue);
                        if isnan(binWidthValue) || binWidthValue <= 0
                            errMessage{end+1} = iGetMessageString('stats:dfittool:error_badBinWidth');
                        end
                    end
                    if strcmp(selectedBinWidthRule, 'BinWidthBinBoundaryAt')
                        if isempty(binBoundaryAtEditFieldValue)
                            errMessage{end+1} = iGetMessageString('stats:dfittool:error_noBinWidthBoundary');
                        else
                            binBoundaryAtValue = str2double(binBoundaryAtEditFieldValue);
                            if isnan(binBoundaryAtValue)
                                errMessage{end+1} = iGetMessageString('stats:dfittool:error_badBinBoundary');
                            end
                        end
                        errMessage = errMessage(~cellfun(@isempty, errMessage));
                    end
                otherwise
                    return;
            end
            if ~isempty(errMessage)
                isValid = false;
            end
        end
    end
end

% helper functions
function tf = iIsInteger(inputValue)
tf = ~mod(inputValue,1);
end

function outputStr = iGetMessageString(inputStr)
outputStr = getString(message(inputStr));
end
