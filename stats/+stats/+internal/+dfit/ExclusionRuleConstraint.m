classdef ExclusionRuleConstraint < handle
    % ExclusionRuleConstraint   Helper methods used to check validity of
    % the parameters used to create an exclusion rule (outlier) object
    
    %   Copyright 2019-2020 The MathWorks, Inc.
    
    methods
        function [errMessage, messageTitle] = getErrorMessage(this, proposedName, lowerLimitTextFieldValue, upperLimitTextFieldValue, isLowerLimitLessEqual, isUpperLimitGreaterEqual)
            trimmedName = strtrim(proposedName);
            if isempty(proposedName)
                errMessage = getString(message('stats:dfittool:error_noName', iGetMessageString('stats:dfittool:label_exclusionRuleLowerCase')));
                messageTitle = getString(message('stats:dfittool:title_noObjectTypeName', iGetMessageString('stats:dfittool:label_exclusionRuleUpperCase')));
                return;
            end
            
            if isempty(trimmedName)
                errMessage = iGetMessageString('stats:dfittool:error_noOnlySpacesAllowed');
                messageTitle = getString(message('stats:dfittool:title_invalidObjectTypeName', iGetMessageString('stats:dfittool:label_exclusionRuleUpperCase')));
                return;
            end
            
            if iHasNoBoundValues(lowerLimitTextFieldValue, upperLimitTextFieldValue)
                errMessage = iGetMessageString('stats:dfittool:error_noExclusions');
                messageTitle = iGetMessageString('stats:dfittool:title_noExclusions');
                return;
            end
            
            [errMessage, messageTitle] = this.getErrorInvalidBoundValues(lowerLimitTextFieldValue, upperLimitTextFieldValue, isLowerLimitLessEqual, isUpperLimitGreaterEqual);
            if ~isempty(errMessage)
                return;
            end
            
            [errMessage, messageTitle] = this.getErrorProposedNameAlreadyExists(trimmedName);
        end
        
        function [errMessage, messageTitle] = getErrorProposedNameAlreadyExists(~, proposedName)
            errMessage = [];
            messageTitle = [];
            
            outlierdb = stats.internal.dfit.getoutlierdb;
            existingExclusionRuleNames = stats.internal.dfit.getDataBaseNames(outlierdb);
            
            if ismember(proposedName, existingExclusionRuleNames)
                errMessage = getString(message('stats:dfittool:error_duplicateOutlierName', proposedName));
                messageTitle = iGetMessageString('stats:dfittool:title_duplicateExclusionRuleName');
                return;
            end
        end
        
        function [errMessage, messageTitle] = getErrorInvalidBoundValues(this, lowerLimitTextFieldValue, upperLimitTextFieldValue, isLowerLimitLessEqual, isUpperLimitGreaterEqual)
            messageTitle = iMessageTitleBoundsError();
            [lowerLimitValue, lowerLimitError] = this.getLowerLimitValueAndErrorMessage(lowerLimitTextFieldValue, isLowerLimitLessEqual);
            [upperLimitValue, upperLimitError] = this.getUpperLimitValueAndErrorMessage(upperLimitTextFieldValue, isUpperLimitGreaterEqual);
            errMessage = {lowerLimitError; upperLimitError};
            errMessage = errMessage(~cellfun(@isempty, errMessage));
            if isempty(errMessage)
                if lowerLimitValue >= upperLimitValue
                    errMessage = iGetMessageString('stats:dfittool:error_excludeYBounds');
                end
            end            
        end
        
        function [lowerLimitValue, errMessage, messageTitle] = getLowerLimitValueAndErrorMessage(~, textFieldValue, ~)
            messageTitle = iMessageTitleBoundsError();
            errMessage = [];
            if isequal(textFieldValue, '-Inf') || isempty(textFieldValue)
                lowerLimitValue = -Inf;
            else
                lowerLimitValue = str2double(textFieldValue);
                if isnan(lowerLimitValue)
                    errMessage = iGetMessageString('stats:dfittool:error_excludeYLowerNotReal');
                end
            end
        end
        
        function [upperLimitValue, errMessage, messageTitle] = getUpperLimitValueAndErrorMessage(~, textFieldValue, ~)
            messageTitle = iMessageTitleBoundsError();
            errMessage = [];
            if isequal(textFieldValue, 'Inf') || isempty(textFieldValue)
                upperLimitValue = Inf;
            else
                upperLimitValue = str2double(textFieldValue);
                if isnan(upperLimitValue)
                    errMessage = iGetMessageString('stats:dfittool:error_excludeYUpperNotReal');
                end
            end
        end
    end
end

% Helpers
function outputStr = iGetMessageString(inputStr)
outputStr = getString(message(inputStr));
end

function tf = iHasNoBoundValues(lowerLimitTextFieldValue, upperLimitTextFieldValue)
tf = isempty(lowerLimitTextFieldValue) && isempty(upperLimitTextFieldValue);
end

function str = iMessageTitleBoundsError()
str = iGetMessageString('stats:dfittool:title_boundsError');
end