classdef EvaluateFunction < handle
    % EvaluateFunction    Interface for Evaluate window's probability
    % function to evaluate fits
    
    %   Copyright 2019 The MathWorks, Inc.
    
    methods(Abstract)
        [errmsg, evaluateTable] = evaluateFits(this,selectedFits, vectorValue, wantBounds, confLevel, plotFun);
        
        label = getAppropriateVectorTextFieldLabel(this);
        
        tf = canHaveConfidenceBoundOptionsEnabled(this); 
    end
end