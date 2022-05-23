function [errmsg,x,values] = evaluatefits(fitNames,x,fun,wantBounds,confLevel,plotFun)
%EVALUATEFITS Evaluate fits for DFITTOOL


%   Copyright 1993-2020 The MathWorks, Inc.

% If the function is the empty string, clear the plot (if there is one) and
% clear any saved data.
if isempty(fitNames)
    stats.internal.dfit.evaluateplot(false); % closes the plot window
    stats.internal.dfit.dfgetset('evaluateResults', []);
    return
end
    
nfits = length(fitNames);
try
    x = sprintf('[ %s ]', x); % allow an unbracketed list of numbers to work
    x = evalin('base',x);
    if ~isnumeric(x)
        error(message('stats:dfittool:BadX'));
    end
    confLevel = evalin('base',confLevel) ./ 100;
    if ~isnumeric(confLevel) || ~isscalar(confLevel) || ~(confLevel>0 && confLevel<1) %#ok<*BDSCI>
        error(message('stats:dfittool:BadConfLevel'));
    end
catch ME
    x = [];
    values = zeros(0, nfits*(1+2*wantBounds));
    errmsg = getString(message('stats:dfstrings:sprintf_InvalidMATLABExpression',ME.message));
    return
end
errmsg = '';

x = x(:);
n = length(x);  

% Output table will have first column for data, then for each fit, one
% column for function, two columns for bounds (if requested).
values = NaN(n, nfits*(1+2*wantBounds));

fitdb = stats.internal.dfit.getfitdb;
for i = 1:nfits
    fit = find(fitdb, 'name', fitNames{i});
    % Cannot compute bounds for kernel smooths and certain parametric fits.
    getBounds = wantBounds && ...
        (~isequal(fit.fittype, 'smooth') && fit.distspec.hasconfbounds);

    % Evaluate the requested function for this fit.
    try
        if getBounds
            [y,ylo,yup] = eval(fit,x,fun,confLevel);
            values(:,3*i-2) = y;
            values(:,3*i-1) = ylo;
            values(:,3*i) = yup;
        else
            y = eval(fit,x,fun);
            if wantBounds
                values(:,3*i-2) = y;
            else
                values(:,i) = y;
            end
        end
    catch ME
       errmsg = getString(message('stats:dfstrings:sprintf_ErrorEvaluatingFit',ME.message));
    end
end

% Save the results for SaveToWorkSpace.  
dfgetset('evaluateResults', [x,values]);

% Save information about the fits that we've evaluated for the plot function.
dfgetset('evaluateFun', fun);
dfgetset('evaluateInfo', struct('fitNames',{fitNames}, 'wantBounds',wantBounds));

stats.internal.dfit.evaluateplot(plotFun);
