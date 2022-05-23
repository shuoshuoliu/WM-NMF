function dfdeleteexrule(names)
%DFDELETEEXRULE GUI helper to delete an exclusion rule


%   Copyright 2003-2011 The MathWorks, Inc.

fitdb = getfitdb;
fit = down(fitdb);
OKtoDelete = true;
fitsToDelete = {};

if ~isempty(fit)
    msg = '';
    for i=1:length(names)
        fitCnt = 0;
        fitnames = '';
        m='';
        fit = down(fitdb);
        while(~isempty(fit))
            if strcmp(names{i}, fit.exclusionrulename)
                fitsToDelete{1, end + 1} = fit.name;
                if fitCnt > 0
                    fitnames = [fitnames, ', '];
                end;
                fitCnt = fitCnt + 1;
                fitnames = [fitnames, fit.name];
            end
            fit = right(fit);
        end
        if fitCnt == 1
            m = getString(message('stats:dfstrings:sprintf_FitWillBeDeleted', names{i}, fitnames));
        elseif fitCnt > 1
            m = getString(message('stats:dfstrings:sprintf_FitsWillBeDeleted', names{i}, fitnames));
        end
        msg = [msg, m]; 
    end
    if length(msg) > 0
        button = questdlg(msg, getString(message('stats:dfstrings:dlg_DeletingExclusionRules')),...
                         getString(message('stats:dfstrings:button_Ok')), ...
                         getString(message('stats:dfstrings:button_Cancel')), ...
                         getString(message('stats:dfstrings:button_Ok')));
        if ~strcmp(button, getString(message('stats:dfstrings:button_Ok')))
            OKtoDelete = false;
        end
    end 
end


if OKtoDelete
    import com.mathworks.toolbox.stats.*;
    if ~isempty(fitsToDelete)
        FitsManager.getFitsManager.deleteFits(fitsToDelete);
    end
    OutliersManager.getOutliersManager.deleteOutliers(names);
end





