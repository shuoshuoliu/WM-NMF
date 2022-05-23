function err = dfcreateexclusionrule(name, dataset, yl, yh, yle, yhe)
%DFCREATEEXCLUSIONRULE Create dfittool exclusion rule


%   Copyright 2003-2004 The MathWorks, Inc.

err = '';
% check for duplicate name 
    outlierdb = getoutlierdb;
    outlier = down(outlierdb);

    while(~isempty(outlier))
        if strcmp(outlier.name, name)
            err = 'outliernamethesame';
            return;
        end
        outlier = right(outlier);
    end

% Make a new exclusion rule
stats.outlier(name, dataset, yl, yh, yle, yhe);



