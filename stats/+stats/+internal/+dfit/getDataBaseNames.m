function names = getDataBaseNames(db)
% getDataBaseNames returns the given database object names

%   Copyright 2019-2020 The MathWorks, Inc.

    numItems = stats.internal.dfit.getNumDataBaseItems(db);

    if numItems > 0
        names = cell(1, numItems);
        c = 1;
        item = down(db);
        
        while(~isempty(item))
            names{1, c} = item.name;
            c = c+1;
            item = right(item);
        end
    else
        names = {};
    end
end
              