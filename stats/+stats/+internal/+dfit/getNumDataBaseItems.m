function count = getNumDataBaseItems(db)
% GETNUMDATABASEITEMS returns the number of items in a distributionFitter database


% Copyright 2019-2020 The MathWorks, Inc.
    
count = 0;
    dbItem = down(db);
    while(~isempty(dbItem))
        count = count + 1;
        dbItem = right(dbItem);
    end
end
