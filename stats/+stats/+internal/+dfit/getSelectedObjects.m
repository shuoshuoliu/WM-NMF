function selectedObjects = getSelectedObjects(db, table, selectedIndices, nameColumn)
% GETSELECTEDOBJECTS returns the udd object that is selected in table, 
% given the selectedIndices. In the table, nameColumn is the column that 
% has the object's name. 

% Copyright 2019-2020 The MathWorks, Inc.
    
    allSelectedRows = selectedIndices(:,1);
    uniqueRows = unique(allSelectedRows);

    for i = length(uniqueRows):-1:1
        selectedNames{i} = table.Data{uniqueRows(i), nameColumn};
    end
    for i = length(selectedNames):-1:1
        selectedObjects(i) = find(db, 'name', selectedNames{i});
    end
end
              