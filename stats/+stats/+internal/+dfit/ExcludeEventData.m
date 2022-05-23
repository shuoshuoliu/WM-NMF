classdef ExcludeEventData < event.EventData
% EXCLUDEEVENTDATA Data for an Exclude event 
 
%   Copyright 2019-2020 The MathWorks, Inc.

    properties
        OldName
        NewName
        DeleteNames
    end
    
    methods
        function this = ExcludeEventData(oldName, newName, deleteNames)
            this.OldName = oldName;
            this.NewName = newName;
            this.DeleteNames = deleteNames;
        end
    end
end