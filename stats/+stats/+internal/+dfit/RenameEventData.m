classdef RenameEventData < event.EventData
    % RenameEventData   Event data object passed to event listeners when
    % a UDD object is being renamed
    
    %   Copyright 2019 The MathWorks, Inc.
    
    properties (SetAccess = 'private', GetAccess = 'public')
        % Type   (string) UDD object class type being renamed
        ClassType 
        
        % OldName   (string) Original name of the UDD object being renamed
        OldName
        
        % NewName   (string) New name of the UDD object being renamed
        NewName
    end
    
    methods
        function this = RenameEventData(classType, oldName, newName)
            this.ClassType = classType;
            this.OldName = oldName;
            this.NewName = newName;
        end
    end
end