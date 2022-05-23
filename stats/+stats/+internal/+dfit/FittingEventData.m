classdef FittingEventData < event.EventData
 % FITTINGEVENTDATA Data for a FitsManager event 
 
 %   Copyright 2019 The MathWorks, Inc.
  
    properties (SetAccess = 'private', GetAccess = 'public')
        Fit;
        FitFrame;
        OldName;
        NewName;
    end
    
    methods
         function data = FittingEventData(f,ff,on,nn)
            % FittingEventData Construct an instance of this class

            data.Fit = f;
            data.FitFrame = ff;
            data.OldName = on;
            data.NewName = nn;
        end
        
    end
end

