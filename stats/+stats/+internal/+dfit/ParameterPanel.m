classdef ParameterPanel < handle
    % ParmeterPanel 
    
    % Copyright 2019 The MathWorks, Inc
    
    properties(Access = protected)
        FitType
        DisplayName
        Panel
    end
    
    methods

        function this = ParameterPanel()           
        end
  
        function ft = getFitType(this)
            ft = this.FitType;
        end
        
        function dName = getDisplayName(this)
            dName = this.DisplayName;
        end
            
        function pp = getParameterPanel(this)
            pp = this.Panel;
        end

    end
end

