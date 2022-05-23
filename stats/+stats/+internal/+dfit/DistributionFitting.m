classdef DistributionFitting < handle
    % DistributionFitting 
    
    % Copyright 2019-2020 The MathWorks, Inc

  properties (Access = private)
      FitTypesArray = {};
   end
   
   events
       SessionClearing
   end
     
   methods (Static)
        function this = getInstance()
            persistent Instance;
            if ~isa(Instance, 'stats.internal.dfit.DistributionFitting')  || ~isvalid(Instance)
                Instance = stats.internal.dfit.DistributionFitting();
            end
            this = Instance;
            mlock; 
        end
   end
    
   methods 
       function triggerSessionClearing(this)
           notify(this, 'SessionClearing');
       end
         
       function fitTypes = getFitTypes(this)
           fitTypes = this.FitTypesArray;
       end
       
       function fitType = getSelectedFitType(this, fullName)
           for i = 1:length(this.FitTypesArray)
               fitType = this.FitTypesArray{i};
               if strcmp(fitType.getDisplayName(), fullName)
                   return;
               end
           end
           % We should never get here.
           assert(false, getString(message('stats:dfittool:assert_getSelectedFitType')));
       end
       
       function clearFitTypes(this)
           this.FitTypesArray = {};
       end
       
       function addFitType(this, categoryName, fullName, codeName, pNames, pDescriptions, pRequired, lolim, uplim, loOK, upOK, censOK, intOnly)
           ft = stats.internal.dfit.FitType(categoryName, fullName, codeName, pNames, pDescriptions, pRequired, lolim, uplim, loOK, upOK, censOK, intOnly);
           this.FitTypesArray{end+1} = ft;        
       end
       
       function displayName = getFitTypeDisplayName(this, codeName)
           displayName = codeName;
           for i = 1:length(this.FitTypesArray)
               ft = this.FitTypesArray{i};
               if strcmp(codeName, ft.getRealName)
                   displayName = ft.getDisplayName();
                   break;
               end
           end
       end
   end
end


