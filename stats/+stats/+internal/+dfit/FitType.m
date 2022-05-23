classdef FitType
    % FitType
    
    % Copyright 2019 The MathWorks, Inc
          
    properties(Access = private)
        ModelType
        DisplayName
        RealName
        LowerLimit
        UpperLimit
        LowerLimitOK
        UpperLimitOK
        CensoredOK
        IntegerOnly
        ParameterNames = {}
        ParameterDescriptions = {}
        ParameterRequired
    end
    
    methods
        function this = FitType(mt, dName, rName, pnames, pdesc, preq, ...
                lolim, uplim, llOK, ulOK, cOK, intOnly)           
            this.ModelType = mt;
            this.DisplayName = dName;
            this.RealName = rName;
            this.LowerLimit = lolim;
            this.UpperLimit = uplim;
            this.LowerLimitOK = llOK;
            this.UpperLimitOK = ulOK;
            this.CensoredOK = cOK;
            this.IntegerOnly = intOnly;
            this.ParameterNames = pnames;
            this.ParameterDescriptions = pdesc;
            this.ParameterRequired = preq;
         end

        function dName = getDisplayName(this)
            dName = this.DisplayName;
        end
        
        function rName = getRealName(this)      
            rName = this.RealName;
        end

        function mt = getModelType(this)
            mt = this.ModelType;
        end

        function ll = getLowerLimit(this)
            ll = this.LowerLimit;
        end
        function ul = getUpperLimit(this)
            ul = this.UpperLimit;
        end
        function intO = getIntegerOnly(this)
            intO = this.IntegerOnly;
        end
        function cOK = getCensoredOK(this)
            cOK = this.CensoredOK;
        end
        function llOK = getLowerLimitOK(this)
            llOK = this.LowerLimitOK;
        end
        function ulOK = getUpperLimitOK(this)
            ulOK = this.UpperLimitOK;
        end
        function pNames = getParameterNames(this)
            pNames = this.ParameterNames;
        end
        function pDescs = getParameterDescriptions(this)
            pDescs = this.ParameterDescriptions;
        end
        function pReq = getParameterRequired(this)
            pReq = this.ParameterRequired;
        end
    end
end

