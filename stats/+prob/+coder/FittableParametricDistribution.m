classdef FittableParametricDistribution < prob.coder.FittableDistribution &...
                                          prob.coder.TruncatableDistribution %#codegen %#internal
 
     
    % This classdef is significantly different from its MATLAB version.
    % This class could borrow its properties and methods from the following
    % classes under the "prob" package:
    % 1. ParametricDistribution
    % 2. ToolboxParametricDistribution
    % 3. ToolboxFittableParametricDistribution
    % 4. FittableParametricDistribution
    % This design is adopted to simplify the architecture for code
    % generation. 
    
    
    % 1. prob.ParametricDistribution properties
    
    properties(GetAccess='public',Constant=true,Abstract=true)
        %NumParameters Number of parameters in the distribution.
        %   NumParameters is a positive integer indicating the number of parameters
        %   in the distribution. This number includes all fixed and estimated
        %   parameters.
        %
        %   See also ParameterValues, ParameterNames, ParameterDescription.

%   Copyright 2019 The MathWorks, Inc.

        NumParameters
   end
    
    properties(GetAccess='public',SetAccess='protected',Abstract=true)
        %ParameterValues Values of the parameters in the distribution.
        %   ParametersValues is a vector containing values of the parameters in
        %   the distribution. The values may represent specified values if the
        %   distribution was created directly, or estimated values if the
        %   distribution was created by fitting to data.
        %
        %   See also ParameterNames, NumParameters, ParameterDescription.
        ParameterValues
    end
         
    % 2. prob.TruncatableParametricDistribution methods
    methods(Access=protected)
        function [varargout] = cdffun(this,x,varargin)
            pcell = convertNumToCell(this);
%             if nargout>=2
%                 covarg = {this.ParameterCovariance};
%             else
                covarg = {};
%             end
            [varargout{1:nargout}] = this.cdffunc(x,pcell{1:min(end,this.NumParameters)},covarg{:},varargin{:});
        end
        function [varargout] = icdffun(this,x,varargin)
            pcell = convertNumToCell(this);
%             if nargout>=2
%                 covarg = {this.ParameterCovariance};
%             else
                covarg = {};
%             end
            [varargout{1:nargout}] = this.invfunc(x,pcell{1:min(end,this.NumParameters)},covarg{:},varargin{:});
        end
        function y = pdffun(this,x)
            pcell = convertNumToCell(this);
            y = this.pdffunc(x,pcell{1:min(end,this.NumParameters)});
        end
    end
    
    methods(Static, Abstract = true)
        y = cdffunc(x,varargin)
        y = invfunc(x,varargin)
        y = pdffunc(x,varargin)
    end
    
    
    % 3. prob.ToolboxFittableParametricDistribution methods
      
    methods(Hidden)
 
        function pcell = convertNumToCell(this)
            pcell = cell(this.NumParameters,1);
            for idx = 1:this.NumParameters
                pcell{idx} = this.ParameterValues(idx);
            end
        end
    end
 
end % classdef


