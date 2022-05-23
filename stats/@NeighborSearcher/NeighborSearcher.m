classdef NeighborSearcher
%NeighborSearcher Neighbor search object.
%   NeighborSearcher is an abstract class used for nearest neighbor search
%   or radius search. You cannot create instances of this class directly.
%   You must create an instance of ExhaustiveSearcher or KDTreeSearcher.
%
%   NeighborSearcher properties:
%       X               - Data used to create the object.
%       Distance        - The distance metric.
%       DistParameter   - Additional parameter for the distance metric.
%
%   NeighborSearcher methods:
%       NeighborSearcher/knnsearch       - An abstract method
%       NeighborSearcher/rangesearch     - An abstract method
%
%   See also ExhaustiveSearcher, KDTreeSearcher, CREATENS, KNNSEARCH,
%   RANGESEARCH.

%   Copyright 2010-2011 The MathWorks, Inc.

    
   properties(GetAccess=protected, SetAccess=protected) 
       PrivDistance ='euclidean';
       PrivDistParameter =[];
   end
   
    properties(GetAccess=public, SetAccess=public, Dependent) 
        %Distance The distance metric.
        %   The Distance property is a string specifying the built-in
        %   distance metric that you used (applies for both
        %   ExhaustiveSearcher and KDTreeSearcher), or a function handle
        %   that you provided (applies only for ExhaustiveSearcher) when
        %   you created the object. This property saves the default
        %   distance metric used when you call the KNNSEARCH method to find
        %   nearest neighbors for future query points.
        %
        %   See also NeighborSearcher. 
        Distance ;
        
        %DistParameter Additional parameter for the distance metric.
        %   The DistParameter property specifies the additional parameter
        %   for the chosen distance metric.
        %                                   Value:
        % 
        %   if 'Distance' is 'minkowski'    A positive scalar indicating
        %                                   the exponent of the minkowski     
        %                                   distance. (Applies for both
        %                                   ExhaustiveSearcher and
        %                                   KDTreeSearcher.)
        %   if 'Distance' is 'mahalanobis'  A positive definite matrix        
        %                                   representing the covariance
        %                                   matrix used for computing the
        %                                   mahalanobis distance. (Applies
        %                                   only for ExhaustiveSearcher.)
        %   if 'Distance' is 'seuclidean'   A vector representing the scale
        %                                   value to use in computing the
        %                                   'seuclidean' distance. (Applies
        %                                   only for ExhaustiveSearcher.)
        %   otherwise                       Empty. 
        %
        %  See also NeighborSearcher. 
        DistParameter;
    end
       
    methods
        function this = set.Distance(this, distMetric)
            this = setDistance(this,distMetric);
        end
        
        function this = set.DistParameter(this, distPara)
            this = setDistParameter(this, distPara);
        end
        function distance = get.Distance(this)
            distance = this.PrivDistance;
        end
        
        function distPara = get.DistParameter(this)
            distPara = this.PrivDistParameter;
        end
    end
    
    methods (Access=protected,Abstract)
        this = setDistance(this, distMetric)
        this = setDistParameter(this, distpara)
    end
    
         
    properties(GetAccess=public, SetAccess=protected) % Public properties
        %X Data used to create the object.
        %  The X property is a matrix used to create the object.
        %
        %  See also NeighborSearcher. 
        X = [];
        
      
    end
   
    methods(Abstract)
         [idx,dist]=knnsearch(obj,y,varargin)
         %KNNSEARCH  A abstract method to find K nearest neighbors.
         
          [idx,dist]=rangesearch(obj,y,varargin)
         %RANGESEARCH  A abstract method for radius search. 
       
    end
    
    methods(Hidden, Static)
        function a = empty(varargin)
            error(message('stats:NeighborSearcher:NoEmptyAllowed'));
        end
    end
  
    methods(Hidden)
        function a = cat(varargin),        throwNoCatError(); end
        function a = horzcat(varargin),    throwNoCatError(); end
        function a = vertcat(varargin),    throwNoCatError(); end
        
%         function [varargout] = subsasgn(varargin)
%             %SUBSASGN Subscripted reference for a NeighborSearcher object.
%             %Subscript assignment is not allowed for a NeighborSearcher object         
%             error(message('stats:NeighborSearcher:subsasgn:NotAllowed'))
%         end
        
        
        function [varargout] = subsref(obj,s)
            %SUBSREF Subscripted reference for a NeighborSearcher object.
            %   B = SUBSREF(OBJ,S) is called for the syntax OBJ(S) when OBJ
            %   is a NeighborSearcher object. S is a structure array with
            %   the fields:
            %       type -- string containing '()', '{}', or '.' specifying
            %       the subscript type.
            %       subs -- Cell array or string containing the actual
            %       subscripts.
            
            
            switch s(1).type
                case '()'
                    error(message('stats:NeighborSearcher:subsref:ArraySubscript', class( obj )));
                case '{}'
                     error(message('stats:NeighborSearcher:subsref:CellReferenceNotAllowed'))    
                case '.'
                    methodsProp = [methods(obj);properties(obj)];
                    %The following is to prevent accessing private methods or
                    %properties
                    if ~any(strcmp(s(1).subs, methodsProp))
                        error(message('stats:NeighborSearcher:subsref:InvalidAccess', s( 1 ).subs, class( obj )));
                    end
                    [varargout{1:nargout}] = builtin('subsref',obj,s);
                  
                otherwise
                   [varargout{1:nargout}] = builtin('subsref',obj,s);
                 
            end
        end
    end
end   %classdef
%----------------------------------------------
function throwNoCatError()
error(message('stats:NeighborSearcher:NoCatAllowed'));
end

