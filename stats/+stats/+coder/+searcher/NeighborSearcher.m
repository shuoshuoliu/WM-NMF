classdef NeighborSearcher
    %#codegen
    %NeighborSearcher Neighbor search object for code generation.
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
    
    %   Copyright 2017 The MathWorks, Inc.
    
    properties (Abstract)
        X
        Distance
        DistParameter
    end
    
    methods(Abstract)
        
        [Idx,D] = knnsearch(obj,Y,varargin)
        %KNNSEARCH  A abstract method to find K nearest neighbors.
        
        [Idx,D] = rangesearch(obj,Y,r,varargin)
        %RANGESEARCH  A abstract method for radius search.
        
        outObj = matlabCodegenToRedirected(inObj)
        %MATLABCODEGENTOREDIRECTED An abstract static method for
        %redirection from MATLAB Class to Codegen class
        
        obj = fromStruct(str)
        %FROMSTRUCT An abstract static method for conversion from a
        %codegen compatible struct to Codegen class
    end
end