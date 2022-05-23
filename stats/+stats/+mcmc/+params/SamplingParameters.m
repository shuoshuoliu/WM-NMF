classdef (Abstract) SamplingParameters < classreg.learning.internal.DisallowVectorOps%
%

%   Copyright 2016 The MathWorks, Inc.
    
    methods
        function out = toStruct(this)
            warning('off','MATLAB:structOnObject');
            out = struct(this);
            warning('on','MATLAB:structOnObject');
        end
    end
    
end