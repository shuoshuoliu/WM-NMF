classdef FittableDistribution %#codegen %#internal
%

%   Copyright 2019 The MathWorks, Inc.
    
    % Code generation version of FittableDistribution
    % Base class for fittable probability distributions.
    % Could serve as the base class for all distributions supported by code
    % generation version of fitdist, including nonparametric distributions
    % such as KernelDistribution
    
    properties( GetAccess='public',SetAccess='protected')
        %InputData Input data used in fitting distribution.
        %    InputData is a struct specifying the data used to fit the
        %    distribution. The struct has the following fields:
        %
        %      'Data'  Data vector 
        %      'Censoring'  Censoring vector, or empty if none
        %      'Frequency'  Frequency vector, or empty if none
        %
        %    See also fitdist.
        InputData = struct('Data',[],'Censoring',[],'Frequency',[]);
        
    end
    
    methods(Static,Hidden,Abstract)
        pd = fit(x,varargin);
    end

end% classdef
