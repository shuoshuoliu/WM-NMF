classdef ChainInfo < classreg.learning.internal.DisallowVectorOps
%

%   Copyright 2016 The MathWorks, Inc.
    
    properties
        StartPoint;                     % Starting point to begin sampling.
        NumStepSizeTuningIterations;    % Number of tuning iterations.
        TargetAcceptanceRatio;          % Target acceptance ratio.
        NumStepsLimit;                  % Maximum NumSteps allowed during stepsize tuning.
        DoNumStepsWarn;                 % Whether to warn if NumSteps is clipped to NumStepsLimit.
        Burnin;                         % Number of burnin iterations (post tuning).
        NumSamples;                     % Number of samples to draw from the Markov chain.
        ThinSize;                       % If ThinSize is M then M-1 out of M values are omitted in the generated Markov chain.
        VerbosityLevel;                 % Non-negative integer specifying the verbosity level.
        NumPrint;                       % Positive integer specifying output frequency for verbose mode.
    end
    
    properties(GetAccess=public,SetAccess=protected)
        NumParameters;          % Number of parameters we are sampling.
    end
    
    methods
        function this = set.StartPoint(this,start)
            start           = checkStart(this,start);
            this.StartPoint = start;
        end
        
        function this = set.NumStepSizeTuningIterations(this,numtuningiter)
            numtuningiter                    = stats.mcmc.utils.InputValidator.validateNumStepSizeTuningIterations(numtuningiter);
            this.NumStepSizeTuningIterations = numtuningiter;
        end
        
        function this = set.TargetAcceptanceRatio(this,targetaccratio)
            targetaccratio             = stats.mcmc.utils.InputValidator.validateTargetAcceptanceRatio(targetaccratio);
            this.TargetAcceptanceRatio = targetaccratio;
        end
        
        function this = set.NumStepsLimit(this,numstepslimit)
            numstepslimit      = stats.mcmc.utils.InputValidator.validateNumStepsLimit(numstepslimit);
            this.NumStepsLimit = numstepslimit;
        end
        
        function this = set.DoNumStepsWarn(this,donumstepswarn)
            donumstepswarn      = stats.mcmc.utils.InputValidator.validateDoNumStepsWarn(donumstepswarn);
            this.DoNumStepsWarn = donumstepswarn;
        end
        
        function this = set.Burnin(this,burnin)
            burnin      = stats.mcmc.utils.InputValidator.validateBurnin(burnin);
            this.Burnin = burnin;
        end
        
        function this = set.NumSamples(this,numsamples)
            numsamples      = stats.mcmc.utils.InputValidator.validateNumSamples(numsamples);
            this.NumSamples = numsamples;
        end
        
        function this = set.ThinSize(this,thinSize)
            thinSize      = stats.mcmc.utils.InputValidator.validateThinSize(thinSize);
            this.ThinSize = thinSize;
        end
        
        function this = set.VerbosityLevel(this,verbose)
            verbose             = stats.mcmc.utils.InputValidator.validateVerbose(verbose);
            this.VerbosityLevel = verbose;
        end
        
        function this = set.NumPrint(this,numprint)
            numprint      = stats.mcmc.utils.InputValidator.validateNumPrint(numprint);
            this.NumPrint = numprint;
        end
    end
    
    methods
        function start = checkStart(this,start)
            p     = this.NumParameters;
            start = stats.mcmc.utils.InputValidator.validateStart(start,p);
        end
    end
    
    methods
        function out = toStruct(this)
            warning('off','MATLAB:structOnObject');
            out = struct(this);
            warning('on','MATLAB:structOnObject');
        end
    end
    
    methods
        function this = ChainInfo(state,varargin)
            % 1.1 Set parameter defaults.
            dfltStartPoint                  = state;
            dfltNumStepSizeTuningIterations = 0;
            dfltTargetAcceptanceRatio       = 0.65;
            dfltNumStepsLimit               = 2000;
            dfltDoNumStepsWarn              = true;
            dfltBurnin                      = 1000;
            dfltNumSamples                  = 1000;
            dfltThinSize                    = 1;
            dfltVerbosityLevel              = 0;
            dfltNumPrint                    = 100;
            
            % 1.2 Parse optional name/value pairs.
            paramNames = {  'StartPoint',   'NumStepSizeTuningIterations',   'TargetAcceptanceRatio',   'NumStepsLimit',   'DoNumStepsWarn',   'Burnin',   'NumSamples',   'ThinSize',   {'Verbose','VerbosityLevel'},   'NumPrint'};
            paramDflts = {dfltStartPoint, dfltNumStepSizeTuningIterations, dfltTargetAcceptanceRatio, dfltNumStepsLimit, dfltDoNumStepsWarn, dfltBurnin, dfltNumSamples, dfltThinSize,             dfltVerbosityLevel, dfltNumPrint};
            [start,...
             numtuningiter,...
             targetaccratio,...
             numstepslimit,...
             donumstepswarn,...
             burnin,...
             numsamples,...
             thin,...
             verbose,...
             numprint] = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});
            
            % 1.3 Validate optional name/value pairs in set methods below.
            this.StartPoint                     = start;
            this.NumStepSizeTuningIterations    = numtuningiter;
            this.TargetAcceptanceRatio          = targetaccratio;
            this.NumStepsLimit                  = numstepslimit;
            this.DoNumStepsWarn                 = donumstepswarn;
            this.Burnin                         = burnin;
            this.NumSamples                     = numsamples;
            this.ThinSize                       = thin;
            this.NumParameters                  = length(start);
            this.VerbosityLevel                 = verbose;
            this.NumPrint                       = numprint;
        end
    end
    
end