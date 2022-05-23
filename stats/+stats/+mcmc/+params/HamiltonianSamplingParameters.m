classdef HamiltonianSamplingParameters < stats.mcmc.params.SamplingParameters
%

%   Copyright 2016 The MathWorks, Inc.
    
    properties(Constant,Hidden)
        SamplerTypeHMC                          = 'hmc';
        SamplerTypeQuasiNewtonHMC               = 'quasi-newton-hmc';
        AllowedSamplerTypes                     = {stats.mcmc.params.HamiltonianSamplingParameters.SamplerTypeHMC,...
                                                   stats.mcmc.params.HamiltonianSamplingParameters.SamplerTypeQuasiNewtonHMC};
        
        JitterMethodNone                        = 'none';
        JitterMethodJitterBoth                  = 'jitter-both';
        JitterMethodJitterNumSteps              = 'jitter-numsteps';
        AllowedJitterMethods                    = {stats.mcmc.params.HamiltonianSamplingParameters.JitterMethodNone,...
                                                   stats.mcmc.params.HamiltonianSamplingParameters.JitterMethodJitterBoth,...
                                                   stats.mcmc.params.HamiltonianSamplingParameters.JitterMethodJitterNumSteps};
                                 
        StepSizeTuningMethodNone                = 'none';
        StepSizeTuningMethodDualAveraging       = 'dual-averaging';
        AllowedStepSizeTuningMethods            = {stats.mcmc.params.HamiltonianSamplingParameters.StepSizeTuningMethodNone,...
                                                   stats.mcmc.params.HamiltonianSamplingParameters.StepSizeTuningMethodDualAveraging};
        
        MassVectorTuningMethodNone              = 'none';
        MassVectorTuningMethodHessian           = 'hessian';
        MassVectorTuningMethodIterativeSampling	= 'iterative-sampling';
        AllowedMassVectorTuningMethods         	= {stats.mcmc.params.HamiltonianSamplingParameters.MassVectorTuningMethodNone,...
                                                   stats.mcmc.params.HamiltonianSamplingParameters.MassVectorTuningMethodHessian,...
                                                   stats.mcmc.params.HamiltonianSamplingParameters.MassVectorTuningMethodIterativeSampling};        
        GammaAuto                             	= 'auto';
    end
    
    properties
        StepSize;               % Maximal step size for Leapfrog steps.
        NumSteps;               % Maximum number of Leapfrog steps.
        MassVector;             % Mass vector.
        JitterMethod;           % Jittering method for StepSize and NumSteps.
        DualAveragingGamma;     % Dual averaging Gamma.
        DualAveragingT0;        % Dual averaging T0.
        DualAveragingKappa;     % Dual averaging Kappa.
        
        SamplerType;            % HMC or quasi-Newton HMC sampler.
        
        HessianHistorySize;     % History size for quasi-Newton HMC sampler.
        Gamma;                  % Parameter for quasi-Newton HMC.
        DoJitterHistorySize;    % Should history size be jittered for quasi-Newton HMC?
        
        StepSizeTuningMethod;   % Method for tuning step size for HMC.
        MassVectorTuningMethod; % Method for tuning mass vector for HMC.
        
        CheckGradient;          % Should gradient of LOGPDF be checked when appropriate?
        
        IterativeSamplingIterations;                % Number of iterations of interleaved step size and mass vector tuning for mass vector tuning method 'iterative-sampling'.
        IterativeSamplingMinSamples;                % Number of samples to draw on the first iteration of mass vector tuning method 'iterative-sampling'.
        IterativeSamplingMaxSamples;                % Maximum number of samples to draw on any iteration of mass vector tuning method 'iterative-sampling'.
        IterativeSamplingInitialMassMultiplier;     % Number multiplying the mass vector to obtain the initial mass vector used in tuning method 'iterative-sampling'.
        IterativeSamplingNumStepsLimit;             % Maximum allowed number of leapfrog steps per HMC proposal during dual averaging while using tuning method 'iterative-sampling'.
        IterativeSamplingRegularizer;               % Amount added to each parameter's empirical variance in mass vector tuning method 'iterative-sampling'.
    end
    
    properties(GetAccess=public,SetAccess=protected)
        NumParameters;          % Number of parameters we are sampling.
        IsGammaAuto;            % True if Gamma is 'auto'.
    end
    
    methods
        function this = set.StepSize(this,stepsize)
            stepsize      = stats.mcmc.utils.InputValidator.validateStepSize(stepsize);
            this.StepSize = stepsize;
        end
        
        function this = set.NumSteps(this,numsteps)
            numsteps      = stats.mcmc.utils.InputValidator.validateNumSteps(numsteps);
            this.NumSteps = numsteps;
        end
        
        function this = set.MassVector(this,massvector)
            massvector      = checkMassVector(this,massvector);
            this.MassVector = massvector;
        end
        
        function this = set.JitterMethod(this,jittermethod)
            allowedJitterMethods = stats.mcmc.params.HamiltonianSamplingParameters.AllowedJitterMethods;
            jittermethod         = stats.mcmc.utils.InputValidator.validateJitterMethod(jittermethod,allowedJitterMethods);
            this.JitterMethod    = jittermethod;
        end
        
        function this = set.HessianHistorySize(this,historysize)
            historysize             = stats.mcmc.utils.InputValidator.validateHessianHistorySize(historysize);
            this.HessianHistorySize = historysize;
        end
        
        function this = set.SamplerType(this,samplertype)
            allowedSamplerTypes = stats.mcmc.params.HamiltonianSamplingParameters.AllowedSamplerTypes;
            samplertype         = stats.mcmc.utils.InputValidator.validateSamplerType(samplertype,allowedSamplerTypes);
            this.SamplerType    = samplertype;
        end
        
        function this = set.Gamma(this,gamma)
            [this,gamma] = checkGamma(this,gamma);            
            this.Gamma   = gamma; 
        end
        
        function this = set.DoJitterHistorySize(this,dojitterhistory)
            dojitterhistory          = stats.mcmc.utils.InputValidator.validateDoJitterHistory(dojitterhistory);
            this.DoJitterHistorySize = dojitterhistory;
        end
        
        function this = set.DualAveragingGamma(this,dualgamma)
            dualgamma               = stats.mcmc.utils.InputValidator.validateDualAveragingGamma(dualgamma);
            this.DualAveragingGamma = dualgamma;
        end
        
        function this = set.DualAveragingT0(this,dualt0)
            dualt0               = stats.mcmc.utils.InputValidator.validateDualAveragingT0(dualt0);
            this.DualAveragingT0 = dualt0;
        end
        
        function this = set.DualAveragingKappa(this,dualkappa)
            dualkappa               = stats.mcmc.utils.InputValidator.validateDualAveragingKappa(dualkappa);
            this.DualAveragingKappa = dualkappa;
        end
        
        function this = set.StepSizeTuningMethod(this,stepsizetuningmethod)
            allowedTuningMethods        = stats.mcmc.params.HamiltonianSamplingParameters.AllowedStepSizeTuningMethods;
            stepsizetuningmethod        = stats.mcmc.utils.InputValidator.validateStepSizeTuningMethod(stepsizetuningmethod,allowedTuningMethods);
            this.StepSizeTuningMethod	= stepsizetuningmethod;
        end
        
        function this = set.MassVectorTuningMethod(this,massvectuningmethod)
            allowedTuningMethods        = stats.mcmc.params.HamiltonianSamplingParameters.AllowedMassVectorTuningMethods;
            massvectuningmethod         = stats.mcmc.utils.InputValidator.validateMassVectorTuningMethod(massvectuningmethod,allowedTuningMethods);
            this.MassVectorTuningMethod = massvectuningmethod;
        end
        
        function this = set.IterativeSamplingIterations(this,iterativesamplingiterations)
            iterativesamplingiterations      = stats.mcmc.utils.InputValidator.validateIterativeSamplingIterations(iterativesamplingiterations);
            this.IterativeSamplingIterations = iterativesamplingiterations;
        end
        
        function this = set.IterativeSamplingMinSamples(this,iterativesamplingminsamples)
            iterativesamplingminsamples      = stats.mcmc.utils.InputValidator.validateIterativeSamplingMinSamples(iterativesamplingminsamples);
            this.IterativeSamplingMinSamples = iterativesamplingminsamples;
        end
        
        function this = set.IterativeSamplingMaxSamples(this,iterativesamplingmaxsamples)
            iterativesamplingmaxsamples      = stats.mcmc.utils.InputValidator.validateIterativeSamplingMaxSamples(iterativesamplingmaxsamples);
            this.IterativeSamplingMaxSamples = iterativesamplingmaxsamples;
        end
        
        function this = set.IterativeSamplingInitialMassMultiplier(this,iterativesamplinginitialmassmultiplier)
            iterativesamplinginitialmassmultiplier      = stats.mcmc.utils.InputValidator.validateIterativeSamplingInitialMassMultiplier(iterativesamplinginitialmassmultiplier);
            this.IterativeSamplingInitialMassMultiplier = iterativesamplinginitialmassmultiplier;
        end
        
        function this = set.IterativeSamplingNumStepsLimit(this,iterativesamplingnumstepslimit)
            iterativesamplingnumstepslimit      = stats.mcmc.utils.InputValidator.validateIterativeSamplingNumStepsLimit(iterativesamplingnumstepslimit);
            this.IterativeSamplingNumStepsLimit = iterativesamplingnumstepslimit;
        end
        
        function this = set.IterativeSamplingRegularizer(this,iterativesamplingregularizer)
            iterativesamplingregularizer      = stats.mcmc.utils.InputValidator.validateIterativeSamplingRegularizer(iterativesamplingregularizer);
            this.IterativeSamplingRegularizer = iterativesamplingregularizer;
        end
        
        function this = set.CheckGradient(this,checkgradient)
            checkgradient      = stats.mcmc.utils.InputValidator.validateCheckGradient(checkgradient);
            this.CheckGradient = checkgradient;
        end
    end
    
    methods
        function massvector = checkMassVector(this,massvector)
            p          = this.NumParameters;
            massvector = stats.mcmc.utils.InputValidator.validateMassVector(massvector,p);
        end
        
        function [this,gamma] = checkGamma(this,gamma)
            gammaAuto           = stats.mcmc.params.HamiltonianSamplingParameters.GammaAuto;
            [gamma,isGammaAuto] = stats.mcmc.utils.InputValidator.validateGamma(gamma,gammaAuto);
            this.IsGammaAuto    = isGammaAuto;
        end
    end
    
    methods
        function out = toStruct(this)
            out   = toStruct@stats.mcmc.params.SamplingParameters(this);
            names = excludedFieldNames(this);
            out   = rmfield(out,names);
        end
        
        function names = excludedFieldNames(this) %#ok<MANU>
            names = {'SamplerTypeHMC',...
                    'SamplerTypeQuasiNewtonHMC',...
                    'AllowedSamplerTypes',...
                    'JitterMethodNone',...
                    'JitterMethodJitterBoth',...
                    'JitterMethodJitterNumSteps',...
                    'AllowedJitterMethods',...
                    'StepSizeTuningMethodNone',...
                    'StepSizeTuningMethodDualAveraging',...
                    'AllowedStepSizeTuningMethods',...
                    'MassVectorTuningMethodNone',...
                    'MassVectorTuningMethodHessian',...
                    'MassVectorTuningMethodIterativeSampling',...
                    'AllowedMassVectorTuningMethods',...
                    'GammaAuto',...
                    'HessianHistorySize',...
                    'Gamma',...
                    'DoJitterHistorySize',...
                    'IsGammaAuto',...
                    'CheckGradient'};
        end
    end
    
    methods
        function this = HamiltonianSamplingParameters(logpdf,start,varargin)                         %#ok<INUSL>
            % 1.1 Set parameter defaults.
            P = length(start);
            
            dfltStepSize                                    = 0.1;
            dfltNumSteps                                    = 50;
            dfltMassVector                                  = ones(P,1);
            dfltJitterMethod                                = stats.mcmc.params.HamiltonianSamplingParameters.JitterMethodJitterBoth;
            dfltHessianHistorySize                          = 15;
            dfltSamplerType                                 = stats.mcmc.params.HamiltonianSamplingParameters.SamplerTypeHMC;
            dfltGamma                                       = 1;
            dfltDoJitterHistorySize                         = false;
            dfltDualAveragingGamma                          = 0.05;
            dfltDualAveragingT0                             = 10;
            dfltDualAveragingKappa                          = 0.75;
            dfltStepSizeTuningMethod                        = stats.mcmc.params.HamiltonianSamplingParameters.StepSizeTuningMethodDualAveraging;
            dfltMassVectorTuningMethod                      = stats.mcmc.params.HamiltonianSamplingParameters.MassVectorTuningMethodIterativeSampling;
            dfltIterativeSamplingIterations                 = 5;
            dfltIterativeSamplingMinSamples                 = 50;
            dfltIterativeSamplingMaxSamples                 = 200;
            dfltIterativeSamplingInitialMassMultiplier      = 100;
            dfltIterativeSamplingNumStepsLimit              = 150;
            dfltIterativeSamplingRegularizer                = eps(1);
            dfltCheckGradient                               = true;
            
            % 1.2 Parse optional name/value pairs.
            paramNames = {'StepSize', ...
                'NumSteps',...
                'MassVector',...
                'JitterMethod',...
                'HessianHistorySize',...
                'SamplerType',...
                'Gamma',...
                'DoJitterHistorySize',...
                'DualAveragingGamma',...
                'DualAveragingT0',...
                'DualAveragingKappa',...
                'StepSizeTuningMethod',...
                'MassVectorTuningMethod',...
                'IterativeSamplingIterations', ...
                'IterativeSamplingMinSamples',...
                'IterativeSamplingMaxSamples',...
                'IterativeSamplingInitialMassMultiplier',...
                'IterativeSamplingNumStepsLimit',...
                'IterativeSamplingRegularizer',...
                'CheckGradient'};
            paramDflts = {dfltStepSize,...
                dfltNumSteps,...
                dfltMassVector,...
                dfltJitterMethod,...
                dfltHessianHistorySize,...
                dfltSamplerType,...
                dfltGamma,...
                dfltDoJitterHistorySize,...
                dfltDualAveragingGamma,...
                dfltDualAveragingT0,...
                dfltDualAveragingKappa,...
                dfltStepSizeTuningMethod,...
                dfltMassVectorTuningMethod, ...
                dfltIterativeSamplingIterations,...
                dfltIterativeSamplingMinSamples,...
                dfltIterativeSamplingMaxSamples,...
                dfltIterativeSamplingInitialMassMultiplier,...
                dfltIterativeSamplingNumStepsLimit,...
                dfltIterativeSamplingRegularizer,...
                dfltCheckGradient};
            [stepsize,...
             numsteps,...
             massvector,...
             jittermethod,...
             historysize,...
             samplertype,...
             gamma,...
             dojitterhistory,...
             dualgamma,...
             dualt0,...
             dualkappa,...
             stepsizetuningmethod,...
             massvectuningmethod,...
             iterativesamplingiterations,...
             iterativesamplingminsamples,...
             iterativesamplingmaxsamples,...
             iterativesamplinginitialmassmultiplier,...
             iterativesamplingnumstepslimit,...
             iterativesamplingregularizer,...
             checkgradient] = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});
            
            % 1.3 Validate optional name/value pairs in set methods below.
            this.NumParameters                              = P;
            this.StepSize                                   = stepsize;
            this.NumSteps                                   = numsteps;
            this.MassVector                                 = massvector;
            this.JitterMethod                               = jittermethod;
            this.HessianHistorySize                         = historysize;
            this.SamplerType                                = samplertype;
            this.Gamma                                      = gamma;
            this.DoJitterHistorySize                        = dojitterhistory;
            this.DualAveragingGamma                         = dualgamma;
            this.DualAveragingT0                            = dualt0;
            this.DualAveragingKappa                         = dualkappa;
            this.StepSizeTuningMethod                       = stepsizetuningmethod;
            this.MassVectorTuningMethod                     = massvectuningmethod;
            this.IterativeSamplingIterations                = iterativesamplingiterations;
            this.IterativeSamplingMinSamples                = iterativesamplingminsamples;
            this.IterativeSamplingMaxSamples                = iterativesamplingmaxsamples;
            this.IterativeSamplingInitialMassMultiplier     = iterativesamplinginitialmassmultiplier;
            this.IterativeSamplingNumStepsLimit             = iterativesamplingnumstepslimit;
            this.IterativeSamplingRegularizer               = iterativesamplingregularizer;
            this.CheckGradient                              = checkgradient;
        end
    end
    
end