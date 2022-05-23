classdef HamiltonianSampler < stats.mcmc.MCMCSampler
%HamiltonianSampler Hamiltonian Monte Carlo (HMC) sampler.
%   HamiltonianSampler is an object that represents the gradient based HMC
%   sampler. It can be used to sample from a PDF p(x) where x is an
%   unconstrained vector, i.e., every element of x can take any value in
%   (-Inf,Inf). HMC sampling requires a specification of log p(x) and its
%   gradient. If the actual density has constrained variables, these
%   constrained variables need to be transformed to unconstrained variables
%   prior to using the HMC sampler.
%
%   HamiltonianSampler properties:
%       StepSize               - Step size for each leapfrog step.
%       NumSteps               - Number of leapfrog steps used to generate new state proposals.
%       MassVector             - Mass vector for HMC sampling.
%       JitterMethod           - Method used to jitter StepSize and NumSteps in leapfrog iterations.
%       StepSizeTuningMethod   - Method used to tune StepSize for HMC.
%       MassVectorTuningMethod - Method used to tune MassVector for HMC.
%       LogPDF                 - Function handle representing the log PDF from which to generate samples.
%       StartPoint             - Start point for the Markov chain.
%       VariableNames          - Cell array of variable names for elements of StartPoint.
%
%   HamiltonianSampler methods:
%       estimateMAP - Compute maximum of LogPDF (MAP).
%       tuneSampler - Tune HMC sampler prior to drawing samples.
%       drawSamples - Draw samples from LogPDF.
%       diagnostics - Assess MCMC convergence.
%
%   See also hmcSampler, slicesample, mhsample.

%   Copyright 2016-2018 The MathWorks, Inc.
    
    properties(Constant,Hidden)
        SamplerName = 'HamiltonianSampler';        
    end
    
    properties(Hidden)
        %Impl - Implementation class for HMC sampler.
        %   The Impl property is of type HamiltonianImpl. It provides an
        %   implementation of HMC sampling. In particular, the Impl class
        %   provides methods for drawing samples and tuning MCMC
        %   parameters.
        Impl;
    end
    
    properties(Dependent)
        %StepSize - Step size for HMC leapfrog steps.
        %   The StepSize property is a scalar specifying the step size to
        %   use for HMC leapfrog steps.
        %
        %   See also NumSteps, MassVector, JitterMethod, StepSizeTuningMethod, tuneSampler, hmcSampler.
        StepSize;
        
        %NumSteps - Number of HMC leapfrog steps.
        %   The NumSteps property is a scalar specifying the number of HMC
        %   leapfrog steps to take when proposing a new state for the
        %   Markov chain.
        %
        %   See also StepSize, MassVector, JitterMethod, tuneSampler, hmcSampler.
        NumSteps;
        
        %MassVector - Mass vector for HMC leapfrog steps.
        %   The MassVector property is a P-by-1 vector specifying the
        %   "mass" of the momentum variables associated with variables of
        %   interest. P is equal to the length of property StartPoint.
        %
        %   See also StepSize, NumSteps, MassVectorTuningMethod, tuneSampler, hmcSampler.
        MassVector;
        
        %JitterMethod - Method to use to jitter StepSize and NumSteps for HMC.
        %   The JitterMethod property is a character vector specifying the
        %   method to use to vary the StepSize and NumSteps for HMC.
        %   Choices are:
        %
        %       JitterMethod         Meaning
        %       'none'            - Do not jitter StepSize and NumSteps.
        %       'jitter-both'     - Jitter both StepSize and NumSteps.
        %       'jitter-numsteps' - Jitter only NumSteps for fixed StepSize.
        %
        %   See also StepSize, NumSteps, hmcSampler.
        JitterMethod;

        %StepSizeTuningMethod - Tuning method for HMC StepSize.
        %   The StepSizeTuningMethod property is a character vector
        %   specifying the method to use to tune HMC StepSize. Choices are:
        %
        %   StepSizeTuningMethod     Meaning    
        %   'dual-averaging'       - StepSize is tuned for a fixed value of
        %                            simulation length (StepSize*NumSteps) 
        %                            to achieve a specified acceptance
        %                            ratio.
        %   'none'                 - StepSize is not tuned.
        %
        %   See also StepSize, MassVectorTuningMethod, tuneSampler, hmcSampler.
        StepSizeTuningMethod;
        
        %MassVectorTuningMethod - Tuning method for MassVector.
        %   The MassVectorTuningMethod property is a character vector
        %   specifying the method to use to tune MassVector. Choices are:
        %
        %   MassVectorTuningMethod   Meaning
        %   'iterative-sampling'   - MassVector is tuned via successive
        %                            approximation by drawing samples using
        %                            a sequence of MassVector estimates.
        %   'hessian'              - MassVector is set equal to the
        %                            negative diagonal Hessian of the
        %                            LogPDF at the specified starting point
        %                            for the Markov chain. For example,
        %                            estimateMAP can be used to compute
        %                            the MAP point and a Markov chain can
        %                            be started from that point with
        %                            MassVectorTuningMethod equal to
        %                            'hessian'.
        %   'none'                 - MassVector is not tuned.
        %
        %   See also MassVector, StepSizeTuningMethod, estimateMAP, tuneSampler, hmcSampler.
        MassVectorTuningMethod;        
    end
    
    methods
        function this = set.StepSize(this,stepsize)
            this.Impl.SamplerParameters.StepSize = stepsize;
        end
        
        function stepsize = get.StepSize(this)
            stepsize = this.Impl.SamplerParameters.StepSize;
        end
        
        function this = set.NumSteps(this,numsteps)
            this.Impl.SamplerParameters.NumSteps = numsteps;
        end
        
        function numsteps = get.NumSteps(this)
            numsteps = this.Impl.SamplerParameters.NumSteps;
        end
        
        function this = set.MassVector(this,massvec)
            this.Impl.SamplerParameters.MassVector = massvec;
        end
        
        function massvec = get.MassVector(this)
            massvec = this.Impl.SamplerParameters.MassVector;
        end
        
        function this = set.JitterMethod(this,jittermethod)
            jittermethod = convertStringsToChars(jittermethod);
            this.Impl.SamplerParameters.JitterMethod = jittermethod;
        end
        
        function jittermethod = get.JitterMethod(this)
            jittermethod = this.Impl.SamplerParameters.JitterMethod;
        end
        
        function this = set.StepSizeTuningMethod(this,stepsizetuningmethod)
            stepsizetuningmethod = convertStringsToChars(stepsizetuningmethod);
            this.Impl.SamplerParameters.StepSizeTuningMethod = stepsizetuningmethod;
        end
        
        function stepsizetuningmethod = get.StepSizeTuningMethod(this)
            stepsizetuningmethod = this.Impl.SamplerParameters.StepSizeTuningMethod;
        end
        
        function this = set.MassVectorTuningMethod(this,massvectuningmethod)
            massvectuningmethod = convertStringsToChars(massvectuningmethod);
            this.Impl.SamplerParameters.MassVectorTuningMethod = massvectuningmethod;
        end
        
        function massvectuningmethod = get.MassVectorTuningMethod(this)
            massvectuningmethod = this.Impl.SamplerParameters.MassVectorTuningMethod;
        end
    end
    
    methods(Hidden)
        function this = HamiltonianSampler(logpdf,start,varargin)
%HamiltonianSampler - Create a Hamiltonian Monte Carlo (HMC) sampler.
%   HMC = HamiltonianSampler(LOGPDF,STARTPOINT) creates a
%   HamiltonianSampler object that can be used for Hamiltonian Monte Carlo
%   (HMC) sampling. LOGPDF is a function handle that evaluates the log
%   probability density (up to an additive constant) of the desired
%   equilibrium distribution. STARTPOINT is a P-by-1 vector representing
%   the initial point from which to start HMC sampling.
%
%   The HamiltonianSampler object can be used to generate a Markov chain
%   using HMC whose equilibrium distribution corresponds to LOGPDF. The
%   function handle LOGPDF should be callable like this:
%
%       [LPDF,GLPDF] = LOGPDF(X)
%
%   where X is a P-by-1 vector, LPDF is a scalar representing the log
%   probability density (up to an additive constant) and GLPDF is a P-by-1
%   vector containing the gradient of LOGPDF at X.
%
%   If UseNumericalGradient name/value pair (see below) is true then LOGPDF
%   should be callable like this:
%
%       LPDF = LOGPDF(X)
%
%   where X is a P-by-1 vector and LPDF is a scalar representing the log
%   probability density (up to an additive constant). In this case, the
%   gradient of LOGPDF at X is approximated numerically. Using numerical
%   gradient makes sampling more expensive compared to using analytical
%   gradient.
%
%   Once you create a HamiltonianSampler object, you can:
%
%   o Tune HMC sampler using the tuneSampler method.
%
%   o Draw samples using the drawSamples method.
%
%   o Estimate the maximum of LOGPDF using the estimateMAP method.
%
%   o Assess convergence using the diagnostics method.
%
%   HMC = HamiltonianSampler(LOGPDF,STARTPOINT,'Name','Value',...)
%   specifies additional name/value pairs like this:
%
%       Name            Value
%       'StepSize'   -  A scalar specifying the step size to use for
%                       Hamiltonian dynamics. Default is 0.1.
%       'NumSteps'   -  Number of steps of Hamiltonian dynamics to use for
%                       generating proposals. Default is 50.
%       'MassVector' -  A P-by-1 vector specifying the "mass" of the
%                       momentum variables associated with variables of
%                       interest. Default is ones(P,1) where P is the
%                       length of STARTPOINT.
%       'JitterMethod'
%                    -  A character vector specifying the method to use to
%                       jitter the StepSize and NumSteps for HMC. Choices
%                       are 'none', 'jitter-both' and 'jitter-numsteps'.
%                       When 'JitterMethod' is 'jitter-both', StepSize and
%                       NumSteps are randomly jittered for each leapfrog
%                       trajectory. When 'JitterMethod' is
%                       'jitter-numsteps', only NumSteps is randomly
%                       jittered. Default is 'jitter-both'.
%       'StepSizeTuningMethod'
%                    -  A character vector specifying the method to use to
%                       tune HMC StepSize. Choices are 'dual-averaging' or
%                       'none'. Default is 'dual-averaging'.
%       'MassVectorTuningMethod'
%                    -  A character vector specifying the method to use to
%                       tune MassVector. Choices are 'iterative-sampling',
%                       'hessian' or 'none'. Default is
%                       'iterative-sampling'.
%       'CheckGradient'
%                    -  A logical scalar that is either true (or 1) or
%                       false (or 0) specifying whether to check the
%                       analytical gradient of LOGPDF at STARTPOINT vs the
%                       numerical gradient or not. Default is true.
%       'VariableNames'
%                    -  A cell array of length P specifying variable names
%                       where P is the length of STARTPOINT. The name for
%                       element i of STARTPOINT is taken from element i of
%                       VariableNames. Default is to use variable names
%                       like x1, x2 etc.
%       'UseNumericalGradient'
%                    -  A logical scalar that is either true (or 1) or
%                       false (or 0) that specifies whether to use the
%                       numerical gradient of LOGPDF or not. If true,
%                       LOGPDF is called like LPDF = LOGPDF(X) so that
%                       LOGPDF need not return gradient as the second
%                       output. Default is false.

            % 1. Parse optional name/value pairs.
            dfltVariableNames        = [];
            dfltUseNumericalGradient = false;

            paramNames = {  'VariableNames',   'UseNumericalGradient'};
            paramDflts = {dfltVariableNames, dfltUseNumericalGradient};

            [variableNames,useNumericalGradient,~,otherArgs] = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});

            % 2. Validate VariableNames.
            if isempty(variableNames)
                variableNames = stats.mcmc.MCMCSampler.getDefaultVariableNames(length(start));
            else
                variableNames = stats.mcmc.utils.InputValidator.validateVariableNames(variableNames,length(start));
            end

            % 3. Validate UseNumericalGradient.
            useNumericalGradient = stats.mcmc.utils.InputValidator.validateUseNumericalGradient(useNumericalGradient);
            
            % 4. Validate inputs logpdf and start. If we are using the
            % numerical gradient, wrap the logpdf to return numerical
            % gradient.
            if useNumericalGradient
                [logpdf,start] = stats.mcmc.HamiltonianSampler.processLogPDFAndStartNumericalGradient(logpdf,start);
            else
                [logpdf,start] = stats.mcmc.HamiltonianSampler.processLogPDFAndStartAnalyticalGradient(logpdf,start);
            end
            
            % 5. Name of this sampler?
            name = stats.mcmc.HamiltonianSampler.SamplerName;

            % 6. Construct an object to hold sampling parameters.
            samplerParameters = stats.mcmc.params.HamiltonianSamplingParameters(logpdf,start,otherArgs{:},'SamplerType',stats.mcmc.params.HamiltonianSamplingParameters.SamplerTypeHMC);

            % 7. This is a gradient based sampler and so we have access to
            % the gradient of logpdf.
            haveGrad = true;
            
            % 8. Make the Impl object.
            this.Impl = stats.mcmc.impl.HamiltonianImpl(logpdf,haveGrad,start,variableNames,name,samplerParameters);
        end
    end
    
    methods
        function [xsmpl,state,accratio,stepSizeTuningInfo] = drawSamples(this,varargin)
%drawSamples - Generate Markov chain using Hamiltonian Monte Carlo (HMC).
%   [XSMPL,ENDPOINT,ACCRATIO] = drawSamples(HMC) takes an object HMC of
%   type HamiltonianSampler and draws samples from it with the default
%   settings. XSMPL is a N-by-P matrix containing samples from the Markov
%   chain where each row of XSMPL is one sample. ENDPOINT is a P-by-1
%   vector containing the final state of the Markov chain after drawing
%   samples. ACCRATIO is a scalar containing the fraction of proposals that
%   were accepted including the burnin period.
%
%   [...] = drawSamples(HMC,'Name','Value',...) specifies additional
%   name/value pairs like this:
%
%       Name            Value
%       'Burnin'     -  A scalar integer specifying the number of initial
%                       samples from the Markov chain that should be
%                       discarded before inclusion in the output XSMPL.
%                       Default is 1000.
%       'NumSamples' -  An integer specifying the number of samples to draw
%                       from the HMC Markov chain after the 'Burnin' period
%                       is over. Default is 1000.
%       'ThinSize'   -  A positive integer specifying the thinning size. If
%                       ThinSize is M then M-1 out of M values are omitted
%                       in the generated Markov chain. Default is 1.
%       'StartPoint' -  A P-by-1 vector specifying the initial point to 
%                       start sampling where P is the length of
%                       HMC.StartPoint. Default is HMC.StartPoint.
%       'VerbosityLevel'
%                    -  A non-negative integer specifying the verbosity
%                       level as follows:
%                          * 0    - no iteration summary is displayed.
%                          * >=1  - iteration summary is displayed on
%                                   screen.
%                       Default is 0.
%       'NumPrint'   -  A positive integer specifying the frequency with
%                       which to output iteration summary when
%                       'VerbosityLevel' is >=1. 'NumPrint' iterations are
%                       processed for every line of iteration summary.
%                       Default is 100.
%
%   See also estimateMAP, tuneSampler, diagnostics, hmcSampler.


%       'NumStepSizeTuningIterations'
%                    -  A scalar integer specifying the number of tuning
%                       iterations before burnin during which the StepSize
%                       parameter is optimized if StepSizeTuningMethod is
%                       set to 'dual-averaging'. Default is 0.
%       'TargetAcceptanceRatio'     
%                    -  A real scalar between 0 and 1 specifying the target
%                       acceptance rate for HMC. If StepSizeTuningMethod is
%                       'dual-averaging' then StepSize is adjusted before
%                       burnin to achieve an MCMC acceptance rate equal to
%                       TargetAcceptanceRatio. Default is 0.65.
%
%   o STEPSIZETUNINGINFO is a structure containing step size tuning
%   information if 'NumStepSizeTuningIterations' > 0.
%   STEPSIZETUNINGINFO has these fields:
%       FieldName           Meaning
%       StepSize          - Tuned StepSize for HMC.
%       NumSteps          - Tuned NumSteps for HMC.
%       StepSizeProfile   - Evolution of StepSize during tuning.
%       AcceptanceRatio   - The final acceptance ratio achieved during
%                           StepSize tuning.
%
%   If 'NumStepSizeTuningIterations' is 0 then StepSizeProfile is empty.
%   More elaborate tuning of HMC parameters can be done using the
%   tuneSampler method.

            % 1.1 Set parameter defaults.            
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            chainInfo = stats.mcmc.params.ChainInfo(this.StartPoint);

            dfltBurnin                      = 1000;
            dfltNumSamples                  = 1000;
            dfltThinSize                    = 1;
            dfltStartPoint                  = this.StartPoint;
            dfltNumStepSizeTuningIterations	= 0;
            dfltTargetAcceptanceRatio       = 0.65;
            dfltVerbosityLevel              = 0;
            dfltNumPrint                    = 100;
            
            % 1.2 Parse optional name/value pairs.
            paramNames = {  'Burnin',   'NumSamples',   'ThinSize',   'StartPoint',   'NumStepSizeTuningIterations',   'TargetAcceptanceRatio',   {'Verbose','VerbosityLevel'},   'NumPrint'};
            paramDflts = {dfltBurnin, dfltNumSamples, dfltThinSize, dfltStartPoint, dfltNumStepSizeTuningIterations, dfltTargetAcceptanceRatio,             dfltVerbosityLevel, dfltNumPrint};
            [burnin,...
             numsamples,...
             thin,...
             start,...
             numtuningiter,...
             targetaccratio,...
             verbose,...
             numprint] = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});
         
            % 2. Update chainInfo object. This object should do validation
            % in its set methods.
            chainInfo.Burnin                        = burnin;
            chainInfo.NumSamples                    = numsamples;
            chainInfo.ThinSize                      = thin;
            chainInfo.StartPoint                    = start;
            chainInfo.NumStepSizeTuningIterations   = numtuningiter;
            chainInfo.TargetAcceptanceRatio         = targetaccratio;
            chainInfo.VerbosityLevel                = verbose;
            chainInfo.NumPrint                      = numprint;
            
            % 3. Call drawSamples on the Impl object.
            [xsmpl,state,accratio,stepSizeTuningInfo] = drawSamples(this.Impl,chainInfo,this.Impl.SamplerParameters);
            
            % 4. Return xsmpl as a N-by-P matrix where N is the number of
            % samples and P is the number of variables.
            xsmpl = xsmpl';
        end
        
        function [this,tuningInfo] = tuneSampler(this,varargin)
%tuneSampler - Tune Hamiltonian Monte Carlo (HMC) parameters.
%   [TUNEDHMC,TUNINGINFO] = tuneSampler(HMC) takes object HMC of type
%   HamiltonianSampler and tunes parameters for subsequent use in
%   drawSamples. The returned object TUNEDHMC is a tuned HMC sampler. In
%   the tuning phase, first the MassVector is tuned according to the
%   specified MassVectorTuningMethod. Then the StepSize is tuned using the
%   specified StepSizeTuningMethod.
%   TUNINGINFO is a structure containing additional tuning information, 
%   with these fields:
%       FieldName            Meaning
%       MassVector           - A P-by-1 vector containing the tuned mass vector.   
%       StepSize             - Tuned StepSize.
%       NumSteps             - Tuned NumSteps.
%       MassVectorTuningInfo - A structure containing information about the
%                              specific MassVector tuning method used. 
%           If MassVectorTuningMethod is 'iterative-sampling', fields are:
%               MassVector   - The tuned mass vector.
%               IterativeSamplingMassVectorProfile 
%                            - A P-by-K matrix of mass vectors used during
%                              the K iterations.
%               IterativeSamplingNumSamples
%                            - A K-by-1 vector of the number of samples
%                              drawn in each of the K iterations.
%           If MassVectorTuningMethod is 'hessian', fields are:
%               MassVector   - The tuned mass vector.
%               NegativeDiagonalHessian 
%                            - A P-by-1 vector containing the negative
%                              diagonal Hessian. This may differ from the
%                              MassVector if some elements are negative.
%               HessianPoint - A P-by-1 vector containing the point at
%                              which the Hessian was evaluated.
%           If MassVectorTuningMethod is 'none' then MassVectorTuningInfo
%           is empty.
%       StepSizeTuningInfo   - A structure containing information about the
%                              specific StepSize tuning method used.
%           If StepSizeTuningMethod is 'dual-averaging', fields are:
%               StepSize     - The final tuned StepSize.
%               NumSteps     - The final tuned NumSteps.
%               StepSizeProfile
%                            - A NumStepSizeTuningIterations-by-1 vector
%                              of StepSize estimates computed while tuning.
%               AcceptanceRatio
%                            - The final acceptance ratio achieved during
%                              StepSize tuning.
%           If 'NumStepSizeTuningIterations' is 0 or 'StepSizeTuningMethod'
%           is 'none' then StepSizeTuningInfo is empty. 
%   TUNEDHMC.StartPoint is a P-by-1 vector containing the final state of
%   the Markov chain after tuning parameters.
%
%   [...] = tuneSampler(HMC,'Name','Value',...) specifies additional
%   name/value pairs like this:
%
%       Name            Value
%       'StartPoint' -  A P-by-1 vector specifying the initial point to 
%                       start tuning where P is the length of
%                       HMC.StartPoint. Default is HMC.StartPoint.
%       'StepSizeTuningMethod'     
%                    -  A character vector or string specifying the method 
%                       to use to tune HMC StepSize. Choices are 
%                       'dual-averaging' or 'none'. Default is 
%                       HMC.StepSizeTuningMethod.
%       'MassVectorTuningMethod'
%                    -  A character vector or string specifying the method
%                       to use to tune MassVector. Choices are 
%                       'iterative-sampling', 'hessian' and 'none'. If 
%                       'iterative-sampling', the MassVector is tuned via 
%                       successive approximation by drawing samples using a 
%                       sequence of MassVector estimates. If 'hessian', the 
%                       MassVector is set equal to the negative diagonal 
%                       Hessian of LogPDF at the value specified in 
%                       'StartPoint'. If 'none', the MassVector is not 
%                       tuned. Default is HMC.MassVectorTuningMethod.
%       'NumStepSizeTuningIterations'     
%                    -  A scalar integer specifying the number of tuning
%                       iterations during which the StepSize parameter is 
%                       adjusted. Applies only if StepSizeTuningMethod is 
%                       'dual-averaging'. Default is 100.
%       'TargetAcceptanceRatio'     
%                    -  A real scalar between 0 and 1 specifying the target
%                       acceptance rate for MCMC during StepSize tuning.
%                       Applies only if StepSizeTuningMethod is
%                       'dual-averaging'. Default is 0.65.
%       'NumStepsLimit'
%                    -  A positive integer specifying the maximum number of
%                       leapfrog steps allowed during StepSize tuning.
%                       Applies only if StepSizeTuningMethod is
%                       'dual-averaging'. Default is 2000.
%       'VerbosityLevel'
%                    -  A non-negative integer specifying the verbosity
%                       level as follows:
%                          * 0    - No iteration summary is displayed.
%                          * 1    - A summary of stepsize tuning iterations
%                                   is displayed on screen.
%                          * 2    - Summaries of mass vector tuning
%                                   iterations (if any) and stepsize tuning
%                                   iterations are displayed on screen.
%                       Default is 0.
%       'NumPrint'   -  A positive integer specifying the frequency with
%                       which to output iteration summary when
%                       'VerbosityLevel' is >=1. 'NumPrint' iterations are
%                       processed for every line of iteration summary.
%                       Default is 100.
%
%   See also estimateMAP, drawSamples, diagnostics, hmcSampler.

            % 1.1 Set parameter defaults.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            samplerParams = this.Impl.SamplerParameters;
            chainInfo     = stats.mcmc.params.ChainInfo(this.StartPoint);

            dfltStartPoint                  = this.StartPoint;
            dfltNumStepSizeTuningIterations	= 100;
            dfltTargetAcceptanceRatio       = 0.65;
            dfltNumStepsLimit               = 2000;
            dfltStepSizeTuningMethod        = samplerParams.StepSizeTuningMethod;
            dfltMassVectorTuningMethod      = samplerParams.MassVectorTuningMethod;
            dfltVerbosityLevel              = 0;
            dfltNumPrint                    = 100;
            
            % 1.2 Parse optional name/value pairs.
            paramNames = {  'StartPoint',   'NumStepSizeTuningIterations',   'TargetAcceptanceRatio',   'NumStepsLimit',   'StepSizeTuningMethod',   'MassVectorTuningMethod',   {'Verbose','VerbosityLevel'},   'NumPrint'};
            paramDflts = {dfltStartPoint, dfltNumStepSizeTuningIterations, dfltTargetAcceptanceRatio, dfltNumStepsLimit, dfltStepSizeTuningMethod, dfltMassVectorTuningMethod,             dfltVerbosityLevel, dfltNumPrint};
            [start,...
             numtuningiter,...
             targetaccratio,...
             numstepslimit,...
             stepsizetuningmethod,...
             massvectuningmethod,...
             verbose,...
             numprint] = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});
         
            % 2. Update chainInfo object. This object should do validation
            % in its set methods.            
            chainInfo.StartPoint                    = start;
            chainInfo.NumStepSizeTuningIterations   = numtuningiter;
            chainInfo.TargetAcceptanceRatio         = targetaccratio;
            chainInfo.NumStepsLimit                 = numstepslimit;
            chainInfo.VerbosityLevel                = verbose;
            chainInfo.NumPrint                      = numprint;
            
            % 3. Update sampling parameters object. This object should do
            % validation in its set methods.
            samplerParams.StepSizeTuningMethod   = stepsizetuningmethod;
            samplerParams.MassVectorTuningMethod = massvectuningmethod;            

            % 4. Call tuneSampler on the Impl object.
            [this.Impl,tuningInfo] = tuneSampler(this.Impl,chainInfo,samplerParams);
        end
    end
    
    methods(Static,Hidden)
        function [logpdf,start] = processLogPDFAndStartAnalyticalGradient(logpdf,start)
            % 1. Validate logpdf.
            logpdf = stats.mcmc.utils.InputValidator.validateLogPDF(logpdf);
            
            % 2. Assume logpdf accepts row vectors if start is a row vector.
            logpdfAcceptsRowVector = isrow(start);
            
            % 3. Validate start and make it into a column vector.
            start = stats.mcmc.utils.InputValidator.validateStart(start);
            
            % 4. Ensure that logpdf can be validly called at start.
            if logpdfAcceptsRowVector
                [~,~,isGradientRowVector] = stats.mcmc.utils.InputValidator.validateLogPDFAndStartAnalyticalGradient(logpdf,start');
            else
                [~,~,isGradientRowVector] = stats.mcmc.utils.InputValidator.validateLogPDFAndStartAnalyticalGradient(logpdf,start);
            end
            
            % 5. Wrap logpdf so it accepts column vectors and returns
            % gradients as column vectors.
            if ( logpdfAcceptsRowVector || isGradientRowVector )
                logpdf = stats.mcmc.HamiltonianSampler.prepareLogPDFAnalyticalGradient(logpdf,logpdfAcceptsRowVector,isGradientRowVector);
            end
        end
        
        function z = prepareLogPDFAnalyticalGradient(logpdf,logpdfAcceptsRowVector,isGradientRowVector)
            % Assume that logpdf can be called like this:
            %   [f,g] = logpdf(start)
            % where:
            %   start = row vector if logpdfAcceptsRowVector is true and column vector otherwise.
            %   f     = scalar output.
            %   g     = row vector if isGradientRowVector is true and column vector otherwise.
            %
            % We would like to transform logpdf such that we can call it
            % like this:
            %   [f,g] = logpdf(start)
            % where:
            %   start = column vector.
            %   f     = scalar output.
            %   g     = column vector.
            
            z = @fun;
            function [f,g] = fun(start)
                if logpdfAcceptsRowVector
                    start = start';
                end
                
                if nargout < 2
                    f = logpdf(start);
                else
                    [f,g] = logpdf(start);
                    if isGradientRowVector
                        g = g(:);
                    end
                end
            end
        end
        
        function [logpdf,start] = processLogPDFAndStartNumericalGradient(logpdf,start)
            % 1. Validate logpdf.
            logpdf = stats.mcmc.utils.InputValidator.validateLogPDF(logpdf);
            
            % 2. Assume logpdf accepts row vectors if start is a row vector.
            logpdfAcceptsRowVector = isrow(start);
            
            % 3. Validate start and make it into a column vector.
            start = stats.mcmc.utils.InputValidator.validateStart(start);
            
            % 4. Ensure that logpdf can be validly called at start. When
            % using numerical gradient, logpdf should return a scalar
            % output.
            if logpdfAcceptsRowVector
                stats.mcmc.utils.InputValidator.validateLogPDFAndStartNumericalGradient(logpdf,start');
            else
                stats.mcmc.utils.InputValidator.validateLogPDFAndStartNumericalGradient(logpdf,start);
            end
            
            % 5. Wrap logpdf so it accepts column vectors.
            if logpdfAcceptsRowVector
                logpdf = stats.mcmc.HamiltonianSampler.prepareLogPDFNumericalGradient(logpdf,logpdfAcceptsRowVector);
            end
            
            % 6. Finally, wrap logpdf to return numerical gradient.
            logpdf = stats.mcmc.HamiltonianSampler.wrapLogPDFToReturnNumericalGradient(logpdf);
        end
        
        function z = prepareLogPDFNumericalGradient(logpdf,logpdfAcceptsRowVector)
            % Assume that logpdf can be called like this:
            %   f = logpdf(start)
            % where:
            %   start = row vector if logpdfAcceptsRowVector is true and column vector otherwise.
            %   f     = scalar output.
            %
            % We would like to transform logpdf such that we can call it
            % like this:
            %   f = logpdf(start)
            % where:
            %   start = column vector.
            %   f     = scalar output.
            
            z = @fun;
            function f = fun(start)
                if logpdfAcceptsRowVector
                    start = start';
                end
                f = logpdf(start);
            end
        end
        
        function z = wrapLogPDFToReturnNumericalGradient(logpdf)
            % Assume that logpdf can be called like this:
            %   f = logpdf(x)
            % where:
            %   x = column vector.
            %   f = scalar output.
            %
            % We would like to transform logpdf such that we can call it
            % like this:
            %   [f,g] = logpdf(x)
            % where:
            %   x = column vector.
            %   f = scalar output.
            %   g = numerical gradient as a column vector.
            
            z = @fun;
            function [f,g] = fun(x)
                % 1. Get the logpdf value.
                f = logpdf(x);
                
                % 2. Get the gradient via central differences (if needed).
                if nargout > 1
                    stepsize = eps(class(x))^(1/3);
                    g        = classreg.learning.fsutils.Solver.getGradient(logpdf,x,stepsize);
                end
            end
        end
    end
end