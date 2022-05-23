classdef HamiltonianImpl < stats.mcmc.impl.MCMCImpl
    
    %% Constants
    properties(Constant)
        JitterMethodNoneCode           = 0;
        JitterMethodJitterBothCode     = 1;
        JitterMethodJitterNumStepsCode = 2;
    end
        
    %% Constructor
    methods
        function this = HamiltonianImpl(logpdf,haveGrad,start,variableNames,name,samplerParameters)
            % 1. Set properties.
            this.LogPDF            = logpdf;
            this.HasGradient       = haveGrad;
            this.StartPoint        = start;
            this.Name              = name;
            this.VariableNames     = variableNames;
            this.SamplerParameters = samplerParameters;
            
            % 2. Check gradient (if needed).
            if samplerParameters.CheckGradient
                stats.mcmc.impl.HamiltonianImpl.validateGradient(logpdf,start);
            end
        end
    end
    
    %% Drawing samples
    methods
        function [xsmpl,state,accratio,stepSizeTuningInfo] = drawSamples(this,chainInfo,samplerParams)
            % 1. Get chain info.
            start          = chainInfo.StartPoint;
            numtuningiter  = chainInfo.NumStepSizeTuningIterations;
            targetaccratio = chainInfo.TargetAcceptanceRatio;
            burnin         = chainInfo.Burnin;
            numsamples     = chainInfo.NumSamples;
            verbose        = chainInfo.VerbosityLevel;
            numprint       = chainInfo.NumPrint;
            thin           = chainInfo.ThinSize;
            numstepslimit  = chainInfo.NumStepsLimit;
            donumstepswarn = chainInfo.DoNumStepsWarn;
            
            % 2. Get sampling parameters.
            stepsize     = samplerParams.StepSize;
            numsteps     = samplerParams.NumSteps;
            massvec      = samplerParams.MassVector;
            jittermethod = samplerParams.JitterMethod;
            dualgamma    = samplerParams.DualAveragingGamma;
            dualt0       = samplerParams.DualAveragingT0;
            dualkappa    = samplerParams.DualAveragingKappa;
            
            % 3. Are we going to tune the step size?
            if strcmpi(samplerParams.StepSizeTuningMethod,stats.mcmc.params.HamiltonianSamplingParameters.StepSizeTuningMethodDualAveraging)
                dotuning = true;
            else
                dotuning = false;
            end
            
            % 4. Do sampling and optionally stepsize tuning.
            if ( verbose > 0 )
                fprintf('\n');
            end
            [xsmpl,accratio,stepSizeTuningInfo,state] = stats.mcmc.impl.HamiltonianImpl.doHMC(start,numsamples,this.LogPDF,burnin,stepsize,numsteps,massvec,...
                                                                                jittermethod,numtuningiter,targetaccratio,dualgamma,dualt0,dualkappa,numstepslimit,...
                                                                                donumstepswarn,dotuning,verbose,numprint,thin);
        end
    end
    
    %% Tuning parameters
    methods
        function [this,tuningInfo,state] = tuneSampler(this,chainInfo,samplerParams)
            % 1. Tune MassVector.
            [samplerParams, massVectorTuningInfo] = tuneMassVector(this,chainInfo,samplerParams);
            
            % 2. Tune StepSize and NumSteps.
            [samplerParams,stepSizeTuningInfo,state] = tuneStepSize(this,chainInfo,samplerParams);
            
            % 3. Store tuning results in tuningInfo
            tuningInfo.MassVector           = samplerParams.MassVector;
            tuningInfo.StepSize             = samplerParams.StepSize;
            tuningInfo.NumSteps             = samplerParams.NumSteps;
            tuningInfo.MassVectorTuningInfo = massVectorTuningInfo;
            tuningInfo.StepSizeTuningInfo   = stepSizeTuningInfo;
            
            % 4. Update the sampling parameters in this.
            this.SamplerParameters = samplerParams;
            
            % 5. Save final state after tuning into this as a StartPoint.
            this.StartPoint = state;
        end
        
        function [samplerParams, massVectorTuningInfo] = tuneMassVector(this,chainInfo,samplerParams)
            if strcmpi(samplerParams.MassVectorTuningMethod,stats.mcmc.params.HamiltonianSamplingParameters.MassVectorTuningMethodHessian)
                [samplerParams, massVectorTuningInfo] = tuneMassVectorHessian(this,chainInfo,samplerParams);
            elseif strcmpi(samplerParams.MassVectorTuningMethod,stats.mcmc.params.HamiltonianSamplingParameters.MassVectorTuningMethodIterativeSampling)
                [samplerParams, massVectorTuningInfo] = tuneMassVectorIterativeSampling(this,chainInfo,samplerParams);
            else
                massVectorTuningInfo = [];
            end
        end
        
        function [samplerParams, massVectorTuningInfo] = tuneMassVectorHessian(this,chainInfo,samplerParams)
            % 1. Display mass vector tuning message.
            if ( chainInfo.VerbosityLevel > 0 )
                fprintf([getString(message('stats:mcmc:impl:HamiltonianImpl:TuningMassVector', samplerParams.MassVectorTuningMethod)),' \n']);
            end
            
            % 2. Set mass vector equal to the negative diagonal Hessian
            % of LogPDF. Warn about negative values in mass vector or
            % values close to 0.
            tol            = 1e-3;
            negdiaghessian = -1*this.getDiagonalHessian(this.LogPDF,chainInfo.StartPoint);
            
            if any(negdiaghessian < 0)
                tolstr = num2str(tol);
                warning(message('stats:mcmc:impl:HamiltonianImpl:TunedMassVectorHasNegativeElements',tolstr));
            elseif any(negdiaghessian < tol)
                tolstr = num2str(tol);
                warning(message('stats:mcmc:impl:HamiltonianImpl:TunedMassVectorHasSmallElements',tolstr,tolstr));
            end
            
            % 3. Mass vector must be >= tol.
            massvec                  = max(tol,negdiaghessian);
            
            % 4. Set return values
            samplerParams.MassVector = massvec;
            massVectorTuningInfo = struct('MassVector', massvec,...
                                          'NegativeDiagonalHessian', negdiaghessian, ...
                                          'HessianPoint', chainInfo.StartPoint);
        end
        
        function [samplerParams, massVectorTuningInfo] = tuneMassVectorIterativeSampling(this,chainInfo,samplerParams)
            % 1. Display mass vector tuning message.
            if ( chainInfo.VerbosityLevel > 0 )
                fprintf([getString(message('stats:mcmc:impl:HamiltonianImpl:TuningMassVector', samplerParams.MassVectorTuningMethod)),' \n']);
            end
            
            % 2. Tune using iterative sampling
            massvecProfile                              = zeros(samplerParams.NumParameters, samplerParams.IterativeSamplingIterations);
            numsamplesProfile                           = zeros(samplerParams.IterativeSamplingIterations, 1);
            samplerParamsNew                            = samplerParams;
            samplerParamsNew.StepSizeTuningMethod       = samplerParams.StepSizeTuningMethodDualAveraging;
            samplerParamsNew.MassVector                 = samplerParams.MassVector * samplerParams.IterativeSamplingInitialMassMultiplier;
            chainInfo.NumSamples                        = samplerParams.IterativeSamplingMinSamples;
            chainInfo.Burnin                            = 0;
            chainInfo.NumStepsLimit                     = samplerParams.IterativeSamplingNumStepsLimit;
            chainInfo.DoNumStepsWarn                    = false;
            IterativeSamplingVerbosityLevel             = chainInfo.VerbosityLevel;
            chainInfo.VerbosityLevel                    = double(chainInfo.VerbosityLevel > 1);
            for tuneIter = 1:samplerParams.IterativeSamplingIterations
                [xsmpl, finalState] = drawSamples(this, chainInfo, samplerParamsNew);
                if IterativeSamplingVerbosityLevel > 0
                    fprintf([getString(message('stats:mcmc:impl:HamiltonianImpl:TuningMassVectorIteration', num2str(tuneIter), num2str(samplerParams.IterativeSamplingIterations))),' \n']);
                end
                newmassvec                  = 1./(var(xsmpl,0,2) + samplerParams.IterativeSamplingRegularizer);
                samplerParamsNew.MassVector = newmassvec;
                massvecProfile(:,tuneIter)  = newmassvec;
                numsamplesProfile(tuneIter) = chainInfo.NumSamples;
                chainInfo.NumSamples        = min(2*chainInfo.NumSamples, samplerParams.IterativeSamplingMaxSamples);
                chainInfo.StartPoint        = finalState;
            end
                        
            % 3. Set return values
            samplerParams.MassVector = samplerParamsNew.MassVector;
            massVectorTuningInfo     = struct('MassVector', samplerParamsNew.MassVector,...
                                              'IterativeSamplingMassVectorProfile', massvecProfile, ...
                                              'IterativeSamplingNumSamples', numsamplesProfile);
        end
        
        function [samplerParamsDA,stepSizeTuningInfo,state] = tuneStepSize(this,chainInfoDA,samplerParamsDA)
            % 1. Set Burnin and NumSamples to 0 since we just want tuning.
            chainInfoDA.Burnin     = 0;
            chainInfoDA.NumSamples = 0;
            
            % 2. Display step size tuning message.
            if ( chainInfoDA.VerbosityLevel > 0 )
                fprintf(['\n', getString(message('stats:mcmc:impl:HamiltonianImpl:TuningStepSize', samplerParamsDA.StepSizeTuningMethod, num2str(chainInfoDA.TargetAcceptanceRatio))),'\n']);
            end
            
            % 3. Do tuning via drawSamples.
            [~,state,accratio,stepSizeTuningInfo] = drawSamples(this,chainInfoDA,samplerParamsDA);
            stepSizeTuningInfo.AcceptanceRatio    = accratio;
            
            % 4. Set optimal StepSize and NumSteps in samplerParamsDA.
            samplerParamsDA.StepSize = stepSizeTuningInfo.StepSize;
            samplerParamsDA.NumSteps = stepSizeTuningInfo.NumSteps;
        end
    end
    
    %% Convergence diagnostics
    methods
        function stats = diagnostics(this,chains,maxlag)
            % 1. Set historysize to 0 for HMC.
            historysize = 0;
            
            % 2. Call superclass method with correct historysize.
            stats = diagnostics@stats.mcmc.impl.MCMCImpl(this,chains,maxlag,historysize);
        end
    end    
  
    %% Gradient check
    methods(Static)
        function validateGradient(logpdf,start)
            [isok,err] = stats.mcmc.utils.checkGradient(logpdf,start);
            if ~isok
                warning(message('stats:mcmc:impl:HamiltonianImpl:BadGradient',num2str(err)));
            end
        end
    end
    
    %% Jittering HMC parameters
    methods(Static)
        function jittercode = getJitterCode(jittermethod)
            if strcmpi(jittermethod,stats.mcmc.params.HamiltonianSamplingParameters.JitterMethodNone)
                
                jittercode = stats.mcmc.impl.HamiltonianImpl.JitterMethodNoneCode;
                
            elseif strcmpi(jittermethod,stats.mcmc.params.HamiltonianSamplingParameters.JitterMethodJitterBoth)
                
                jittercode = stats.mcmc.impl.HamiltonianImpl.JitterMethodJitterBothCode;
                
            elseif strcmpi(jittermethod,stats.mcmc.params.HamiltonianSamplingParameters.JitterMethodJitterNumSteps)
                
                jittercode = stats.mcmc.impl.HamiltonianImpl.JitterMethodJitterNumStepsCode;
            end
        end
        
        function [stepsize,numsteps] = hmcParameters(maxstepsize,maxnumsteps,jittercode)
            switch jittercode
                case stats.mcmc.impl.HamiltonianImpl.JitterMethodNoneCode
                    
                    stepsize = maxstepsize;
                    numsteps = maxnumsteps;
                    
                case stats.mcmc.impl.HamiltonianImpl.JitterMethodJitterBothCode
                    
                    if rand < 0.5
                        stepsize = unifrnd(0,maxstepsize,1,1);
                    else
                        stepsize = 0.9*maxstepsize;
                    end
                    numsteps = randi([1,maxnumsteps],1,1);
                    
                case stats.mcmc.impl.HamiltonianImpl.JitterMethodJitterNumStepsCode
                    
                    stepsize = maxstepsize;
                    numsteps = randi([1,maxnumsteps],1,1);
            end
        end        
    end
    
    %% LeapFrog iterations
    methods(Static)
% Here's another way to do LeapFrog iterations:
%
%         function [currx,lpdfcurrx,glpdfcurrx,currz] = doLeapFrog(logpdf,currx,lpdfcurrx,glpdfcurrx,currz,stepsize,numsteps,massvector)
%             % Half momentum step.
%             dUdx  = -glpdfcurrx;
%             currz = currz - (stepsize/2)*dUdx;
%             
%             % Full position and full momentum steps.
%             for r = 1:(numsteps-1)
%                 dKdz  = currz./massvector;
%                 currx = currx + stepsize*dKdz;
%                 
%                 [~,glpdfcurrx] = logpdf(currx);
%                 dUdx           = -glpdfcurrx;
%                 currz          = currz - stepsize*dUdx;
%             end
%             
%             % Full position and half momentum step.
%             dKdz  = currz./massvector;
%             currx = currx + stepsize*dKdz;
%             
%             [lpdfcurrx,glpdfcurrx] = logpdf(currx);
%             dUdx                   = -glpdfcurrx;
%             currz                  = currz - (stepsize/2)*dUdx;
%         end

        function [currx,lpdfcurrx,glpdfcurrx,currz] = doLeapFrog(logpdf,currx,lpdfcurrx,glpdfcurrx,currz,stepsize,numsteps,massvector)
            for r = 1:numsteps
                % x(t), z(t) -> x(t), z(t + epsilon/2)
                dUdx  = -glpdfcurrx;
                currz = currz - (stepsize/2)*dUdx;
                
                % x(t), z(t + epsilon/2) -> x(t + epsilon), z(t + epsilon/2)
                dKdz  = currz./massvector;
                currx = currx + stepsize*dKdz;
                
                % x(t + epsilon), z(t + epsilon/2) -> x(t + epsilon), z(t + epsilon)
                [lpdfcurrx,glpdfcurrx] = logpdf(currx);
                dUdx                   = -glpdfcurrx;
                currz                  = currz - (stepsize/2)*dUdx;
            end
        end        
    end
    
    %% StepSize tuning by dual-averaging
    methods(Static)
        function stepsize = chooseInitialStepSize(logpdf,x,lpdfx,glpdfx,massvector)
            % 1. Find acceptance probability of a Langevin proposal using
            % the specified mass vector and an initial stepsize = 1.
            p        = length(x);
            sigma    = sqrt(massvector);
            stepsize = 1;
            numsteps = 1;
            
            zstar           = sigma .* randn(p,1);
            [~,~,~,~,aprob] = stats.mcmc.impl.HamiltonianImpl.doHMCProposal(logpdf,x,lpdfx,glpdfx,zstar,stepsize,numsteps,massvector);
            
            % 2. If acceptance probability is > 0.5 already, we want to
            % increase the stepsize and if acceptance probability is < 0.5
            % we want to decrease the stepsize. This process doubles or
            % halves the stepsize until the acceptance probability of the
            % Langevin proposal crosses 0.5.
            a = 2*(aprob > 0.5) - 1;
            
            while (aprob^a > 0.5^a)
                stepsize        = (2^a)*stepsize;
                [~,~,~,~,aprob] = stats.mcmc.impl.HamiltonianImpl.doHMCProposal(logpdf,x,lpdfx,glpdfx,zstar,stepsize,numsteps,massvector);
            end
        end
        
        function [stepsizeBar,numstepsBar,x,lpdfx,glpdfx,nacc,stepSizeBarProfile] = tuneHMCDualAveraging(x0,logpdf,numtune,maxstepsize,maxnumsteps,...
                                                                                                         massvector,jittermethod,delta,gamma,t0,kappa,...
                                                                                                         numstepslimit,donumstepswarn,verbose,numprint)
            % 1. What is our adaptation method for stepsize and numsteps?
            % NOTE: This function is currently being called only with
            % JitterMethod = 'none'. It may be the only sensible way to call
            % this function. In that case, we could remove the jittermethod
            % input argument from this function.
            jittercode = stats.mcmc.impl.HamiltonianImpl.getJitterCode(jittermethod);
            
            % 2. sigma is the p-by-1 standard deviation vector for sampling
            % momentum vector.
            p     = length(x0);
            sigma = sqrt(massvector);
            
            % 3. Initial logpdf and it's gradient.
            %    nacc = number of steps accepted so far.
            %    ndiv = number of divergent steps so far.
            %     x   = current position of the chain.
            %   lpdfx = log pdf at x.
            %  glpdfx = gradient of log pdf at x.
            nacc           = 0;
            ndiv           = 0;
            x              = x0;
            [lpdfx,glpdfx] = logpdf(x);
            
            % 4. Get initial stepsize.
            stepsize = stats.mcmc.impl.HamiltonianImpl.chooseInitialStepSize(logpdf,x,lpdfx,glpdfx,massvector);
            
            if ~isfinite(stepsize)
                warning(message('stats:mcmc:impl:HamiltonianImpl:BadInitialStepSizeForDualAveraging'));
                stepsize = 1;
            end
            
            if ( verbose > 0 )
                fprintf([getString(message('stats:mcmc:impl:HamiltonianImpl:InitialStepSizeForDualAveraging',num2str(stepsize))),'\n']);
            end            
            
            % 5. Randomly draw a simulation length as per jitter method.
            % Then keeping the same simulation length, compute numsteps for
            % stepsize. Finally, clip numSteps in [1,numstepslimit]. This
            % may shorten the simulation length.
            [stepsize0,numsteps0] = stats.mcmc.impl.HamiltonianImpl.hmcParameters(maxstepsize,maxnumsteps,jittercode);
            lambda0               = stepsize0*numsteps0;
            [numsteps, clipped]   = stats.mcmc.impl.HamiltonianImpl.clipNumSteps(round(lambda0/stepsize), numstepslimit, donumstepswarn);
            if clipped
                donumstepswarn = false;
            end
            
            % 6. Initialize dual averaging quantities. Note that we
            % optimize (delta - accratio) to be close to 0 as a function of
            % log(stepsize).
            stepsizeBar        = 1;
            gBar               = 0;
            mu                 = log(10*stepsize);            
            stepSizeBarProfile = zeros(numtune,1);
            
            % 7. ncalls is the number of calls to display utility.
            if ( verbose > 0 )
                ncalls = 0;
            end
            
            % 8. Dual averaging iterations.
            for i = 1:numtune
                zstar = sigma .* randn(p,1);
                
                [x,lpdfx,glpdfx,accepted,aprob,isdivergent] = stats.mcmc.impl.HamiltonianImpl.doHMCOneStep(logpdf,x,lpdfx,glpdfx,zstar,stepsize,...
                                                                                                           numsteps,massvector);                                
                
                if accepted
                    nacc = nacc + 1;
                end
                
                if isdivergent
                    ndiv = ndiv + 1;
                end
                
                g     = delta - aprob;
                eta_g = 1/(t0 + i);
                gBar  = eta_g*g + (1-eta_g)*gBar;
                
                stepsize = exp(mu - (sqrt(i)/gamma)*gBar);
                
                eta_stepsize   = i^(-kappa);
                logStepSizeBar = eta_stepsize*log(stepsize) + (1-eta_stepsize)*log(stepsizeBar);
                stepsizeBar    = exp(logStepSizeBar);
                                                
                [stepsize0,numsteps0] = stats.mcmc.impl.HamiltonianImpl.hmcParameters(maxstepsize,maxnumsteps,jittercode);
                lambda0               = stepsize0*numsteps0;
                [numsteps, clipped]   = stats.mcmc.impl.HamiltonianImpl.clipNumSteps(round(lambda0/stepsize), numstepslimit, donumstepswarn);
                if clipped
                    donumstepswarn = false;
                end
                
                stepSizeBarProfile(i) = stepsizeBar;
                
                if ( verbose > 0 && rem(i,numprint) == 0 )
                    if ( rem(ncalls,20) == 0 )
                        doprintheader = true;
                    else
                        doprintheader = false;
                    end
                    ncalls = ncalls + 1;
                    
                    numstepsBar = stats.mcmc.impl.HamiltonianImpl.clipNumSteps(round(lambda0/stepsizeBar), numstepslimit, false);
                    stats.mcmc.impl.HamiltonianImpl.displaySamplerInfo(i,lpdfx,stepsizeBar,numstepsBar,nacc/i,ndiv,doprintheader);
                end
            end

            % 9. Set step size and number of steps based on stepsizeBar.
            lambda0     = maxstepsize*maxnumsteps;
            numstepsBar = stats.mcmc.impl.HamiltonianImpl.clipNumSteps(round(lambda0/stepsizeBar), numstepslimit, donumstepswarn);
        end
        
        function [numsteps, clipped] = clipNumSteps(numsteps, numstepslimit, dowarn)
            if numsteps > numstepslimit
                clipped = true;
                if dowarn
                    warning(message('stats:mcmc:impl:HamiltonianImpl:LargeNumStepsInDualAveraging', num2str(numsteps)));
                end
            else
                clipped = false;
            end
            numsteps = max(1,min(numstepslimit,numsteps));
        end
    end
    
    %% Convergence summary
    methods(Static)
        function displaySamplerInfo(iter,lpdf,stepsize,numsteps,accratio,ndiv,doprintheader)
        %displaySamplerInfo - Helper function to display iteration info.
        %   displaySamplerInfo(iter,lpdf,stepsize,numsteps,accratio,ndiv,doprintheader)
        %   accepts several inputs (described below) and prints out a
        %   diagnostic summary of progress made by the sampler.
        %
        %   INPUTS:
        %
        %   iter          = iteration number.
        %   lpdf          = log PDF.
        %   stepsize      = step size for HMC.
        %   numsteps      = number of LeapFrog steps for HMC.
        %   accratio      = currently realized acceptance ratio.
        %   ndiv          = number of divergent steps so far.
        %   doprintheader = true to print header and false otherwise.
        %
        %   We will display convergence info like this:
        %
        % |==================================================================================|
        % |   ITER   |    LOG PDF    |  STEP SIZE  |  NUM STEPS  |  ACC RATIO  |  DIVERGENT  |
        % |==================================================================================|
        % |      100 |  1.915580e+01 |   8.877e-02 |          57 |   6.200e-01 |           0 |
        % |      200 |  1.531223e+01 |   1.044e-01 |          49 |   6.450e-01 |           0 |
        % |      300 |  1.505835e+00 |   1.031e-01 |          49 |   6.300e-01 |           0 |
        % |      400 |  1.701939e+01 |   1.064e-01 |          48 |   6.300e-01 |           0 |
        % |      500 |  1.970405e+00 |   1.093e-01 |          47 |   6.420e-01 |           0 |
        % |      600 |  2.149306e+01 |   1.076e-01 |          47 |   6.417e-01 |           0 |
        % |      700 |  7.537130e+00 |   1.128e-01 |          45 |   6.386e-01 |           0 |
        % |      800 |  1.395544e+01 |   1.101e-01 |          46 |   6.425e-01 |           0 |
        % |      900 |  1.165858e+01 |   1.113e-01 |          46 |   6.467e-01 |           0 |
        % |     1000 |  1.275114e+01 |   1.129e-01 |          45 |   6.480e-01 |           0 |

        
            % 1. Display header if required.
            if doprintheader
                fprintf('\n');
                fprintf('|==================================================================================|\n');
                fprintf('|   ITER   |    LOG PDF    |  STEP SIZE  |  NUM STEPS  |  ACC RATIO  |  DIVERGENT  |\n');
                fprintf('|==================================================================================|\n');
            end

            % 2. Display iteration wise info.
            fprintf('|%9d |%14.6e |%12.3e |%12d |%12.3e |%12d |\n', iter, lpdf, stepsize, numsteps, accratio, ndiv);
        end
    end
    
    %% HMC sampling
    methods(Static)                
        function [xsmpl,accratio,stepSizeTuningInfo,state] = doHMC(x0,numsamples,logpdf,burnin,maxstepsize,maxnumsteps,massvector,jittermethod,...
                                                           numtuningiter,targetaccratio,dualgamma,dualt0,dualkappa,numstepslimit,donumstepswarn,dotuning,verbose,numprint,thin)
            % 1. What is our adaptation method for stepsize and numsteps?
            jittercode = stats.mcmc.impl.HamiltonianImpl.getJitterCode(jittermethod);
            
            % 2. Prepare to start sampling. sigma is the p-by-1 standard
            % deviation vector for sampling momentum vector.
            p     = length(x0);
            sigma = sqrt(massvector);
            xsmpl = nan(p,numsamples);
            
            % 3. Initialize variables and tune parameters by dual-averaging
            % if requested.
            %
            %   nincl  = number of samples included for output so far beyond the burnin period.
            %   ndiv   = number of divergent steps so far.
            %   ncalls = the number of calls to display utility if (verbose > 0).
            %   nacc   = number of steps accepted so far.
            %     x    = current position of the chain.
            %   lpdfx  = log pdf at x.
            %  glpdfx  = gradient of log pdf at x.
            nincl = 0;
            ndiv  = 0;
            
            if ( verbose > 0 )
                ncalls = 0;
            end
            
            if ( dotuning && (numtuningiter > 0) )
               [maxstepsize,maxnumsteps,x,lpdfx,glpdfx,nacc,stepsizeprofile] = stats.mcmc.impl.HamiltonianImpl.tuneHMCDualAveraging(x0,logpdf,...
                                                                                    numtuningiter,maxstepsize,maxnumsteps,massvector,...
                                                                                    stats.mcmc.params.HamiltonianSamplingParameters.JitterMethodNone,...
                                                                                    targetaccratio,dualgamma,dualt0,dualkappa,numstepslimit,donumstepswarn,verbose,numprint);
               
               stepSizeTuningInfo.StepSize        = maxstepsize;
               stepSizeTuningInfo.NumSteps        = maxnumsteps;
               stepSizeTuningInfo.StepSizeProfile = stepsizeprofile;
            else
                nacc           = 0;
                x              = x0;
                [lpdfx,glpdfx] = logpdf(x);
                
                stepSizeTuningInfo.StepSize        = maxstepsize;
                stepSizeTuningInfo.NumSteps        = maxnumsteps;
                stepSizeTuningInfo.StepSizeProfile = [];
            end
            
            % 4. Freeze tuned parameters and begin the sampling phase -
            % this includes the burnin period.
            maxiter = burnin + numsamples*thin;
            
            for i = 1:maxiter
                % 4.1 Get zstar.
                zstar = sigma .* randn(p,1);
                
                % 4.2 Randomly select stepsize and numsteps based on chosen
                % jittering method.
                [stepsize,numsteps] = stats.mcmc.impl.HamiltonianImpl.hmcParameters(maxstepsize,maxnumsteps,jittercode);
                
                % 4.3 Take one HMC step.
                [x,lpdfx,glpdfx,accepted,~,isdivergent] = stats.mcmc.impl.HamiltonianImpl.doHMCOneStep(logpdf,x,lpdfx,glpdfx,zstar,stepsize,numsteps,massvector);
                
                % 4.4 Update acceptance counter.
                if accepted
                    nacc = nacc + 1;
                end
                
                if isdivergent
                    ndiv = ndiv + 1;
                end
                
                % 4.5 Start saving states if burnin period is over.
                % Account for thinning while saving states.
                if ( i > burnin && rem(i-burnin,thin) == 0 )
                    nincl = nincl + 1;
                    
                    xsmpl(:,nincl) = x;
                end
                
                % 4.6 Print diagnostic info.
                if ( verbose > 0 && rem(i,numprint) == 0 )
                    if ( rem(ncalls,20) == 0 )
                        doprintheader = true;
                    else
                        doprintheader = false;
                    end
                    ncalls = ncalls + 1;
                    
                    stats.mcmc.impl.HamiltonianImpl.displaySamplerInfo(i,lpdfx,stepsize,numsteps,nacc/(i+numtuningiter),ndiv,doprintheader);
                end
            end
            
            % 5. Acceptance ratio including the burnin and tuning period.
            accratio = nacc/(maxiter + numtuningiter);
            
            % 6. Save final state of the chain.
            state = x;
        end
        
        function [currx,lpdfcurrx,glpdfcurrx,currz,aprob,isdivergent] = doHMCProposal(logpdf,x,lpdfx,glpdfx,zstar,stepsize,numsteps,massvector)            
            % 1. Do LeapFrog iterations. Proposal for (x,z) is (xp,zp) = (currx,-currz).
            [currx,lpdfcurrx,glpdfcurrx,currz] = stats.mcmc.impl.HamiltonianImpl.doLeapFrog(logpdf,x,lpdfx,glpdfx,zstar,stepsize,numsteps,massvector);
            
            % 2. Compute acceptance probability and mark divergence.
            % Expzp   = E(xp,zp) = E(currx,-currz)
            % Exzstar = E(x,zstar)
            %
            % Mark the proposal as divergent if expterm below turns out to
            % be NaN. Results from '* are not exactly reproducible across
            % machines. So the next 2 lines are rewritten differently
            % without using '*.
            % Expzp   = -lpdfcurrx + 0.5*currz'*(currz./massvector);
            % Exzstar = -lpdfx     + 0.5*zstar'*(zstar./massvector);
            Expzp   = -lpdfcurrx + 0.5*sum(currz.*(currz./massvector));
            Exzstar = -lpdfx     + 0.5*sum(zstar.*(zstar./massvector));
            expterm = exp(Exzstar - Expzp);
            
            if isnan(expterm)
                aprob       = 0;
                isdivergent = true;
            else
                aprob       = min(1,expterm);
                isdivergent = false;
            end
        end
        
        function [x,lpdfx,glpdfx,accepted,aprob,isdivergent] = doHMCOneStep(logpdf,x,lpdfx,glpdfx,zstar,stepsize,numsteps,massvector)
            % 1. Run numsteps of Hamiltonian dynamics with step size stepsize.
            currx      = x;
            lpdfcurrx  = lpdfx;
            glpdfcurrx = glpdfx;
            currz      = zstar;
            
            % 2. Do LeapFrog iterations. Proposal is (xp,zp)=(currx,-currz)
            % but we don't care about the z proposal.
            [currx,lpdfcurrx,glpdfcurrx,~,aprob,isdivergent] = stats.mcmc.impl.HamiltonianImpl.doHMCProposal(logpdf,currx,lpdfcurrx,glpdfcurrx,currz,...
                                                                                                             stepsize,numsteps,massvector);
            
            % 3. Update state of the chain or remain at the same point.
            U = rand;
            if U <= aprob
                % Accept.
                x      = currx;
                lpdfx  = lpdfcurrx;
                glpdfx = glpdfcurrx;
                
                % Mark the transition as accepted.
                accepted = true;
            else
                % Mark the transition as not accepted.
                accepted = false;
            end
        end                        
    end
    
    %% Numerical approximation to diagonal Hessian
    methods(Static)
        function Hdiag = getDiagonalHessian(fun,theta)
        % Get the diagonal Hessian of the function fun at theta.
        %
        % Hdiag(i) = (fun(theta + 2*h*ei) - 2*fun(theta) + fun(theta - 2*h*ei))/(4*h^2)
        %
        % where ei is a vector of all zeros except a 1 at position i and h
        % is a small step size.
            
            % 1. Set step size.
            step = eps^(1/4);
            
            % 2. Initialize output.
            p     = length(theta);
            Hdiag = zeros(p,1);
            
            % 3. Evaluate fun(theta) once.
            funtheta = fun(theta);
            
            % 4. Get the diagonal elements of Hdiag.
            for i = 1:p
                % Use central differences.
                theta2    = theta;
                theta2(i) = theta2(i) + 2*step;
                
                theta1    = theta;
                theta1(i) = theta1(i) - 2*step;
                
                Hdiag(i)  = (fun(theta2) + fun(theta1) - 2*funtheta)/(4*step*step);
            end
        end
    end
    
end






