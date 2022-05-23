classdef InputValidator < classreg.learning.internal.DisallowVectorOps
%

%   Copyright 2016 The MathWorks, Inc.
    
    methods(Static)
        function stepsize = validateStepSize(stepsize)
            [isok,stepsize] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(stepsize,1);
            isok = isok && (stepsize > 0);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadStepSize'));
            end
        end
        
        function numsteps = validateNumSteps(numsteps)
            [isok,numsteps] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(numsteps,1);
            isok = isok && internal.stats.isIntegerVals(numsteps,1);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadNumSteps'));
            end
        end
        
        function massvector = validateMassVector(massvector,p)
            if isempty(massvector)
                massvector = ones(p,1);
            else
                [isok,massvector] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(massvector,p);
                isok = isok && all(massvector > 0);
                if ~isok
                    error(message('stats:mcmc:utils:InputValidator:BadMassVector',p));
                end
                massvector = massvector(:);
            end
        end
        
        function jittermethod = validateJitterMethod(jittermethod,allowedJitterMethods)
            jittermethod = internal.stats.getParamVal(jittermethod,allowedJitterMethods,'JitterMethod');
        end
        
        function historysize = validateHessianHistorySize(historysize)
            [isok,historysize] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(historysize,1);
            isok = isok && internal.stats.isIntegerVals(historysize,0);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadHessianHistorySize'));
            end
        end
        
        function samplertype = validateSamplerType(samplertype,allowedSamplerTypes)
            samplertype = internal.stats.getParamVal(samplertype,allowedSamplerTypes,'SamplerType');
        end
        
        function [gamma,isGammaAuto] = validateGamma(gamma,gammaAuto)
            if strcmpi(gamma,gammaAuto)
                isGammaAuto = true;
            else
                isGammaAuto = false;
                [isok,gamma] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(gamma,1);
                isok = isok && (gamma > 0);
                if ~isok
                    error(message('stats:mcmc:utils:InputValidator:BadGamma'));
                end
            end
        end
        
        function dojitterhistory = validateDoJitterHistory(dojitterhistory)
            [isok,dojitterhistory] = stats.mcmc.utils.InputValidator.isTrueFalseZeroOne(dojitterhistory);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadDoJitterHistory'));
            end
        end
        
        function dualgamma = validateDualAveragingGamma(dualgamma)
            [isok,dualgamma] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(dualgamma,1);
            isok = isok && (dualgamma > 0);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadDualAveragingGamma'));
            end
        end
        
        function dualt0 = validateDualAveragingT0(dualt0)
            [isok,dualt0] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(dualt0,1);
            isok = isok && (dualt0 > 0);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadDualAveragingT0'));
            end
        end
        
        function dualkappa = validateDualAveragingKappa(dualkappa)
            [isok,dualkappa] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(dualkappa,1);
            isok = isok && (dualkappa > 0.5) && (dualkappa <= 1);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadDualAveragingKappa'));
            end
        end
        
        function tuningmethod = validateStepSizeTuningMethod(tuningmethod,allowedTuningMethods)
            tuningmethod = internal.stats.getParamVal(tuningmethod,allowedTuningMethods,'StepSizeTuningMethod');
        end
        
        function tuningmethod = validateMassVectorTuningMethod(tuningmethod,allowedTuningMethods)
            tuningmethod = internal.stats.getParamVal(tuningmethod,allowedTuningMethods,'MassVectorTuningMethod');
        end
        
        function IterativeSamplingIterations = validateIterativeSamplingIterations(IterativeSamplingIterations)
            [isok,IterativeSamplingIterations] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(IterativeSamplingIterations,1);
            isok = isok && internal.stats.isIntegerVals(IterativeSamplingIterations,1);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadIterativeSamplingIterations'));
            end
        end
        
        function IterativeSamplingMinSamples = validateIterativeSamplingMinSamples(IterativeSamplingMinSamples)
            [isok,IterativeSamplingMinSamples] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(IterativeSamplingMinSamples,1);
            isok = isok && internal.stats.isIntegerVals(IterativeSamplingMinSamples,1);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadIterativeSamplingMinSamples'));
            end
        end
        
        function IterativeSamplingMaxSamples = validateIterativeSamplingMaxSamples(IterativeSamplingMaxSamples)
            [isok,IterativeSamplingMaxSamples] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(IterativeSamplingMaxSamples,1);
            isok = isok && internal.stats.isIntegerVals(IterativeSamplingMaxSamples,1);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadIterativeSamplingMaxSamples'));
            end
        end
        
        function IterativeSamplingInitialMassMultiplier = validateIterativeSamplingInitialMassMultiplier(IterativeSamplingInitialMassMultiplier)
            [isok,IterativeSamplingInitialMassMultiplier] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(IterativeSamplingInitialMassMultiplier,1);
            isok = isok && IterativeSamplingInitialMassMultiplier>0;
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadIterativeSamplingInitialMassMultiplier'));
            end
        end
        
        function IterativeSamplingNumStepsLimit = validateIterativeSamplingNumStepsLimit(IterativeSamplingNumStepsLimit)
            [isok,IterativeSamplingNumStepsLimit] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(IterativeSamplingNumStepsLimit,1);
            isok = isok && internal.stats.isIntegerVals(IterativeSamplingNumStepsLimit,1);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadIterativeSamplingNumStepsLimit'));
            end
        end
        
        function IterativeSamplingRegularizer = validateIterativeSamplingRegularizer(IterativeSamplingRegularizer)
            [isok,IterativeSamplingRegularizer] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(IterativeSamplingRegularizer,1);
            isok = isok && IterativeSamplingRegularizer>0;
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadIterativeSamplingRegularizer'));
            end
        end
        
        function checkgradient = validateCheckGradient(checkgradient)
            [isok,checkgradient] = stats.mcmc.utils.InputValidator.isTrueFalseZeroOne(checkgradient);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadCheckGradient'));
            end
        end
        
        function usenumericalgradient = validateUseNumericalGradient(usenumericalgradient)
            [isok,usenumericalgradient] = stats.mcmc.utils.InputValidator.isTrueFalseZeroOne(usenumericalgradient);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadUseNumericalGradient'));
            end
        end
    end
    
    methods(Static)
        function numsamples = validateNumSamples(numsamples)
            [isok,numsamples] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(numsamples,1);
            isok = isok && internal.stats.isIntegerVals(numsamples,0);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadNumSamples'));
            end
        end
        
        function thinSize = validateThinSize(thinSize)
            [isok,thinSize] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(thinSize,1);
            isok = isok && internal.stats.isIntegerVals(thinSize,1);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadThinSize'));
            end
        end
        
        function logpdf = validateLogPDF(logpdf)
            isok = isa(logpdf,'function_handle');
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadLogPDF'));
            end
        end
        
        function burnin = validateBurnin(burnin)
            [isok,burnin] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(burnin,1);
            isok = isok && internal.stats.isIntegerVals(burnin,0);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadBurnin'));
            end
        end
        
        function numtuningiter = validateNumStepSizeTuningIterations(numtuningiter)
            [isok,numtuningiter] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(numtuningiter,1);
            isok = isok && internal.stats.isIntegerVals(numtuningiter,0);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadNumStepSizeTuningIterations'));
            end
        end
        
        function targetaccratio = validateTargetAcceptanceRatio(targetaccratio)
            [isok,targetaccratio] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(targetaccratio,1);
            isok = isok && (targetaccratio > 0) && (targetaccratio < 1);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadTargetAcceptanceRatio'));
            end
        end
        
        function numstepslimit = validateNumStepsLimit(numstepslimit)
            [isok,numstepslimit] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(numstepslimit,1);
            isok = isok && internal.stats.isIntegerVals(numstepslimit,1);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadNumStepsLimit'));
            end
        end
        
        function donumstepswarn = validateDoNumStepsWarn(donumstepswarn)
            [isok,donumstepswarn] = stats.mcmc.utils.InputValidator.isTrueFalseZeroOne(donumstepswarn);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadDoNumStepsWarn'));
            end
        end
    end
    
    methods(Static)
        function start = validateStart(start,p)
            if ( nargin < 2 )
                p = [];
            end
            [isok,start] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(start,p);
            if ~isok
                if isempty(p)
                    error(message('stats:mcmc:utils:InputValidator:BadStart'));
                else
                    error(message('stats:mcmc:utils:InputValidator:BadStartLength',p));
                end
            end
        end
        
        function variableNames = validateVariableNames(variableNames,p)
            isok = isvector(variableNames) && iscellstr(variableNames) && (length(variableNames) == p);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadVariableNames',p));
            end
            variableNames = variableNames(:);
        end
        
        function haveGrad = validateHasGradient(haveGrad)
            [isok,haveGrad] = stats.mcmc.utils.InputValidator.isTrueFalseZeroOne(haveGrad);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadHasGradient'));
            end
        end
        
        function iterationLimit = validateIterationLimit(iterationLimit)
            [isok,iterationLimit] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(iterationLimit,1);
            isok = isok && internal.stats.isIntegerVals(iterationLimit,1);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadIterationLimit'));
            end
        end
        
        function gradientTolerance = validateGradientTolerance(gradientTolerance)
            [isok,gradientTolerance] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(gradientTolerance,1);
            isok = isok && (gradientTolerance > 0);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadGradientTolerance'));
            end
        end
        
        function stepTolerance = validateStepTolerance(stepTolerance)
            [isok,stepTolerance] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(stepTolerance,1);
            isok = isok && (stepTolerance > 0);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadStepTolerance'));
            end
        end
        
        function verbose = validateVerbose(verbose)
            [isok,verbose] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(verbose,1);
            isok = isok && internal.stats.isIntegerVals(verbose,0);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadVerbose'));
            end
        end
        
        function numprint = validateNumPrint(numprint)
            [isok,numprint] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(numprint,1);
            isok = isok && internal.stats.isIntegerVals(numprint,1);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadNumPrint'));
            end
        end
        
        function chains = validateChains(chains,p)
            if ~iscell(chains)
                chains = {chains};
            end
            
            okfun = @(x) isfloat(x) && isreal(x) && ismatrix(x) && (size(x,2) == p) && ~isempty(x);
            isok  = cellfun(okfun,chains);
            
            if ~all(isok)
                error(message('stats:mcmc:utils:InputValidator:BadChains',p));
            end
            
            chains = cellfun(@(z)z',chains,'UniformOutput',false);
        end
        
        function maxlag = validateMaxLag(maxlag)
            [isok,maxlag] = stats.mcmc.utils.InputValidator.isNumericRealVectorNoNaNInf(maxlag,1);
            isok = isok && internal.stats.isIntegerVals(maxlag,1);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadMaxLag'));
            end
        end
    end
    
    methods(Static)
        function [logpdf,start] = validateLogPDFAndStartNumericalGradient(logpdf,start)
            try
                lpdf = logpdf(start);
            catch causeException
                baseException = MException(message('stats:mcmc:utils:InputValidator:BadLogPDFStartCombination'));
                baseException = addCause(baseException,causeException);
                throw(baseException);
            end
            
            isok = isfloat(lpdf) && isreal(lpdf) && isscalar(lpdf);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadLogPDFFirstOutput'));
            end
        end
        
        function [logpdf,start,isgradrowvec] = validateLogPDFAndStartAnalyticalGradient(logpdf,start)
            try
                [lpdf,glpdf] = logpdf(start);
            catch causeException
                baseException = MException(message('stats:mcmc:utils:InputValidator:BadLogPDFStartCombination'));
                baseException = addCause(baseException,causeException);
                throw(baseException);
            end
            
            isok = isfloat(lpdf) && isreal(lpdf) && isscalar(lpdf);
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadLogPDFFirstOutput'));
            end
            
            [p,q]    = size(glpdf);
            lenstart = length(start);
            
            isgradcolumnvec = (p == lenstart) && (q == 1);
            isgradrowvec    = (q == lenstart) && (p == 1);
            isgradvec       = isgradcolumnvec || isgradrowvec;
            
            isok = isfloat(glpdf) && isreal(glpdf) && isgradvec;
            if ~isok
                error(message('stats:mcmc:utils:InputValidator:BadLogPDFSecondOutput',lenstart));
            end
        end
    end
    
    methods(Static)
        function [isok,x] = isNumericRealVectorNoNaNInf(x,N)
        % INPUTS:
        %   x = a potential numeric, real vector.
        %   N = length of x or [] if length of x is not known.
        % OUTPUT:
        %   isok = true if x is a numeric real vector of length N. NaN and 
        %          Inf values are not allowed. If N is empty x can be of 
        %          any length. 
        %      x = validated value of x as a column vector if isok is true.
        %
        % NOTE: If x contains integer values such as int8, uint8 etc. then
        %       x is cast to a double.
            
            isok = isnumeric(x) && isreal(x) && isvector(x) && ~any(isnan(x)) && ~any(isinf(x));            
            if ( isempty(N) )
                % x can be of any length.
            else
                % x must be of length N.
                isok = isok && (length(x) == N);
            end
            if ( isok && (size(x,1) == 1) )
                % Make into column vector.
                x = x';
            end
            if ( isok && isinteger(x) )
                x = cast(x,'double');
            end
        end
        
        function [isok,x] = isTrueFalseZeroOne(x)
        % INPUTS:
        %   x = a potential 0/1 or true/false value.
        % OUTPUTS:
        %   isok = true if x is valid.
        %      x = validated value of x as a logical if isok is true.
            
            if ( isscalar(x) && islogical(x) )
                isok = true;
                return;
            end
            
            isint = internal.stats.isScalarInt(x);
            if ( isint )
                if  ( x == 1 )
                    isok = true;
                    x    = true;
                elseif ( x == 0 )
                    isok = true;
                    x    = false;
                else
                    isok = false;
                end
            else
                isok = false;
            end
        end        
    end
    
end