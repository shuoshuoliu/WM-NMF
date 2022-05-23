classdef (Abstract) MCMCImpl < classreg.learning.internal.DisallowVectorOps
%

%   Copyright 2016 The MathWorks, Inc.
    
    properties
        LogPDF;                 % Function handle to compute LogPDF.
        StartPoint;             % Start point for the Markov chain.
        Name;                   % Name of the sampler.
        SamplerParameters;      % Sampling parameters.
        HasGradient;            % True if LogPDF can compute gradient.
        VariableNames;          % Cell array of length P-by-1 containing variable names where P is the length of StartPoint.
    end
    
    methods(Abstract)
        drawSamples(this,varargin);        
        tuneSampler(this,varargin);
    end
    
    % This is a default implementation based on LogPDF property. For
    % stochastic gradient samplers, this method should be overridden to
    % compute the negative log PDF based on log likelihood and prior.
    methods
        function z = makeNegativeLogPDFForMAP(this)
            logPDF = this.LogPDF;
            z = @negLogPDF;
            function [fun,grad] = negLogPDF(x)
                if ( nargout < 2 )
                    fun = logPDF(x);
                    fun = -fun;
                else
                    [fun,grad] = logPDF(x);
                    
                    fun  = -fun;
                    grad = -grad;
                end
            end
        end
    end
    
    methods       
        function [xHat,fitInfo] = estimateMAP(this,x0,haveGrad,iterationLimit,verbose,gradTol,stepTol)
            % 1. Get negative log posterior.
            fun = makeNegativeLogPDFForMAP(this);

            % 2. Create a solver object.
            solver                   = classreg.learning.fsutils.Solver(1);
            solver.HaveGradient      = haveGrad;
            solver.Verbose           = verbose;
            solver.IterationLimit    = iterationLimit;
            solver.GradientTolerance = gradTol;
            solver.StepTolerance     = stepTol;

            history      = struct();
            history.fval = [];
            history.iter = [];
            
            % 3. Set up output function.
            function stop = outfun(~,optimValues,~)
                history.iter = [history.iter; optimValues.iteration];
                history.fval = [history.fval; optimValues.fval];
                stop         = false;
            end
            
            % 4. Do minimization.
            results = doMinimization(solver,fun,x0,1,'OutputFcn',@outfun);
            
            % 5. Set xHat.
            xHat = results.xHat;
            
            % 6. Set fitInfo.
            fitInfo           = struct();
            fitInfo.Iteration = history.iter;
            fitInfo.Objective = history.fval;
            fitInfo.Gradient  = results.gHat;
        end

        function statstbl = diagnostics(this,chains,maxlag,historysize)            
            % 1. Set default historysize. Note that historysize is > 0 only
            % for those samplers that generate an order K Markov chain with
            % K > 1 - for example, the quasi-Newton HMC sampler.
            if ( nargin < 4 )
                historysize = 0;
            end
            
            % 2. How many variable? How many chains?
            P = length(this.StartPoint);
            M = length(chains);
            
            % 3. rhat below is a P-by-1 vector.
            rHat = stats.mcmc.utils.computeRHat(chains);
            
            % 4. Compute effective sample size for each variable for a
            % given chain and then add up the values across chains.
            ess = zeros(P,M);
            
            for i = 1:M
                ess(:,i) = stats.mcmc.utils.computeESS(chains{i},maxlag,historysize);
            end
            ess = sum(ess,2);
            
            % 5. Get posterior mean and variance using all samples.
            chains     = chains(:)';
            samples    = cell2mat(chains);
            muHat      = mean(samples,2);
            sigma2Hat  = var(samples,0,2);
            quantile5  = quantile(samples,0.05,2);
            quantile95 = quantile(samples,0.95,2);
            
            % 6. Estimate the monte carlo variance - this is the variance
            % of muHat.
            varMuHat = sigma2Hat./ess;
            
            % 7. Put stuff into a table.
            names = this.VariableNames;
            
%             fprintf('%s\n',getString(message('stats:mcmc:impl:MCMCImpl:NumberOfChains',M)));
            
            statstbl      = table();
            statstbl.Name = names;
            statstbl.Mean = muHat;
            statstbl.MCSE = sqrt(varMuHat);
            statstbl.SD   = sqrt(sigma2Hat);
            statstbl.Q5   = quantile5;
            statstbl.Q95  = quantile95;
            statstbl.ESS  = ess;
            statstbl.RHat = rHat;
        end
    end
    
end