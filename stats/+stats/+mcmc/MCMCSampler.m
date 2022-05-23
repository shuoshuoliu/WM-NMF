classdef (Abstract) MCMCSampler < classreg.learning.internal.DisallowVectorOps
%

%   Copyright 2016-2018 The MathWorks, Inc.    
    
    properties(GetAccess=public,SetAccess=protected,Dependent)
        %LogPDF - Function handle to evaluate the log posterior.
        %   The LogPDF property is a function handle that should be
        %   callable like this:
        %
        %       LPDF = LOGPDF(X)         (for non-gradient based samplers)
        %
        %   or like this:
        %
        %       [LPDF,GLPDF] = LOGPDF(X) (for gradient based samplers)
        %
        %   where X is a P-by-1 vector, LPDF is a scalar representing the
        %   log probability density (up to an additive constant) and GLPDF
        %   is a P-by-1 vector containing the gradient of LOGPDF at X. The
        %   length of input X accepted by LOGPDF is the same as the length
        %   of property StartPoint.
        %
        %   See also StartPoint, VariableNames, estimateMAP, hmcSampler.
        LogPDF;
        
        %VariableNames - Variable names.
        %   The VariableNames property is a cell array of length P-by-1
        %   where P is the length of StartPoint property. The name for
        %   element i of StartPoint is taken from element i of
        %   VariableNames.
        %
        %   See also LogPDF, StartPoint, diagnostics, hmcSampler.
        VariableNames;
    end
    
    properties(Dependent)
        %StartPoint - Start point for the sampler.
        %   The StartPoint property is a P-by-1 vector of variables
        %   accepted by property LogPDF. It represents the start point of
        %   the Markov chain created by the sampler.
        %
        %   See also LogPDF, VariableNames, estimateMAP, hmcSampler.
        StartPoint;
    end
    
    properties(GetAccess=public,SetAccess=protected,Dependent,Hidden)
        %HasGradient - True if LogPDF can compute gradient.
        %   The HasGradient property is a logical scalar that is true if
        %   LogPDF can compute gradient of the log posterior and false
        %   otherwise.
        HasGradient;
        
        %Name - Name of the sampler.
        %   The Name property is a character vector representing the name
        %   of the sampler.
        Name;
        
        %SamplerParameters - Parameters for the sampler.
        %   The SamplerParameters property is a structure containing the
        %   tuning parameters for the sampler.
        SamplerParameters;
    end
    
    properties(Abstract,Hidden)
        %Impl - Implementation class for MCMC.
        %   The Impl property is of type MCMCImpl. It provides an
        %   implementation of MCMC. In particular, the Impl class should
        %   provide methods for drawing samples and tuning sampler
        %   parameters.
        Impl;
    end
    
    methods
        function logPDF = get.LogPDF(this)
            logPDF = this.Impl.LogPDF;
        end
        
        function hasGradient = get.HasGradient(this)
            hasGradient = this.Impl.HasGradient;
        end
        
        function start = get.StartPoint(this)
            start = this.Impl.StartPoint;
        end
        
        function this = set.StartPoint(this,start)
            this = setStart(this,start);
        end
        
        function name = get.Name(this)
            name = this.Impl.Name;
        end
        
        function samplerParameters = get.SamplerParameters(this)
            samplerParameters = toStruct(this.Impl.SamplerParameters);
        end
        
        function variableNames = get.VariableNames(this)
            variableNames = this.Impl.VariableNames;
        end
    end
       
    methods(Hidden)
        function this = setStart(this,start)
            p                    = length(this.StartPoint);
            start                = stats.mcmc.utils.InputValidator.validateStart(start,p);
            this.Impl.StartPoint = start;
        end
    end
    
    methods(Hidden,Static)
        function names = getDefaultVariableNames(numVariables)
            names = cell(numVariables,1);
            for i = 1:numVariables
                names{i} = ['x',num2str(i)];
            end
        end
    end
    
    methods(Abstract)
        drawSamples(this,varargin);        
        tuneSampler(this,varargin);                        
    end
        
    methods
        function [xHat,fitInfo] = estimateMAP(this,varargin)
%estimateMAP - Computes the maximum of the posterior log PDF.
%   [XHAT,FITINFO] = estimateMAP(SMP) takes an object SMP of type
%   MCMCSampler, maximizes the posterior log PDF and returns the estimated
%   solution in XHAT and fitting information in FITINFO.
%
%   o Output FITINFO is a structure with these fields:
%       FieldName           Meaning
%       Iteration         - Iteration index.
%       Objective         - Negative log PDF.
%       Gradient          - Final gradient of negative log PDF.
%
%   [...] = estimateMAP(SMP,'Name','Value',...) accepts additional
%   name/value pairs like this:
%
%       Name                Value
%       'StartPoint'     -  A P-by-1 vector specifying the initial point
%                           to start the optimization where P is the length
%                           of SMP.StartPoint. Default is SMP.StartPoint.
%       'IterationLimit' -  A scalar integer specifying the maximum number
%                           of iterations for optimization. Default is
%                           1000.
%       'VerbosityLevel' -  A non-negative integer specifying the verbosity
%                           level as follows:
%                             * 0    - no convergence summary is displayed.
%                             * >=1  - convergence summary is displayed on
%                                      screen.
%                           Default is 0.
%       'GradientTolerance' 
%                         - A positive real scalar specifying the relative
%                           convergence tolerance on the norm of the
%                           objective function gradient. Default is 1e-6.
%       'StepTolerance'   - A positive real scalar specifying the absolute
%                           convergence tolerance on the step size. Default
%                           is 1e-6.
%
%   See also diagnostics, hmcSampler.

%       'HasGradient'    -  A logical scalar indicating whether log
%                           posterior used to create this object can return
%                           gradient information (true) or not (false).
%                           Default is SMP.HasGradient.

            % 1.1 Set parameter defaults.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            dfltStartPoint        = this.StartPoint;
            dfltHasGradient       = this.HasGradient;
            dfltIterationLimit    = 1000;
            dfltVerbosityLevel    = 0;
            dfltGradientTolerance = 1e-6;
            dfltStepTolerance     = 1e-6;
            
            % 1.2 Parse optional name/value pairs.
            paramNames = {  'StartPoint',   'HasGradient',   'IterationLimit',   {'Verbose','VerbosityLevel'},   'GradientTolerance',   'StepTolerance'};
            paramDflts = {dfltStartPoint, dfltHasGradient, dfltIterationLimit,             dfltVerbosityLevel, dfltGradientTolerance, dfltStepTolerance};
            [x0,haveGrad,iterationLimit,verbose,gradTol,stepTol] = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});
            
            % 1.3 Validate x0, haveGrad, iterationLimit, verbose.
            x0             = stats.mcmc.utils.InputValidator.validateStart(x0,length(dfltStartPoint));
            haveGrad       = stats.mcmc.utils.InputValidator.validateHasGradient(haveGrad);
            iterationLimit = stats.mcmc.utils.InputValidator.validateIterationLimit(iterationLimit);
            verbose        = stats.mcmc.utils.InputValidator.validateVerbose(verbose);
            gradTol        = stats.mcmc.utils.InputValidator.validateGradientTolerance(gradTol);
            stepTol        = stats.mcmc.utils.InputValidator.validateStepTolerance(stepTol);
            
            % 2. Call estimateMAP on the Impl class.
            [xHat,fitInfo] = estimateMAP(this.Impl,x0,haveGrad,iterationLimit,verbose,gradTol,stepTol);
        end
        
        function statstbl = diagnostics(this,chains,varargin)
%DIAGNOSTICS - Computes convergence diagnostics for MCMC.
%   TBL = DIAGNOSTICS(SMP,CHAINS) computes MCMC diagnostics using the
%   matrix or cell array CHAINS. If CHAINS is a matrix, it is of size
%   N-by-P where N is the number of samples in P dimensions where P is
%   equal to the length of SMP.StartPoint. If CHAINS is a cell array,
%   CHAINS{i} is a Ni-by-P matrix of Ni samples in P dimensions. TBL is a
%   table containing MCMC diagnostics.
%
%   [...] = DIAGNOSTICS(SMP,CHAINS,'Name','Value',...) accepts additional
%   name/value pairs like this:
%
%       Name                Value
%       'MaxLag'         -  A positive integer specifying the maximum
%                           number of lags to use when computing effective
%                           sample sizes. When CHAINS is a matrix, the
%                           maximum allowable lag for N samples is (N-1).
%                           When CHAINS is a cell array, the maximum
%                           allowable lag for CHAINS{i} with Ni samples is
%                           (Ni-1). Default is 100.
%
%   See also estimateMAP, hmcSampler.

            % 1. Validate chains. If chains is a matrix on input then
            % chains below will be a length 1 cell array. Input chains has
            % elements of size Ni-by-P but chains below has elements of
            % size P-by-Ni.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            p      = length(this.StartPoint);
            chains = stats.mcmc.utils.InputValidator.validateChains(chains,p);
            
            % 2.1 Set parameter defaults.
            dfltMaxLag = 100;                        
            % 2.2 Parse optional name/value pairs.
            paramNames = {  'MaxLag'};
            paramDflts = {dfltMaxLag};
            [maxlag] = internal.stats.parseArgs(paramNames,paramDflts,varargin{:});            
            % 2.3 Validate maxlag.
            maxlag = stats.mcmc.utils.InputValidator.validateMaxLag(maxlag);
            
            % 3. Call diagnostics on the Impl class.
            statstbl = diagnostics(this.Impl,chains,maxlag);
        end
    end
    
end