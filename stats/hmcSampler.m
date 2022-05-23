function hmc = hmcSampler(logpdf,start,varargin)
%hmcSampler - Create a Hamiltonian Monte Carlo (HMC) sampler.
%   hmc = hmcSampler(LOGPDF,STARTPOINT) creates a HamiltonianSampler object
%   that can be used for Hamiltonian Monte Carlo (HMC) sampling. LOGPDF is
%   a function handle that evaluates the log probability density (up to an
%   additive constant) of the desired equilibrium distribution. STARTPOINT
%   is a P-by-1 vector representing the initial point from which to start
%   HMC sampling.
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
%   o Tune the HMC sampler using the tuneSampler method.
%
%   o Draw samples using the drawSamples method.
%
%   o Estimate the maximum of LOGPDF using the estimateMAP method.
%
%   o Assess convergence using the diagnostics method.
%
%   hmc = hmcSampler(LOGPDF,STARTPOINT,'Name','Value',...) specifies
%   additional name/value pairs like this:
%
%       Name            Value 'StepSize'   -  A scalar specifying the step
%       size to use for
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
%                       analytical gradient of LOGPDF at STARTPOINT vs. the
%                       numerical gradient or not. Default is true.
%       'VariableNames'
%                    -  A string array or cell array of length P specifying
%                    variable names
%                       where P is the length of STARTPOINT. The name for
%                       element i of STARTPOINT is taken from element i of
%                       VariableNames. Default is to use variable names x1,
%                       x2, etc.
%       'UseNumericalGradient'
%                    -  A logical scalar that is either true (or 1) or
%                       false (or 0) that specifies whether to use the
%                       numerical gradient of LOGPDF or not. If true,
%                       LOGPDF is called like LPDF = LOGPDF(X) so that
%                       LOGPDF need not return gradient as the second
%                       output. Default is false.
%
%   See also stats.mcmc.HamiltonianSampler, slicesample, mhsample.

%   Copyright 2016 The MathWorks, Inc.

if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

hmc = stats.mcmc.HamiltonianSampler(logpdf,start,varargin{:});
end