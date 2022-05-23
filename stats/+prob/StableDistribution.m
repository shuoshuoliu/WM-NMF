classdef StableDistribution < prob.ToolboxFittableParametricDistribution
    %StableDistribution Stable probability distribution.
    %    An object of the StableDistribution class represents a stable
    %    probability distribution with specific values for the ALPHA, BETA,
    %    GAM and DELTA parameters. This distribution object can be created
    %    directly using the MAKEDIST function or fit to data using the
    %    FITDIST function.
    %
    %    StableDistribution methods:
    %       cdf                   - Cumulative distribution function
    %       icdf                  - Inverse cumulative distribution function
    %       iqr                   - Interquartile range
    %       mean                  - Mean
    %       median                - Median
    %       negloglik             - Negative log likelihood function
    %       paramci               - Confidence intervals for parameters
    %       pdf                   - Probability density function
    %       proflik               - Profile likelihood function
    %       random                - Random number generation
    %       std                   - Standard deviation
    %       truncate              - Truncation distribution to an interval
    %       var                   - Variance
    %
    %    StableDistribution properties:
    %       DistributionName      - Name of the distribution
    %       alpha                 - Value of the alpha parameter (first shape parameter)
    %       beta                  - Value of the beta parameter (second shape parameter)
    %       gam                   - Value of the gam parameter (scale)
    %       delta                 - Value of the delta parameter (location)
    %       NumParameters         - Number of parameters
    %       ParameterNames        - Names of parameters
    %       ParameterDescription  - Descriptions of parameters
    %       ParameterValues       - Vector of values of parameters
    %       Truncation            - Two-element vector indicating truncation limits
    %       IsTruncated           - Boolean flag indicating if distribution is truncated
    %       ParameterCovariance   - Covariance matrix of estimated parameters
    %       ParameterIsFixed      - Four-element boolean vector indicating fixed parameters
    %       InputData             - Structure containing data used to fit the distribution
    %
    %    See also fitdist, makedist.
    
    %    Copyright 2015-2019 The MathWorks, Inc.
    
    properties(Dependent=true)
        %ALPHA Value of ALPHA (shape) parameter
        %      The ALPHA property represents the first of two shape parameters of the
        %      stable distribution. It represents the characteristic exponent of the stable
        %      distribution.
        %
        %      See also BETA, GAM, DELTA.
        alpha
        
        %BETA Value of BETA (shape) parameter
        %     The BETA property represents the second of two shape parameters of the
        %     stable distribution. It represents the skewness of the stable
        %     distribution.
        %
        %     See also ALPHA, GAM, DELTA.
        beta
        
        %GAM Value of GAM (scale) parameter
        %    The GAM property represents the scale parameter of the stable
        %    distribution.
        %
        %    See also ALPHA, BETA, DELTA.
        gam
        
        %DELTA Value of DELTA (location) parameter
        %      The DELTA property represents the location parameter of the stable
        %      distribution.
        %
        %      See also ALPHA, BETA, GAM.
        delta
    end
    properties(GetAccess='public',Constant=true)
        %DistributionName Distribution name.
        %    The DistributionName property indicates the name of the probability
        %    distribution.
        %
        %    See also ParameterNames, ParameterValues.
        DistributionName = getString(message('stats:dfittool:NameStable'));
        
        %NumParameter Number of parameters.
        %    NumParameters is the number of parameters in the distribution.
        %
        %    See also ParameterValues.
        NumParameters = 4;
        
        %ParameterNames Parameter names.
        %    ParameterNames is a cell array of strings containing the names of the
        %    parameters of the probability distribution.
        %
        %    See also ParameterValues, ParameterDescription.
        ParameterNames = {'alpha' 'beta' 'gam' 'delta'};
        
        %ParameterDescription Parameter description.
        %    ParameterNames is a cell array of strings containing short
        %    descriptions of the parameters of the probability distribution.
        %
        %    See also ParameterNames, ParameterValues.
        ParameterDescription = {getString(message('stats:probdists:ParameterDescriptionShape1')) ...
            getString(message('stats:probdists:ParameterDescriptionShape2')) ...
            getString(message('stats:probdists:ParameterDescriptionScale')) ...
            getString(message('stats:probdists:ParameterDescriptionLocation'))};
    end
    properties(GetAccess='public',SetAccess='protected')
        %ParameterValues Parameter values.
        %    ParameterVales is a vector containing the values of the parameters of
        %    the probability distribution.
        %
        %    See also ALPHA, BETA, GAM, DELTA.
        ParameterValues
    end
    methods(Hidden)
        function pd = StableDistribution(alpha,beta,gam,delta)
            if nargin==0 % default parameter values give the Normal distribution N(0,2)
                alpha = 2;
                beta = 0;
                gam = 1;
                delta = 0;
            end
            checkargs(alpha,beta,gam,delta)
            
            pd.ParameterValues = [alpha beta gam delta];
            pd.ParameterIsFixed = [true true true true];
            pd.ParameterCovariance = zeros(pd.NumParameters);
        end
    end
    methods
        function m = mean(this)
            requireScalar(this)
            if this.IsTruncated
                m = truncatedMoment(this,1);
                return
            end
            m = NaN;
            % mean exists when 1< alpha <=2, proposition 1.13 in J. P. Nolan (2015)
            if (this.alpha>1) && (this.alpha<=2)
                m = this.delta-this.beta.*this.gam.*tan(pi*this.alpha/2);
            end
        end
        function v = var(this)
            requireScalar(this)
            if this.IsTruncated
                v = truncatedMoment(this,2);
                return
            end
            v = NaN;
            % variance exists only when alpha == 2, which is the Normal
            % distribution.
            if this.alpha == 2
                v = 2.*this.gam.^2;
            end
        end
        function this = set.alpha(this,alpha)
            checkargs(alpha,this.beta,this.gam,this.delta)
            this.ParameterValues(1) = alpha;
            this = invalidateFit(this);
        end
        function this = set.beta(this,beta)
            checkargs(this.alpha,beta,this.gam,this.delta)
            this.ParameterValues(2) = beta;
            this = invalidateFit(this);
        end
        function this = set.gam(this,gam)
            checkargs(this.alpha,this.beta,gam,this.delta)
            this.ParameterValues(3) = gam;
            this = invalidateFit(this);
        end
        function this = set.delta(this,delta)
            checkargs(this.alpha,this.beta,this.gam,delta)
            this.ParameterValues(4) = delta;
            this = invalidateFit(this);
        end
        function alpha = get.alpha(this)
            alpha = this.ParameterValues(1);
        end
        function beta = get.beta(this)
            beta = this.ParameterValues(2);
        end
        function gam = get.gam(this)
            gam = this.ParameterValues(3);
        end
        function delta = get.delta(this)
            delta = this.ParameterValues(4);
        end
    end
    methods(Static,Hidden)
        function pd = fit(varargin)
            %FIT Fit distribution to data.
            %    FIT is a static method that fits the stable distribution to data.
            %    Fitting requires the Statistics and Machine Learning Toolbox. You
            %    should call the FITDIST function instead of calling this method directly.
            %
            % 	 The FIT function of stable distribution interpolates a precomputed
            %    density table for ALPHA>=0.4 instead of the direct integration method to
            %	 find the maximum likelihood estimation. Estimates from McCulloch's quantile
            %	 method is used as an initial values of the optimization procedure.
            %
            % Reference:
            %   [1] J. P. Nolan (2001) "Maximum likelihood estimation and diagnostics for
            %       stable distributions." Lévy processes. Birkhäuser Boston, p379-400.
            %
            %    See also FITDIST.
            
            %    The density table used in the stable fitting function is calculated
            %    by internal.stats.StablePdfLookupTable.
            
            [xOriginal,cens,freq,opt] = prob.ToolboxFittableParametricDistribution.processFitArgs(varargin{:});
            x = prob.ToolboxFittableParametricDistribution.removeCensoring(xOriginal,cens,freq,'stable');
            params = stablefit(x,0.05,opt);
            [nll,cov] = stablelike(params,x);
            pd = prob.StableDistribution.makeFitted(params,nll,cov,xOriginal,cens,freq);
        end
        function varargout = likefunc(varargin)
            [varargout{1:nargout}] = stablelike(varargin{:});
        end
        function varargout = cdffunc(varargin)
            [varargout{1:nargout}] = stablecdf(varargin{:});
        end
        function varargout = pdffunc(varargin)
            [varargout{1:nargout}] = stablepdf(varargin{:});
        end
        function varargout = invfunc(varargin)
            [varargout{1:nargout}] = stableinv(varargin{:});
        end
        function varargout = randfunc(varargin)
            [varargout{1:nargout}] = stablernd(varargin{:});
        end
        function pd = makeFitted(params,nll,cov,x,cens,freq)
            c = num2cell(params);
            pd = prob.StableDistribution(c{:});
            pd.NegativeLogLikelihood = nll;
            pd.ParameterCovariance = cov;
            pd.ParameterIsFixed = [false false false false];
            pd.InputData = struct('data',x,'cens',cens,'freq',freq);
        end
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.StableDistribution');
            info.name = prob.StableDistribution.DistributionName;
            info.code = 'stable';
            info.optimopts = true;
            info.islocscale = false;
            info.logci = [false false false false];
            info.plim = [0 -1 0 -Inf;2 1 Inf Inf];
        end
    end
end % classdef

function checkargs(alpha,beta,gam,delta)
if ~(isscalar(alpha) && isnumeric(alpha) && isreal(alpha) && isfinite(alpha) && alpha>0 && alpha<=2)
    error(message('stats:probdists:PositiveParameterLE2','ALPHA'))
end
if ~(isscalar(beta) && isnumeric(beta) && isreal(beta) && isfinite(beta) && beta>=-1 && beta<=1)
    error(message('stats:probdists:ParameterN1ToP1','BETA'))
end
if ~(isscalar(gam) && isnumeric(gam) && isreal(gam) && isfinite(gam) && gam>0)
    error(message('stats:probdists:PositiveParameter','GAM'))
end
if ~(isscalar(delta) && isnumeric(delta) && isreal(delta) && isfinite(delta))
    error(message('stats:probdists:ScalarParameter','DELTA'))
end
end

%===========Stable Distribution Functions=============

function p = stablecdf(x,alpha,beta,gam,delta,uflag)
%STABLECDF Stable cumulative distribution function (cdf).
%
% Reference:
%       [1] J. P. Nolan (1997) "Numerical Calculation of Stable Densities
%           and Distribution Functions", Communications in statistics,
%           Stochastic models, 13(4), 759-774.

[err, x, alpha, beta, gam, delta] = distchck(5, x, alpha, beta, gam, delta);

if err > 0
    error(message('stats:addstable:InputSizeMismatch'));
end

% Return NaN for illegal parameter values.
alpha(alpha<=0 | alpha>2) = NaN;
beta(beta<-1 | beta>1) = NaN;
gam(gam <= 0) = NaN;

% Initialize P to zero.
p = zeros(size(x));

% Standardize data
z = (x - delta)./gam;

% Simple cases
Gau = alpha == 2;
Cau = (alpha == 1) & (beta == 0);
Levy = (alpha == 0.5) & (abs(beta) == 1);
if nargin > 5
    % Upper tail calculation
    if ~strcmpi(uflag, 'upper')
        error(message('stats:cdf:UpperTailProblem'));
    else
        if any(Gau(:))  % Gaussian distribution
            p(Gau) = 0.5*erfc(z(Gau)/2);
        end
        if any(Cau(:))  % Cauchy distribution
            p(Cau) = tcdf(z(Cau),1,'upper');
        end
        if any(Levy(:)) % Levy distribution
            LbetaE1 = Levy & (beta==1);
            if any(LbetaE1(:))
                % S(0.5,1,1,0;0) = Levy(-1,1) which is of mean=-1 and var = 1
                % so, subtract mean (-1) from data
                z(LbetaE1) = z(LbetaE1) + 1;
                idx1 = LbetaE1 & (z>0);
                idx2 = LbetaE1 & (z<=0);
                p(idx1) = erf(sqrt(1./(2*z(idx1))));
                p(idx2) = 1;
            end
            LbetaN1 = Levy & (beta==-1);
            if any(LbetaN1(:))
                % p(LbetaN1) = 1 - stablecdf(-z(LbetaN1),0.5,1,1,0); % Reflection property
                % To avoid truncation error, use exact formula instead.
                z(LbetaN1) = -z(LbetaN1) + 1;
                idx1 = LbetaN1 & (z>0);
                idx2 = LbetaN1 & (z<=0);
                p(idx1) = erfc(sqrt(1./(2*z(idx1))));
                p(idx2) = 0;
            end
        end
        others = (~Gau) & (~Cau) & (~Levy);
        if any(others(:))
            p(others) = stablecdf(-z(others),alpha(others),-beta(others),1,0); % Reflection property
        end
    end
else
    if any(Gau(:))  % Gaussian distribution
        p(Gau) = 0.5*(1 + erf(z(Gau)/2));
    end
    
    if any(Cau(:))  % Cauchy distribution
        p(Cau) = tcdf(z(Cau),1); % Cauchy(0,1) is t-distribution with d.f. = 1
    end
    
    if any(Levy(:)) % Levy distribution
        LbetaE1 = Levy & (beta==1);
        if any(LbetaE1(:))
            % S(0.5,1,1,0;0) = Levy(-1,1) which is of mean=-1 and var = 1
            % so, subtract mean (-1) from data
            z(LbetaE1) = z(LbetaE1) + 1;
            idx1 = LbetaE1 & (z>0);
            idx2 = LbetaE1 & (z<=0);
            p(idx1) = erfc(sqrt(1./(2*z(idx1))));
            p(idx2) = 0;
        end
        LbetaN1 = Levy & (beta==-1);
        if any(LbetaN1(:))
            % p(LbetaN1) = 1 - stablecdf(-z(LbetaN1),0.5,1,1,0); % Reflection property
            % To avoid truncation error, use exact formula instead.
            z(LbetaN1) = -z(LbetaN1) + 1;
            idx1 = LbetaN1 & (z>0);
            idx2 = LbetaN1 & (z<=0);
            p(idx1) = erf(sqrt(1./(2*z(idx1))));
            p(idx2) = 1;
        end
    end
    
    % When alpha is very close to 1 but not equal to 1, integral in the
    % distribution function change very rapidly and it is difficult to
    % be approximated accurately. Thus, Alpha with values |1-alpha|<0.02+eps
    % is rounded to 1.
    alphaN1 = (abs(alpha - 1) >= 0.02+eps) & (~Gau) & (~Cau) & (~Levy);
    alphaR1 = (abs(alpha - 1) < 0.02+eps) & alpha~=1;
    if any(alphaR1(:))
        alpha(alphaR1) = 1;
    end
    
    % General Case: alpha ~= 1, i.e. abs(alpha-1)>=0.02+eps
    if any(alphaN1(:))
        zeta = -beta .* tan(pi*alpha/2);
        theta0 = (1./alpha) .* atan(-zeta);
        % Make sure theta0 is between -pi/2 and pi/2 caused by roundoff
        theta0(theta0>pi/2) = pi/2;
        theta0(theta0<-pi/2) = -pi/2;
        
        % There are three cases z>zeta, z=zeta, and z<zeta. In practice, z
        % can be very close to zeta but not exactly equal to zeta. Thus,
        % we consider z with values abs(z-zeta)<=4*eps (an empirical value)
        % as z=zeta.
        zEzeta = (abs(z-zeta) <= 4*eps) & alphaN1;
        zGzeta = (z - zeta > 4*eps) & alphaN1;
        zLzeta = (zeta - z > 4*eps) & alphaN1;
        
        if any(zEzeta(:))
            p(zEzeta) = 0.5 - theta0(zEzeta)/pi;
        end
        
        if any(zGzeta(:))
            % Define V function. Use change of variables, theta = (pi/2 + theta0)*theta - theta0,
            % to make the integral goes from 0 to 1:
            A1 = alpha.*theta0;
            A2 = cos(A1).^(1./(alpha-1));
            c1 = (alpha > 1 & alphaN1) + (alpha < 1 & alphaN1).*(0.5 - theta0/pi);
            exp1 = alpha./(alpha-1);
            zshift = (z- zeta).^(alpha./(alpha - 1));
            thetashift1 = pi/2 + theta0(zGzeta);
            V = @(theta) real(A2(zGzeta) .* ( cos(theta) ./ sin( alpha(zGzeta).*(theta + theta0(zGzeta)) ) ).^exp1(zGzeta)).*...% in case of complex number
                cos( A1(zGzeta) + (alpha(zGzeta)-1).*theta ) ./ cos(theta);
            % Calculate cumulative densities by integration
            p( zGzeta ) = c1(zGzeta) + sign(1-alpha(zGzeta))/pi.* thetashift1.* ...
                integral(@(theta) exp(-zshift(zGzeta).*V(thetashift1 .* theta -  theta0(zGzeta))),...
                0,1,'ArrayValued',1,'AbsTol',1e-6,'RelTol',1e-4); % Default accuracy of INTEGRAL is sometimes hard to reach
        end
        
        if any(zLzeta(:))
            p(zLzeta) = 1 - stablecdf(-z(zLzeta),alpha(zLzeta),-beta(zLzeta),1,0); % Reflection property
        end
    end
    
    % General Case: alpha = 1 and beta > 0
    betaG0 = (alpha == 1) & (beta > 0);
    if any(betaG0(:))
        % Define logV instead of V function to avoid overflow/underflow.
        % Use change of variables, theta = pi/2*theta, to make the integral goes from -1 to 1
        logV = @(theta) log(2/pi*((pi/2 + beta(betaG0).*theta)./cos(theta))) + ...
            ( 1./beta(betaG0).*(pi/2 + beta(betaG0).*theta) .* tan(theta) );
        % Calculate cumulative densities by integration
        p(betaG0) = 0.5*integral(@(theta) exp(-exp(-pi*z(betaG0)./(2*beta(betaG0))+logV(pi/2*theta))),...
            -1,1,'ArrayValued',1,'AbsTol',1e-6,'RelTol',1e-4); % Default accuracy of INTEGRAL is sometimes hard to reach
    end
    
    % General Case: alpha = 1 and beta < 0. Note: alpha=1 and beta = 0 is the special case: Cauchy
    betaL0 = (alpha == 1) & (beta < 0);
    if any(betaL0(:))
        p(betaL0) = 1 - stablecdf(-z(betaL0),alpha(betaL0),-beta(betaL0),1,0); % Reflection property
    end
    
end

% In case of small imaginary or negative results from INTEGRAL
p = max(real(p),0);
end

function y = stablepdf(x,alpha,beta,gam,delta,varargin)
%STABLEPDF Stable probability density function (pdf).
%
% Reference:
%       [1] J. P. Nolan (1997) "Numerical Calculation of Stable Densities
%           and Distribution Functions", Communications in statistics,
%           Stochastic models, 13(4), 759-774.

[err, x, alpha, beta, gam, delta] = distchck(5, x, alpha, beta, gam, delta);

if err > 0
    error(message('stats:addstable:InputSizeMismatch'));
end

% Return NaN for illegal parameter values.
alpha(alpha<=0 | alpha>2) = NaN;
beta(beta<-1 | beta>1) = NaN;
gam(gam <= 0) = NaN;

% Initialize Y to zero.
y = zeros(size(x));

% Standardize data
z = (x - delta);

% Round small values to 0 for stability
z(abs(z) <= sqrt(eps)) = 0;
z = z./gam;

% Check for simple cases
Gau = alpha == 2;  % Gaussian distribution
if any(Gau(:))
    y(Gau) = 1/sqrt(4*pi).*exp( -.25 * z(Gau).^2 );
end

Cau = (alpha==1) & (beta == 0);  % Cauchy distribution
if any(Cau(:))
    y(Cau) = tpdf(z(Cau),1); % Cauchy(0,1) is t-distribution with d.f. = 1
end

Levy = (alpha ==0.5) & (abs(beta) == 1); % Levy distribution
if any(Levy(:))
    LbetaP1 = Levy & (beta == 1);
    LbetaN1 = Levy & (beta ==-1);
    
    if any(LbetaP1(:))
        % S(0.5,1,1,0;0) = Levy(-1,1) which is of mean=-1 and var = 1
        % so, subtract mean (-1) from data
        z(LbetaP1) = z(LbetaP1) + 1;
        idx1 = LbetaP1 & (z>0);
        idx2 = LbetaP1 & (z<=0);
        y(idx1) = sqrt(1/(2*pi)) .* exp(-0.5./z(idx1)) ./ (z(idx1).^1.5);
        y(idx2) = 0;
    end
    
    if any(LbetaN1(:))
        y(LbetaN1) = stablepdf(-z(LbetaN1),alpha(LbetaN1),1,1,0); % reflection property
    end
end

% When alpha is very close to 1 or 0, the approximation to integral may not
% be accurate.  When alpha is very close to 1 but not equal to 1, integral in
% the distribution function change very rapidly and it is difficult to
% be approximated accurately. Thus, Alpha with values |1-alpha|<0.02+eps
% is rounded to 1.
alphaN1 = (abs(alpha - 1) >= 0.02+eps) & (~Gau) & (~Cau) & (~Levy);
alphaR1 = (abs(alpha - 1) < 0.02+eps) & alpha~=1;
if any(alphaR1(:))
    alpha(alphaR1) = 1;
end

% General Case: alpha ~= 1
if any(alphaN1(:))
    zeta = -beta .* tan(pi*alpha/2);
    theta0 = (1./alpha) .* atan(-zeta);
    % Make sure theta0 is between -pi/2 and pi/2 caused by roundoff
    theta0(theta0>pi/2) = pi/2;
    theta0(theta0<-pi/2) = -pi/2;
    % There are three cases z>zeta, z=zeta, and z<zeta. In practice, z
    % can be very close to zeta but not exactly equal to zeta. Thus,
    % we consider z with values abs(z-zeta)<=sqrt(eps) as z=zeta.
    zEzeta = (abs(z-zeta)<= 4*eps) & alphaN1;
    zGzeta = (z - zeta > 4*eps) & alphaN1;
    zLzeta = (zeta - z > 4*eps) & alphaN1;
    
    if any(zGzeta(:))
        A1 = alpha.*theta0;
        A2 = cos(A1).^(1./(alpha-1));
        exp1 = alpha./(alpha-1);
        zshift = (z - zeta).^exp1;
        
        % In some cases, the integrand can be very peaked, so use a
        % bisectionSlover to locate the peak and then evaluate the integral
        % in two pieces increase the accuracy of the integral.
        
        % First, define Vcol function with columned variables to be used in bisectionSolver
        Vcol = @(theta) A2(:) .* ( cos(theta) ./ sin( alpha(:).*(theta + theta0(:)) ) ).^exp1(:).*...
            cos( A1(:) + (alpha(:)-1).*theta ) ./ cos(theta);
        
        % Second, find the location of peak, theta2, of the integrand by solving g(theta)-1 = 0.
        g = @(theta) zshift(:).*Vcol(theta) - 1;
        R = [-theta0(:), pi/2*ones(size(theta0(:)))];
        theta2 = bisectionSolver(g,R,alpha(:));
        theta2 = reshape(theta2,size(theta0));
        
        % Third, define V function used in the integration. Change variables
        % so the two integrals go from 0 to 1/2 and 1/2 to 1.
        theta2shift1 = 2*(theta2(zGzeta) + theta0(zGzeta));
        theta2shift2 = 2*(pi/2 - theta2(zGzeta));
        V = @(theta) real(A2(zGzeta) .* ( cos(theta) ./ sin( alpha(zGzeta).*(theta + theta0(zGzeta)) ) ).^exp1(zGzeta))...% in case of complex number
            .*cos( A1(zGzeta) + (alpha(zGzeta)-1).*theta ) ./ cos(theta);
        g1 = @(theta)  zshift(zGzeta) .* V(theta2shift1 .* theta - theta0(zGzeta));
        g2 = @(theta)  zshift(zGzeta) .* V(theta2shift2 .* (theta - .5) + theta2(zGzeta));
        zexpz = @(z) max(0,z .* exp(-z)); % use max in case of NaN
        
        % Fourth, calculate the density by integration
        c2 = alpha(zGzeta) ./ (pi * abs(alpha(zGzeta) - 1) .* ( z(zGzeta) - zeta(zGzeta)) );
        y(zGzeta) = c2 .* (theta2shift1 .* integral(@(theta) zexpz( g1(theta) ),0,0.5, 'ArrayValued',1,'AbsTol',1e-10,'RelTol',1e-4) ...
            + theta2shift2 .* integral(@(theta) zexpz( g2(theta) ),0.5, 1,'ArrayValued',1,'AbsTol',1e-10,'RelTol',1e-4) ); % Default accuracy of INTEGRAL is sometimes hard to reach
    end
    
    if any(zEzeta(:))
        y(zEzeta) = gamma(1+1./alpha(zEzeta)).*cos(theta0(zEzeta))./ (pi.*(1+zeta(zEzeta).^2).^(1./(2*alpha(zEzeta))));
    end
    
    if any(zLzeta(:))
        y(zLzeta) = stablepdf(-z(zLzeta),alpha(zLzeta),-beta(zLzeta),1,0);
    end
end

% General Case: alpha = 1 and beta ~= 0
aE1bN0 = (alpha==1) & (beta~=0);
if any(aE1bN0(:))
    zexp = pi*z./(2*beta);
    
    % Find the peak of the integrand, use logs to avoid overflow/underflow
    logVcol = @(theta) log(2/pi*((pi/2 + beta(:).*theta)./cos(theta))) + ...
        ( 1./beta(:).*(pi/2 + beta(:).*theta) .* tan(theta) );
    logg = @(theta) -zexp(:) + logVcol(theta);
    R = repmat([-pi/2, pi/2],size(zexp(:),1),1);
    theta2 = bisectionSolver(logg,R,1-beta(:));
    theta2 = reshape(theta2,size(zexp));
    
    % Change variables so the two integrals go from 0 to 1/2 and 1/2 to 1.
    theta2shift1 = 2*(theta2(aE1bN0) + pi/2);
    theta2shift2 = 2*(pi/2 - theta2(aE1bN0));
    logV = @(theta) log(2/pi*((pi/2 + beta(aE1bN0).*theta)./cos(theta))) + ...
        ( 1./beta(aE1bN0).*(pi/2 + beta(aE1bN0).*theta) .* tan(theta) );
    logg1 = @(theta)  -zexp(aE1bN0) + logV(theta2shift1 .* theta - pi/2);
    logg2 = @(theta)  -zexp(aE1bN0) + logV(theta2shift2 .* (theta - 0.5) + theta2(aE1bN0));
    zexpz = @(z) max(0,exp(z) .* exp(-exp(z))); % use max in case of NaN
    y(aE1bN0) = 1./(2*abs(beta(aE1bN0))) .* ...
        (theta2shift1 .* integral(@(theta) zexpz( logg1(theta) ),0 , 0.5, 'ArrayValued',1,'AbsTol',1e-10,'RelTol',1e-4) ...
        + theta2shift2 .* integral(@(theta) zexpz( logg2(theta) ),0.5 , 1, 'ArrayValued',1,'AbsTol',1e-10,'RelTol',1e-4) );  % Default accuracy of INTEGRAL is sometimes hard to reach
    
end
% In case of small imaginary or negative results from INTEGRAL
y = max(real(y),0);
% Rescale from a standard distribution
y = y./gam;
end


function X = bisectionSolver(f,R,alpha,varargin)
%BISECTIONSOLVER A vectorized bisection root solver.
% f is a function handle, R is a N by 2 matrix such that f takes opposite
% sign at each column.

tol = 1e-6;
maxiter = 30;

[N, M] = size(R);
if M ~= 2
    error(message('stats:addstable:BadInputBisec'));
end

a = R(:,1);
b = R(:,2);
X = (a+b)/2;

try
    val = f(X);
catch ME
    error(message('stats:addstable:BadFunBisec'));
end

if size(val,1) ~= N
    error(message('stats:addstable:BadOutBisec'));
end

% Main loop
val = inf;
iter = 0;

while( max(abs(val)) > tol && iter < maxiter )
    X = (a + b)/2;
    val = f(X);
    l = (val > 0);
    alphaG1 = alpha>1;
    if any(alphaG1(:))
        l(alphaG1) = 1-l(alphaG1);
    end
    a = a.*l + X.*(1-l);
    b = X.*l + b.*(1-l);
    iter = iter + 1;
end
end

function x = stableinv(p,alpha,beta,gam,delta)
%STABLEINV Inverse of the stable cumulative distribution function (cdf).
% When alpha is very close to 1 or 0, the approximation to integral in
% density functions may not be accurate. Alpha with values 0<|1-alpha|<0.02+eps
% is rounded to 1 in the calculation of pdf and cdf.

[err, p, alpha, beta, gam, delta] = distchck(5, p, alpha, beta, gam, delta);

if err > 0
    error(message('stats:addstable:InputSizeMismatch'));
end

% Return NaN for out of range parameters or probabilities.
alpha(alpha<=0 | alpha>2) = NaN;
beta(beta<-1 | beta>1) = NaN;
gam(gam <= 0) = NaN;
p(p < 0 | 1 < p) = NaN;

% For Newton's Method.  In general, the Newton method will converge within
% 10 iterations. Non-convergence may happen for very small/large p's.
itermax = 100;
tol = 1e-10;

% Pre-allocation
x = zeros(size(p));

% Gaussian distribution
Gau = alpha == 2;
if any(Gau(:))
    x(Gau) = - 2*erfcinv(2*p(Gau)).*gam(Gau) + delta(Gau);
end

% Cauchy distribution
Cau = (alpha==1) & (beta == 0);
if any(Cau(:))
    % Cauchy(0,1) is t-distribution with d.f. = 1
    x(Cau) =  tinv(p(Cau),1).*gam(Cau) + delta(Cau);
end

% Levy distribution
Levy = (alpha == 0.5) & (abs(beta) == 1);
if any(Levy(:))
    LbetaE1 = Levy & (beta == 1);
    LbetaN1 = Levy & (beta == -1);
    
    if any(LbetaE1(:))
        % S(0.5,1,1,0;0) = Levy(-1,1) which is of mean=-1 and var = 1
        x(LbetaE1) = (0.5*beta(LbetaE1)./(erfcinv(p(LbetaE1))).^2 - 1).*gam(LbetaE1) + delta(LbetaE1);
    end
    
    if any(LbetaN1(:))
        x(LbetaN1) = -stableinv(1-p(LbetaN1),alpha(LbetaN1),-beta(LbetaN1),...
            gam(LbetaN1),-delta(LbetaN1)); % reflection property
    end
end

others = (~Gau) & (~Cau) & (~Levy);
if any(others(:))
    F = zeros(size(p));
    x(others) = intStableinv(alpha(others),beta(others),p(others));
    % Transform from standard distribution
    x(others) = x(others).*gam(others) + delta(others);
    
    % Newton's Method
    F(others) = stablecdf(x(others),alpha(others),beta(others),gam(others),...
        delta(others)) - p(others);
    diff = max(abs(F),0); % max in case of NaNs
    badOthers = (diff > tol) & others;
    iter = 1;
    Fold = F;
    xold = x;
    
    while any(badOthers(:)) && iter < itermax
        % Perform Newton step. If Fprime = 0, step closer to origin instead
        Fprime = stablepdf(x(badOthers),alpha(badOthers),beta(badOthers),...
            gam(badOthers),delta(badOthers));
        x(badOthers) = x(badOthers) - F(badOthers) ./ Fprime;
        blowup = isinf(x(badOthers)) | isnan(x(badOthers));
        if any(blowup(:)==1)
            idx = (isinf(x) | isnan(x)) & badOthers;
            x(idx) = xold(idx) / 2;
        end
        
        % When ALPHA <1, BETA = 1 and -1, the support of the distribution
        % has lower and upper bounds respectively. If Newton method uses
        % values beyond the bound, update the new x by the bound to avoid
        % divergence.
        bdID = (alpha<1) & (abs(beta)==1) & badOthers;
        if any(bdID(:))
            s = sign(beta);
            bd = delta - s.*gam.*tan(pi*alpha/2);
            beyond = (s.*(x-bd) < 0) & badOthers;
            x(beyond) = bd(beyond)+ 1e-10*s(beyond);
        end
        
        F(badOthers) = stablecdf(x(badOthers),alpha(badOthers),beta(badOthers),...
            gam(badOthers),delta(badOthers)) - p(badOthers);
        
        % Make sure we are getting closer; if not, do bisections until we do.
        notCvg = (abs(F) > 1.1*abs(Fold)) & badOthers;
        maxbisecs =(log(abs(x(notCvg)-xold(notCvg)))-log(tol))/log(2);
        maxbisecs = min(max(maxbisecs(:)),400); % use 400 iterations if maxbiscs is too large
        bisecs = 0;
        while any(notCvg(:)) && (bisecs < maxbisecs)
            x(notCvg) = .5*(x(notCvg) + xold(notCvg));
            F(notCvg) = stablecdf(x(notCvg),alpha(notCvg),beta(notCvg),...
                gam(notCvg),delta(notCvg)) - p(notCvg);
            notCvg = (abs(F) > 1.1*abs(Fold)) & badOthers;
            bisecs = bisecs + 1;
        end
        
        % Test for convergence
        diff = max(abs(F),0); % max in case of NaNs
        badOthers = (diff > tol) & others;
        
        % Save for next iteration
        xold = x;
        Fold = F;
        iter = iter + 1;
    end
    if (iter == itermax) && (any(badOthers(:)))
        warning(message('stats:addstable:NotConvStableinv'));
    end
end
end


function X0 = intStableinv(alpha,beta,u)
%INTSTABLEINV Compute initial values used in inverse of the stable cdf.
% Find quantiles of standard stable distributions by interpolation. ALPHA
% and BETA are arrays of same size of the two shape parameters of standard
% stable distribution. U is an array of the same size as ALPHA and BETA
% containing probabilities.
%
% If 0.1 < u < 0.9, interpolates tabulated values to obtain initial value
% If u < 0.1 or u > 0.9, uses asymptotic formulas to make a starting value

% Columnwise to apply interplation
utemp = u(:);
alphaCol = alpha(:);
betaCol = beta(:);
% Pre-allocation
X0 = zeros(size(utemp));

% Asymptotic formula doesn't apply to upper tail if beta=-1 and also
% doesn't apply to lower tail if beta=1.
utemp((betaCol==-1) & (utemp>0.9)) = 0.9;
utemp((betaCol==1) & (utemp<0.1)) = 0.1;

high = (utemp > 0.9);
low = (utemp < 0.1);
middle = (~high) & (~low);

% Use asymptotic formulas to guess upper and lower tails, Thm 1.12 in J. P. Nolan (2015)
if any(high | low)
    C = sin(pi*alphaCol/2).*gamma(alphaCol)/pi;
    X0(high) = ( (1-utemp(high))./(C(high) .* (1 + betaCol(high))) ).^(-1./alphaCol(high));
    X0(low)  = -(utemp(low)./(C(low) .* (1 - betaCol(low)))).^(-1./alphaCol(low));
end

% Use pre-calculated lookup table
if any(middle)
    % Bring ALPHA to table range
    alphaCol(middle) = max(alphaCol(middle),0.1); % Lower limit of the searching table is 0.1
    [Alp,Bet,P] = meshgrid(0.1:0.1:2 , 0:0.2:1 , 0.1:0.1:0.9 );
    stblfrac = zeros(6,20,9);
    stblfrac(:,1:5,1) = ...
        [-1.890857122067030e+006   -1.074884919696010e+003   -9.039223076384694e+001   -2.645987890965098e+001   -1.274134564492298e+001;...
        -1.476366405440763e+005   -2.961237538429159e+002   -3.771873580263473e+001   -1.357404219788403e+001   -7.411052003232824e+000;...
        -4.686998894118387e+003   -5.145071882481552e+001   -1.151718246460839e+001   -5.524535336243413e+000   -3.611648531595958e+000;...
        -2.104710824345458e+001   -3.379418096823576e+000   -1.919928049616870e+000   -1.508399002681057e+000   -1.348510542803496e+000;...
        -1.267075422596289e-001   -2.597188113311268e-001   -4.004811495862077e-001   -5.385024279816432e-001   -6.642916520777534e-001;...
        -1.582153175255304e-001   -3.110425775503970e-001   -4.383733961816599e-001   -5.421475800719634e-001   -6.303884905318050e-001];
    stblfrac(:,6:10,1) = ...
        [-7.864009406553024e+000   -5.591791397752695e+000   -4.343949435866958e+000   -3.580521076832391e+000   -3.077683537175253e+000;...
        -4.988799898398770e+000   -3.787942909197120e+000   -3.103035515608863e+000   -2.675942594722292e+000   -2.394177022026705e+000;...
        -2.762379160216148e+000   -2.313577186902494e+000   -2.052416861482463e+000   -1.893403771865641e+000   -1.796585983161395e+000;...
        -1.284465355994317e+000   -1.267907903071982e+000   -1.279742001004255e+000   -1.309886183701422e+000   -1.349392554642457e+000;...
        -7.754208907962602e-001   -8.732998811318613e-001   -9.604322013853581e-001   -1.039287445657237e+000   -1.111986321525904e+000;...
        -7.089178961038225e-001   -7.814055112235459e-001   -8.502117698317242e-001   -9.169548634355569e-001   -9.828374636178471e-001];
    stblfrac(:,11:15,1) = ...
        [-2.729262880847457e+000   -2.479627528870857e+000   -2.297138304998905e+000   -2.162196365947914e+000   -2.061462692277420e+000;...
        -2.202290611202918e+000   -2.070075681428623e+000   -1.979193969170630e+000   -1.917168989568703e+000   -1.875099179801364e+000;...
        -1.740583121589162e+000   -1.711775396141753e+000   -1.700465158047576e+000   -1.700212465596452e+000   -1.707238269631509e+000;...
        -1.391753942957071e+000   -1.434304119387730e+000   -1.476453646904256e+000   -1.518446568503842e+000   -1.560864595722380e+000;...
        -1.180285915835185e+000   -1.245653509438976e+000   -1.309356535558631e+000   -1.372547245869795e+000   -1.436342854982504e+000;...
        -1.048835660976022e+000   -1.115815771583362e+000   -1.184614345408666e+000   -1.256100352867799e+000   -1.331235978799527e+000];
    stblfrac(:,16:20,1) = ...
        [-1.985261982958637e+000   -1.926542865732525e+000   -1.880296841910385e+000   -1.843044812063057e+000   -1.812387604873646e+000;...
        -1.846852935880107e+000   -1.828439745755405e+000   -1.817388844989596e+000   -1.812268962543248e+000   -1.812387604873646e+000;...
        -1.719534615317151e+000   -1.736176665562027e+000   -1.756931455967477e+000   -1.782079727531726e+000   -1.812387604873646e+000;...
        -1.604464355709833e+000   -1.650152416312346e+000   -1.699029550621646e+000   -1.752489822658308e+000   -1.812387604873646e+000;...
        -1.501904088536648e+000   -1.570525854475943e+000   -1.643747672313277e+000   -1.723509779436442e+000   -1.812387604873646e+000;...
        -1.411143947581252e+000   -1.497190629447853e+000   -1.591104422133556e+000   -1.695147748117837e+000   -1.812387604873646e+000];
    
    stblfrac(:,1:5,2) = ...
        [-4.738866777987500e+002   -1.684460387562537e+001   -5.619926961081743e+000   -3.281734135829228e+000   -2.397479160864619e+000;...
        -2.185953347160669e+001   -3.543320127025984e+000   -1.977029667649595e+000   -1.507632281031653e+000   -1.303310228044346e+000;...
        -2.681009914911080e-001   -4.350930213152404e-001   -5.305212880041126e-001   -6.015232065896753e-001   -6.620641788021128e-001;...
        -9.503065419472154e-002   -1.947070824738389e-001   -2.987136341021804e-001   -3.973064532664002e-001   -4.838698271554803e-001;...
        -1.264483719244014e-001   -2.437377726529247e-001   -3.333750988387906e-001   -4.016893641684894e-001   -4.577316520822721e-001;...
        -1.526287733702501e-001   -2.498255243669921e-001   -3.063859169446500e-001   -3.504924054764082e-001   -3.911254396222550e-001];
    stblfrac(:,6:10,2) = ...
        [-1.959508008521143e+000   -1.708174380583835e+000   -1.550822278332538e+000   -1.447013328833974e+000   -1.376381920471173e+000;...
        -1.199548019673933e+000   -1.144166826374866e+000   -1.115692821970145e+000   -1.103448361903579e+000   -1.101126400280696e+000;...
        -7.174026993828067e-001   -7.694003004766365e-001   -8.178267862332173e-001   -8.615585464741182e-001   -9.003104216523169e-001;...
        -5.579448431371428e-001   -6.215822273361273e-001   -6.771753949313707e-001   -7.267793058476849e-001   -7.720164852674839e-001;...
        -5.069548741156986e-001   -5.523620701546919e-001   -5.956554729327528e-001   -6.378655338388568e-001   -6.796745661620428e-001;...
        -4.309657384679277e-001   -4.709130419301468e-001   -5.113624096299824e-001   -5.525816075847192e-001   -5.948321009341774e-001];
    stblfrac(:,11:15,2) = ...
        [-1.327391983207241e+000   -1.292811209009340e+000   -1.267812588403031e+000   -1.249132310044230e+000   -1.234616432819130e+000;...
        -1.104531584444055e+000   -1.110930462397609e+000   -1.118760810700929e+000   -1.127268239360369e+000   -1.136171639806347e+000;...
        -9.347554970493899e-001   -9.658656088352816e-001   -9.945788535033495e-001   -1.021718797792234e+000   -1.048005562158225e+000;...
        -8.141486817740096e-001   -8.541760575495752e-001   -8.929234555236560e-001   -9.311104141820112e-001   -9.694099704722252e-001;...
        -7.215886443544494e-001   -7.640354693071291e-001   -8.074261467088205e-001   -8.522003643607233e-001   -8.988670244927735e-001;...
        -6.384119892912432e-001   -6.836776839822375e-001   -7.310612144698296e-001   -7.810921001396979e-001   -8.344269070778757e-001];
    stblfrac(:,16:20,2) = ...
        [-1.222879780072203e+000   -1.213041554808853e+000   -1.204541064608597e+000   -1.197016952370690e+000   -1.190232162899989e+000;...
        -1.145449097190615e+000   -1.155224344271089e+000   -1.165719407748303e+000   -1.177246763148178e+000   -1.190232162899989e+000;...
        -1.074094694885961e+000   -1.100624477495892e+000   -1.128270402039747e+000   -1.157812818875688e+000   -1.190232162899989e+000;...
        -1.008502023575024e+000   -1.049129636922346e+000   -1.092166845038550e+000   -1.138712425453996e+000   -1.190232162899989e+000;...
        -9.480479125009214e-001   -1.000533792677121e+000   -1.057363229272293e+000   -1.119941850176443e+000   -1.190232162899989e+000;...
        -8.918931068397437e-001   -9.545526172382969e-001   -1.023797332562095e+000   -1.101496412960141e+000   -1.190232162899989e+000];
    
    stblfrac(:,1:5,3) = ...
        [-1.354883142615948e+000   -8.855778500552980e-001   -7.773858277863260e-001   -7.357727812399337e-001   -7.181850957003714e-001;...
        -5.193811327974376e-002   -1.633949875159595e-001   -2.617724006156590e-001   -3.392619822712012e-001   -4.018554923458003e-001;...
        -6.335376612981386e-002   -1.297738965263227e-001   -1.985319371835911e-001   -2.624863717000360e-001   -3.174865471926985e-001;...
        -9.460338726038994e-002   -1.756165596280472e-001   -2.282691311262980e-001   -2.638458905915733e-001   -2.918110046315503e-001;...
        -1.158003423724520e-001   -1.620942232133271e-001   -1.790483132028017e-001   -1.937097725890709e-001   -2.109729530977958e-001;...
        -5.695213481951577e-002   -2.485009114767256e-002   -2.455774348005581e-002   -4.243720620421176e-002   -6.906960852184874e-002];
    stblfrac(:,6:10,3) = ...
        [ -7.120493514301658e-001   -7.121454153857569e-001   -7.157018373526386e-001   -7.209253714350538e-001   -7.265425280053609e-001;...
        -4.539746445467862e-001   -4.979328472153985e-001   -5.348184073267474e-001   -5.654705188376931e-001   -5.909430146259388e-001;...
        -3.637544360366539e-001   -4.030045272659678e-001   -4.369896090801292e-001   -4.671253359013797e-001   -4.944847533335236e-001;...
        -3.167744873288179e-001   -3.408290016876749e-001   -3.649204420006245e-001   -3.894754728525021e-001   -4.146904022890949e-001;...
        -2.311198638992638e-001   -2.537077422985343e-001   -2.783252370301364e-001   -3.047045003309861e-001   -3.327092628454751e-001;...
        -1.000745485866474e-001   -1.334091111747126e-001   -1.681287272131953e-001   -2.038409527302062e-001   -2.404547731975402e-001];
    stblfrac(:,11:15,3) = ...
        [-7.317075569303094e-001   -7.359762286696208e-001   -7.392122467978279e-001   -7.414607677550720e-001   -7.428480570989012e-001;...
        -6.123665499489599e-001   -6.307488506465194e-001   -6.469130897780404e-001   -6.615145568123281e-001   -6.750798357120451e-001;...
        -5.198770070249209e-001   -5.439265161390062e-001   -5.671356857543234e-001   -5.899325077218274e-001   -6.127077038151078e-001;...
        -4.406707089221509e-001   -4.675033009839270e-001   -4.952960990683358e-001   -5.242037261193876e-001   -5.544463409264927e-001;...
        -3.623063449447594e-001   -3.935470145089454e-001   -4.265595391976379e-001   -4.615525703717921e-001   -4.988293297210071e-001;...
        -2.780623638274261e-001   -3.168837529800063e-001   -3.572466721186688e-001   -3.995862986780706e-001   -4.444626893956575e-001];
    stblfrac(:,16:20,3) = ...
        [-7.435216571211187e-001   -7.436225251216279e-001   -7.432733099840527e-001   -7.425762029730668e-001   -7.416143171871161e-001;...
        -6.880470899358724e-001   -7.008026232247697e-001   -7.137148222421971e-001   -7.271697520465581e-001   -7.416143171871161e-001;...
        -6.358474023877762e-001   -6.597648782206755e-001   -6.849381555866478e-001   -7.119602076523737e-001   -7.416143171871161e-001;...
        -5.863313160876512e-001   -6.202819599064874e-001   -6.568811178840162e-001   -6.969403639254603e-001   -7.416143171871159e-001;...
        -5.388134824040952e-001   -5.820906647738434e-001   -6.294732446564461e-001   -6.821024214831549e-001   -7.416143171871159e-001;...
        -4.925935308416445e-001   -5.449092276644302e-001   -6.026377433551201e-001   -6.674379829825384e-001   -7.416143171871159e-001];
    
    stblfrac(:,1:5,4) = ...
        [-4.719005698760254e-003   -5.039419714218448e-002   -1.108600074872916e-001   -1.646393852283324e-001   -2.088895889525075e-001;...
        -3.167687806490741e-002   -6.488347295237770e-002   -9.913854730442322e-002   -1.306663969875579e-001   -1.574578108363950e-001;...
        -6.256908981229170e-002   -1.058190431028687e-001   -1.215669874255146e-001   -1.261149689648148e-001   -1.284283108027729e-001;...
        -7.132464704948761e-002   -5.885471032381771e-002   -3.846810486653290e-002   -2.801768649688129e-002   -2.615407079824540e-002;...
        1.186775035989228e-001    1.847231744541209e-001    1.899666578065291e-001    1.756596652192159e-001    1.538218851318199e-001;...
        1.359937191266603e+000    7.928324704017256e-001    6.068350758065271e-001    4.949176895753282e-001    4.117787224185477e-001];
    stblfrac(:,6:10,4) = ...
        [-2.445873831127209e-001   -2.729819770922066e-001   -2.951510874462016e-001   -3.121233685073350e-001   -3.249196962329062e-001;...
        -1.797875581290475e-001   -1.986122400020671e-001   -2.148458045681510e-001   -2.292024720743768e-001   -2.422125650878785e-001;...
        -1.318108373643454e-001   -1.372885008966837e-001   -1.450218673440198e-001   -1.548461140242879e-001   -1.664940537646226e-001;...
        -3.037902421859952e-002   -3.894619676380785e-002   -5.076849313651704e-002   -6.518223105549245e-002   -8.178056142331483e-002;...
        1.287679439328719e-001    1.022243387982872e-001    7.488543991005173e-002    4.698265181928261e-002    1.852002327642577e-002;...
        3.435869264264112e-001    2.844376471729288e-001    2.312306852681522e-001    1.820841981890349e-001    1.357181057787019e-001];
    stblfrac(:,11:15,4) = ...
        [-3.344714240325961e-001   -3.415532212363377e-001   -3.467713617249639e-001   -3.505859000173167e-001   -3.533413466958321e-001;...
        -2.542699931601989e-001   -2.656748454748664e-001   -2.766656461455947e-001   -2.874428940341864e-001   -2.981872822548070e-001;...
        -1.796994139325742e-001   -1.942454974557965e-001   -2.099854734361004e-001   -2.268483937252861e-001   -2.448403779828917e-001;...
        -1.003134231215546e-001   -1.206343411798188e-001   -1.426762955132322e-001   -1.664453845103147e-001   -1.920257997377931e-001;...
        -1.062008675791458e-002   -4.062891141128176e-002   -7.175196683590498e-002   -1.042870733773311e-001   -1.385948877988075e-001;...
        9.117291945474759e-002    4.766184332000264e-002    4.481886485253039e-003   -3.904933750228177e-002   -8.364689014849616e-002];
    stblfrac(:,16:20,4) = ...
        [-3.552947623689004e-001   -3.566384591258251e-001   -3.575167387322836e-001   -3.580387843935552e-001   -3.582869092425832e-001;...
        -3.090746307371333e-001   -3.202900038682522e-001   -3.320450798333745e-001   -3.445973947956370e-001   -3.582869092425832e-001;...
        -2.640470286750166e-001   -2.846415660837839e-001   -3.069024734642628e-001   -3.312464672828315e-001   -3.582869092425832e-001;...
        -2.195942670864279e-001   -2.494428999135824e-001   -2.820166786810741e-001   -3.179740384308457e-001   -3.582869092425832e-001;...
        -1.751227987938045e-001   -2.144432379167035e-001   -2.573138196343415e-001   -3.047716553689650e-001   -3.582869092425832e-001;...
        -1.301133939768983e-001   -1.794049920724848e-001   -2.327202766583559e-001   -2.916310469293936e-001   -3.582869092425832e-001];
    
    stblfrac(:,1:5,5) = ...
        [                      0                         0                         0                         0                         0;...
        -2.998229841415443e-002   -3.235136568035350e-002   -1.058934315424071e-002    1.472786013654386e-002    3.649529125352272e-002;...
        -4.911181618214269e-004    7.928758678692660e-002    1.295711243349632e-001    1.575625247967377e-001    1.726794061650541e-001;...
        6.444732609572413e-001    5.412205715497974e-001    4.864603927210872e-001    4.457073928551408e-001    4.118964225372133e-001;...
        4.884639795042095e+000    1.686842470765597e+000    1.132342494635284e+000    8.944978064032267e-001    7.538011200000044e-001;...
        2.410567057697245e+001    4.005534670805399e+000    2.144263118197206e+000    1.518214626927320e+000    1.198109338317733e+000];
    stblfrac(:,6:10,5) = ...
        [                      0                         0                         0                         0                         0;...
        5.320761222262883e-002    6.497369053185199e-002    7.235439352353751e-002    7.603800885095309e-002    7.671459793802817e-002;...
        1.799982238321182e-001    1.821699713013862e-001    1.806145618464317e-001    1.761248753943454e-001    1.691770293512301e-001;...
        3.823074983529713e-001    3.554905959697276e-001    3.305043126978712e-001    3.066571802106021e-001    2.834017043112906e-001;...
        6.558265419066330e-001    5.806408912949470e-001    5.191065509143589e-001    4.663489244354866e-001    4.194539705064985e-001;...
        9.966378800612080e-001    8.532685386168033e-001    7.427048697651345e-001    6.524693172360032e-001    5.756299950589361e-001];
    stblfrac(:,11:15,5) = ...
        [                      0                         0                         0                         0                         0;...
        7.500001602159387e-002    7.139599669434762e-002    6.628276247821394e-002    5.992932695316782e-002    5.250925428603021e-002;...
        1.600901411017374e-001    1.491003610537801e-001    1.363865273697878e-001    1.220722641614886e-001    1.062191001109524e-001;...
        2.602853501366307e-001    2.369238065872132e-001    2.129824521942899e-001    1.881563959610275e-001    1.621474808586950e-001;...
        3.765099860312678e-001    3.361566147323812e-001    2.973499640484341e-001    2.592283952427927e-001    2.210255604589869e-001;...
        5.079606300067100e-001    4.466711396792393e-001    3.897746494263863e-001    3.357416130711989e-001    2.832892169418335e-001];
    stblfrac(:,16:20,5) = ...
        [                      0                         0                         0                         0                         0;...
        4.411421669339249e-002    3.476266163507976e-002    2.439917920106283e-002    1.289010976694223e-002                         0;...
        8.881586460416716e-002    6.976629777350905e-002    4.886974404989612e-002    2.578932638717129e-002                         0;...
        1.346349888095220e-001    1.052403813710735e-001    7.348119932151805e-002    3.870673240105876e-002                         0;...
        1.820030836908522e-001    1.413881485626739e-001    9.829989964989198e-002    5.165115573609639e-002                         0;...
        2.312355801087936e-001    1.783807793433976e-001    1.233869208812706e-001    6.463145748462040e-002    9.714451465470120e-017];
    
    stblfrac(:,1:5,6) = ...
        [ 4.719005698760275e-003    5.039419714218456e-002    1.108600074872919e-001    1.646393852283322e-001    2.088895889525074e-001;...
        1.944613194060750e-001    3.117984496788369e-001    3.615078716560812e-001    3.879646155737581e-001    4.042606354602197e-001;...
        3.045958300133999e+000    1.315675725057089e+000    9.757973307352019e-001    8.294361410388060e-001    7.456405896421690e-001;...
        2.339312510820383e+001    3.858569195402605e+000    2.091507439545032e+000    1.515362821077606e+000    1.231804842218289e+000;...
        1.231812404655975e+002    9.151933726881032e+000    3.856468345925451e+000    2.470027172456050e+000    1.862167039303084e+000;...
        5.049829135345403e+002    1.890722475322573e+001    6.427275565975617e+000    3.715903402980179e+000    2.636417882085815e+000];
    stblfrac(:,6:10,6) = ...
        [ 2.445873831127209e-001    2.729819770922065e-001    2.951510874462016e-001    3.121233685073347e-001    3.249196962329060e-001;...
        4.152379986226543e-001    4.229018705591941e-001    4.280900470005300e-001    4.311273812611276e-001    4.321442286112657e-001;...
        6.900226415397631e-001    6.495436520935480e-001    6.180526887451320e-001    5.921654464012007e-001    5.697923159645174e-001;...
        1.060749495885882e+000    9.442937075476816e-001    8.583603822642385e-001    7.911221543980916e-001    7.360251815557063e-001;...
        1.521067254392224e+000    1.300039377551776e+000    1.142711537858461e+000    1.023045102736937e+000    9.273664178094935e-001;...
        2.065989355542487e+000    1.711228437455139e+000    1.466088158475343e+000    1.283765226486882e+000    1.140575450959062e+000];
    stblfrac(:,11:15,6) = ...
        [3.344714240325963e-001    3.415532212363379e-001    3.467713617249641e-001    3.505859000173170e-001    3.533413466958320e-001;...
        4.312423594533669e-001    4.285591238013830e-001    4.242644840754073e-001    4.185310514289916e-001    4.115050794489342e-001;...
        5.495326577846258e-001    5.304020801294532e-001    5.116943409858906e-001    4.928954730588648e-001    4.736165965702772e-001;...
        6.890778676134198e-001    6.476526200515113e-001    6.099033923678876e-001    5.744600864566568e-001    5.402514096915735e-001;...
        8.477633920324498e-001    7.792812067953944e-001    7.185943530039393e-001    6.633207377171386e-001    6.116407715135426e-001;...
        1.023262411940948e+000    9.237922892835746e-001    8.369566524681974e-001    7.591595457820644e-001    6.877508180861301e-001];
    stblfrac(:,16:20,6) = ...
        [ 3.552947623689000e-001    3.566384591258254e-001    3.575167387322835e-001    3.580387843935554e-001    3.582869092425831e-001;...
        4.032875933324668e-001    3.939222836649399e-001    3.833860261287606e-001    3.715758694363207e-001    3.582869092425831e-001;...
        4.535361612745278e-001    4.323485980953122e-001    4.097162006469898e-001    3.852184728042033e-001    3.582869092425835e-001;...
        5.063904595668142e-001    4.720865286037160e-001    4.365637761840112e-001    3.989743423180101e-001    3.582869092425835e-001;...
        5.620594176198462e-001    5.132627179036522e-001    4.639774715385669e-001    4.128508865888630e-001    3.582869092425835e-001;...
        6.206265009880273e-001    5.559603894356728e-001    4.919976875425384e-001    4.268552022160075e-001    3.582869092425835e-001];
    
    stblfrac(:,1:5,7) = ...
        [ 1.354883142615939e+000    8.855778500552969e-001    7.773858277863266e-001    7.357727812399328e-001    7.181850957003700e-001;...
        2.264297017396562e+001    3.703766301758638e+000    2.034998948698223e+000    1.510923485095245e+000    1.265729978744353e+000;...
        1.955956459466261e+002    1.118917023817671e+001    4.357570503031440e+000    2.718083521990130e+000    2.041945502327640e+000;...
        1.131527106972301e+003    2.742019413138009e+001    8.094356141096943e+000    4.405625422851678e+000    3.045873292912599e+000;...
        4.991370610374878e+003    5.832596523112534e+001    1.361736440227531e+001    6.617793943005997e+000    4.277065691957527e+000;...
        1.808482789458792e+004    1.120299053944505e+002    2.131886896428897e+001    9.395528700779570e+000    5.735282952993835e+000];
    stblfrac(:,6:10,7) = ...
        [ 7.120493514301658e-001    7.121454153857567e-001    7.157018373526382e-001    7.209253714350531e-001    7.265425280053608e-001;...
        1.126910935459891e+000    1.039315711942880e+000    9.801156996469297e-001    9.380990288559633e-001    9.070002633955093e-001;...
        1.682687096145072e+000    1.462088170281394e+000    1.313508264506275e+000    1.206803763884095e+000    1.126395471042167e+000;...
        2.368493556832589e+000    1.968378518204384e+000    1.704951233806636e+000    1.518043793772535e+000    1.377948007790416e+000;...
        3.176211386678905e+000    2.549432728119129e+000    2.146593646702069e+000    1.865193645178458e+000    1.656315874739094e+000;...
        4.099439855675913e+000    3.198582996879541e+000    2.632582798272859e+000    2.243339709179312e+000    1.957469852365064e+000];
    stblfrac(:,11:15,7) = ...
        [ 7.317075569303093e-001    7.359762286696208e-001    7.392122467978273e-001    7.414607677550722e-001    7.428480570989009e-001;...
        8.829463516299942e-001    8.633779161543368e-001    8.465599716104961e-001    8.313215935120923e-001    8.168794983145117e-001;...
        1.063360967480519e+000    1.012144436660489e+000    9.690437805764626e-001    9.314651792280744e-001    8.975270882378618e-001;...
        1.268363069256580e+000    1.179563109954373e+000    1.105319244270462e+000    1.041384485194864e+000    9.846979577532636e-001;...
        1.493891969504980e+000    1.362797559741365e+000    1.253624580847262e+000    1.160149469096889e+000    1.078008118654219e+000;...
        1.736744887299007e+000    1.559416515511960e+000    1.412280239489399e+000    1.286729855523644e+000    1.176933895080190e+000];
    stblfrac(:,16:20,7) = ...
        [ 7.435216571211178e-001    7.436225251216276e-001    7.432733099840527e-001    7.425762029730666e-001    7.416143171871158e-001;...
        8.027015701907034e-001    7.884022863227798e-001    7.736657968963813e-001    7.581862145381915e-001    7.416143171871158e-001;...
        8.658237613571567e-001    8.352619776464638e-001    8.049334692839693e-001    7.740056420537431e-001    7.416143171871158e-001;...
        9.329399521299938e-001    8.842632875709708e-001    8.371061471443788e-001    7.900396709438159e-001    7.416143171871157e-001;...
        1.003953952010710e+000    9.354146255148074e-001    8.702022492276336e-001    8.062927602676150e-001    7.416143171871157e-001;...
        1.078670034479511e+000    9.886802003678273e-001    9.042295460529033e-001    8.227686378257326e-001    7.416143171871157e-001];
    
    stblfrac(:,1:5,8) = ...
        [4.738866777987514e+002    1.684460387562540e+001    5.619926961081758e+000    3.281734135829232e+000    2.397479160864624e+000;...
        4.841681688643794e+003    5.491635522391771e+001    1.256979234254407e+001    6.069209132601843e+000    3.940274296039883e+000;...
        3.154616792561625e+004    1.420805372229245e+002    2.403953052063284e+001    9.998426062380954e+000    5.930362539243756e+000;...
        1.520631636586534e+005    3.148956061770992e+002    4.132943146104890e+001    1.518515134801384e+001    8.367182529059960e+000;...
        5.901656732159231e+005    6.246491282963873e+002    6.581680474603525e+001    2.173557079848703e+001    1.125045444319795e+001;...
        1.944624278667431e+006    1.139848804168331e+003    9.894809619823921e+001    2.974824391888133e+001    1.458002371721213e+001];
    stblfrac(:,6:10,8) = ...
        [1.959508008521145e+000    1.708174380583837e+000    1.550822278332539e+000    1.447013328833976e+000    1.376381920471174e+000;...
        2.963447020215305e+000    2.423693540860402e+000    2.089182215079736e+000    1.865572849084425e+000    1.708118159360888e+000;...
        4.190132768594454e+000    3.268280841745006e+000    2.710662024401290e+000    2.341995909523891e+000    2.082469140437107e+000;...
        5.624308785058203e+000    4.226708866462347e+000    3.402197103627229e+000    2.865360079281767e+000    2.490393899977397e+000;...
        7.254212029229660e+000    5.287806421003054e+000    4.154585933912857e+000    3.428194997160839e+000    2.925780747207696e+000;...
        9.070365685373144e+000    6.442950257298201e+000    4.960971490178073e+000    4.025088868546689e+000    3.384287797654701e+000];
    stblfrac(:,11:15,8) = ...
        [ 1.327391983207241e+000    1.292811209009341e+000    1.267812588403031e+000    1.249132310044230e+000    1.234616432819130e+000;...
        1.593041126030172e+000    1.506471132927683e+000    1.439628954887186e+000    1.386580264484466e+000    1.343153406231364e+000;...
        1.891158929781140e+000    1.745070641877115e+000    1.630251730907927e+000    1.537630629971792e+000    1.460938380853296e+000;...
        2.214464603850502e+000    2.003098342270666e+000    1.835905829230373e+000    1.700021765831942e+000    1.586823477367793e+000;...
        2.557944985263177e+000    2.276562749626175e+000    2.053593165082403e+000    1.871725504345519e+000    1.719630879614922e+000;...
        2.918103805585008e+000    2.562588803694463e+000    2.281050180010934e+000    2.051085944176459e+000    1.858294826115218e+000];
    stblfrac(:,16:20,8) = ...
        [ 1.222879780072204e+000    1.213041554808854e+000    1.204541064608597e+000    1.197016952370690e+000    1.190232162899990e+000;...
        1.306371038922589e+000    1.274091491606534e+000    1.244744203398707e+000    1.217124809801410e+000    1.190232162899990e+000;...
        1.395630981221581e+000    1.338301797693731e+000    1.286320343916442e+000    1.237570697847646e+000    1.190232162899990e+000;...
        1.490188322141933e+000    1.405530485165501e+000    1.329245194088195e+000    1.258353899045780e+000    1.190232162899990e+000;...
        1.589489775546923e+000    1.475587597649461e+000    1.373481210080780e+000    1.279472666002594e+000    1.190232162899990e+000;...
        1.692973560150181e+000    1.548256386823049e+000    1.418980226656540e+000    1.300924242481222e+000    1.190232162899990e+000];
    
    stblfrac(:,1:5,9) = ...
        [1.890857122067037e+006    1.074884919696010e+003    9.039223076384690e+001    2.645987890965103e+001    1.274134564492299e+001;...
        1.434546473316804e+007    2.987011338973518e+003    1.804473474220022e+002    4.487048929338575e+001    1.960113433547389e+001;...
        7.716266115204613e+007    6.969521346220721e+003    3.196657990381036e+002    6.941784107578008e+001    2.798990029407097e+001;...
        3.253192550565641e+008    1.437315176424486e+004    5.205876769957880e+002    1.006582035946658e+002    3.790739646062081e+001;...
        1.143638705833100e+009    2.703823367877713e+004    7.964291266167923e+002    1.391051003571698e+002    4.935349274736288e+001;...
        3.492208269966229e+009    4.737075925045248e+004    1.161019167208514e+003    1.852377745522907e+002    6.232811767701676e+001];
    stblfrac(:,6:10,9) = ...
        [7.864009406553027e+000    5.591791397752693e+000    4.343949435866960e+000    3.580521076832391e+000    3.077683537175252e+000;...
        1.132727408868559e+001    7.671280872680232e+000    5.732691330034323e+000    4.573075545294608e+000    3.818589092027862e+000;...
        1.533578393202605e+001    9.991349773961725e+000    7.243609507849516e+000    5.634462725204553e+000    4.601857009791827e+000;...
        1.985701175129152e+001    1.252691966593449e+001    8.859346059355138e+000    6.752431092162364e+000    5.418366793527828e+000;...
        2.486500490402286e+001    1.525895955988075e+001    1.056731639206889e+001    7.918478700695184e+000    6.262067266019560e+000;...
        3.033836510475647e+001    1.817240938152932e+001    1.235792736188858e+001    9.126360342186048e+000    7.128676006881803e+000];
    stblfrac(:,11:15,9) = ...
        [2.729262880847459e+000    2.479627528870858e+000    2.297138304998906e+000    2.162196365947915e+000    2.061462692277420e+000;...
        3.297130126832188e+000    2.920640582387343e+000    2.640274592919582e+000    2.426998377788287e+000    2.262233765245289e+000;...
        3.893417077901593e+000    3.382471282597797e+000    2.999860062957988e+000    2.705234908082859e+000    2.473610569743775e+000;...
        4.510891038980249e+000    3.859051363710381e+000    3.370702720510665e+000    2.992693808833481e+000    2.692636527934335e+000;...
        5.144949915764652e+000    4.346635348399592e+000    3.749599176843221e+000    3.286641675099088e+000    2.917178603817272e+000;...
        5.792462636377325e+000    4.842756648977701e+000    4.134472567050430e+000    3.585273662390985e+000    3.145733197974777e+000];
    stblfrac(:,16:20,9) = ...
        [1.985261982958638e+000    1.926542865732524e+000    1.880296841910385e+000    1.843044812063057e+000    1.812387604873647e+000;...
        2.133064562958712e+000    2.029912595114798e+000    1.945516531961286e+000    1.874392545595589e+000    1.812387604873647e+000;...
        2.288441176274372e+000    2.137883347336651e+000    2.012884307837858e+000    1.906295529437326e+000    1.812387604873647e+000;...
        2.449737610939970e+000    2.249772716121334e+000    2.082221357100924e+000    1.938735806854783e+000    1.812387604873647e+000;...
        2.615585030563546e+000    2.364937633815368e+000    2.153342270485199e+000    1.971693892149562e+000    1.812387604873647e+000;...
        2.784907129216124e+000    2.482804054400846e+000    2.226062706102394e+000    2.005149380181030e+000    1.812387604873647e+000];
    
    % Interpolate to find initial value of stableinv, set 0 values in case of extrapolation
    betaL0 = betaCol<0 & middle;
    if any(betaL0)
        utemp(betaL0) = 1-utemp(betaL0);
    end
    X0(middle) = sign(betaCol(middle)).*interp3(Alp,Bet,P,stblfrac,alphaCol(middle), abs(betaCol(middle)), utemp(middle), 'linear',0);
end

X0 = reshape(X0,size(u));
end

function r = stablernd(alpha,beta,gam,delta,varargin)
%STABLERND Random arrays from the stable distribution.
%
% References:
%       [1] A. Weron and R. Weron (1995), "Computer Simulation of Levy
%       alpha-Stable Variables and Processes", Lecture Notes in Physics,
%       Springer-Verlag.
%       [2] J.P. Nolan (2015), "Stable Distributions - Models for Heavy
%       Tailed Data", Birkhauser, Boston. In progress, Chapter 1 online at
%       academic2.american.edu/~jpnolan

[err, sizeOut, numelOut] = internal.stats.statsizechk(4, alpha, beta, gam, delta, varargin{:});
if err > 0
    error(message('stats:addstable:InputSizeMismatch'));
end

if numelOut ~= 1
    if isscalar(alpha), alpha = repmat(alpha,sizeOut); end
    if isscalar(beta), beta = repmat(beta,sizeOut); end
    if isscalar(gam), gam = repmat(gam,sizeOut); end
    if isscalar(delta), delta = repmat(delta,sizeOut); end
end

% Return NaN for illegal parameter values.
alpha(alpha<=0 | alpha>2) = NaN;
beta(beta<-1 | beta>1) = NaN;
gam(gam <= 0) = NaN;

% Preallocation
outType = internal.stats.dominantType(alpha, beta, gam, delta);
r = zeros(sizeOut,"like",outType);

% Generate the STANDARD stable distribution, gam = 1 and delta = 0.
%   Use corresponding method for special cases (Normal and Cauchy),
%   otherwise use general method.
Gau = alpha == 2;
if any(Gau(:)) % Gaussian distribution
    r(Gau) = sqrt(2) * randn([sum(Gau(:)),1],"like",r);
end

Cau = (alpha == 1) & (beta == 0);
if any(Cau(:))
    % Cauchy(0,1) is t-distribution with d.f. = 1
    r(Cau) = trnd(1,[sum(Cau(:)),1]);
end

Levy = (alpha == 0.5) & (abs(beta)==1);
if any(Levy(:)) % Levy distribution
    % Levy(0,1) = S(0.5,1,1,1;0), convert it to S(0.5,1,1,0;0) by subtracting -1
    r(Levy) = 1./ randn([sum(Levy(:)),1],"like",r).^2-1;
    r(Levy) = beta(Levy).*r(Levy);
end

other = (~Cau) & (~Gau) & (~Levy);
if any(other(:))
    U = pi.*(rand(sizeOut,"like",r)-0.5); % Uniform distribution on (-pi/2,pi/2)
    W = -log(rand(sizeOut,"like",r)); % Exponential distribution with mean 1
    
    alphaE1 = (alpha == 1) & other;
    alphaN1 = (alpha ~= 1) & other;
    if any(alphaE1(:))
        r(alphaE1) = (2/pi) .* ( (pi/2 + beta(alphaE1) .*U(alphaE1)).*tan(U(alphaE1)) - ...
            beta(alphaE1).*log( (pi/2.*W(alphaE1).*cos(U(alphaE1)))./(pi/2 + beta(alphaE1).*U(alphaE1)) ));
    end
    if any(alphaN1(:))
        Btan = beta(alphaN1).*tan(pi*alpha(alphaN1)/2);
        aUB = alpha(alphaN1).*(U(alphaN1) + atan(Btan)./alpha(alphaN1));
        r(alphaN1) = (1 + Btan.^2).^(1./(2*alpha(alphaN1))).*...
            sin( aUB ) ./ (cos(U(alphaN1)).^(1./alpha(alphaN1))) .* ...
            ( cos(U(alphaN1)-aUB) ./ W(alphaN1) ).^((1-alpha(alphaN1))./alpha(alphaN1));
        % shift when alpha~=1
        r(alphaN1) = r(alphaN1) - Btan;
    end
end

% Scale and shift from the STANDARD stable distribution
r = gam.*r + delta;

end

function [nlogL,acov] = stablelike(params, data)
%STABLELIKE Negative log-likelihood for the stable distribution.
% NLOGL and ACOV are computed using interpolated densities.

% Reference:
%   [1] J. P. Nolan (2001) "Maximum Likelihood Estimation and Diagnostics
%       for Stable Distributions". Lévy processes, Birkhäuser Boston, 379-400.

n = numel(data);
alpha = params(1);
beta = params(2);
gam = params(3);
delta = params(4);
z = (data-delta)./gam;
nlogL = sum(neglog_pdf(z,alpha,beta)) + n*log(gam);

atBoundary = false;
if abs(2-alpha)<eps^(1/3)
    warning(message('stats:addstable:ConvergedToBoundary1'));
    atBoundary = true;
end

if abs(abs(beta)-1)<eps^(1/3)
    warning(message('stats:addstable:ConvergedToBoundary2'));
    atBoundary = true;
end

if nargout > 1
    if ~atBoundary
        % Get the integrand in the Fisher matrix, F, as defined in [1] by
        %  F_{i,j} = \int_{-\infty}^{infty} df/d\theta_i * df/d\theta_j * 1/f dx
        F = zeros(4,4);
        step = eps^(1/4);
        
        % There might be a bunch warnings in the integration, turn them off temporarily
        ws1 = warning('off','MATLAB:integral:MaxIntervalCountReached');
        ws2 = warning('off','MATLAB:integral:NonFiniteValue');
        ws3 = warning('off','MATLAB:integral:MinStepSize');
        
        for j = 1:4
            for i = 1:j
                % Default accuracy of INTEGRAL is sometimes difficult to
                % reach, because the pdf values of tails in the pdf table
                % are from fewer grid points than those in the middle. Even
                % at this accuracy setting, it may sometimes have warnings
                % from the integration procedure.
                F(i,j) = integral(@(x)infoMtxCal(x,params,step,i,j),...
                    -Inf,Inf,'AbsTol',1e-6,'RelTol',1e-4);
            end
        end
        
        % Restore the original warning state
        cleanupObj = onCleanup(@()warning(ws1));
        cleanupObj = onCleanup(@()warning(ws2));
        cleanupObj = onCleanup(@()warning(ws3));
        
        % Find the asymptotic covariance matrix.
        F = triu(F,1)' + F;
        %Make sure the Hessian is pos def, refuse to compute the cov matrix if not.
        [~,s] = chol(F);
        if s > 0
            warning(message('stats:addstable:NonPosDefHessian'));
            acov = NaN(4);
        else
            acov = inv(F)./n;
        end
    else
        %If the shape parameter estimates are at boundary, cov is not reliable.
        acov = NaN(4);
    end
end
end


function integ = infoMtxCal(x, params, step, i ,j)
% INFOMTXCAL Calculate the Fisher information matrix using tabulated pdf.

theta1 = params;
theta1(i) = theta1(i) - step;
theta2 = params;
theta2(i) = theta2(i) + step;
dfthetai = (tabpdf((x-theta2(4))/theta2(3),theta2(1),theta2(2))./theta2(3) -...
    tabpdf((x-theta1(4))/theta1(3),theta1(1),theta1(2))./theta1(3))/2/step;
theta3 = params;
theta3(j) = theta3(j) - step;
theta4 = params;
theta4(j) = theta4(j) + step;
dfthetaj = (tabpdf((x-theta4(4))/theta4(3),theta4(1),theta4(2))./theta4(3) -...
    tabpdf((x-theta3(4))/theta3(3),theta3(1),theta3(2))./theta3(3))/2/step;
f = tabpdf((x-params(4))/params(3),params(1),params(2))./params(3);
integ = dfthetai.*dfthetaj./f;

% It is possible that variable x after standardization is out of the
% range of x in the pdf table, which causes extrapolation resulting 0 pdf
% values (or sometimes the pdf values itself is 0). This can make the
% quotient above NaN's or Inf's. If it happens, the integrand at those
% points are redefined as zero.
integ(f==0) = 0;

end


%===========Stable Fitting Functions===========

function [parmhat,parmci] = stablefit(x,alpha,options)
%STABLEFIT Parameter estimates and confidence intervals for stable data.
% This method is only suitable for the shape parameter ALPHA >= 0.4.

% Reference:
%   [1] J. P. Nolan (2001) "Maximum likelihood estimation and diagnostics for
%       stable distributions." Lévy processes. Birkhäuser Boston, p379-400.

if ~isvector(x)
    error(message('stats:addstable:VectorRequired'));
end

if nargin < 2 || isempty(alpha)
    alpha = 0.05;
end

% Compute initial estimates by McCulloch quantile method to speed up the
% optimization procedure.
params0 = intMle(x);

% If the initial estimates of ALPHA or BETA are near or at the boundary,
% fminsearch may not make any update on the parameter estimates. So, try to
% reassign some reasonable closer values to the initial estimates of ALPHA
% and BETA.
if (2-params0(1) < 1e-3)
    params0(1) = 1.95;
end
if (1-abs(params0(2)) < 1e-3)
    params0(2)= sign(params0(2))*0.95;
end

% Transform the parameters to eliminate bounds in order to apply fminsearch
phi0 = varTrans(params0,'forward');

% The default options include turning fminsearch's display off.  This
% function gives its own warning/error messages, and the caller can turn
% display on to get the text output from fminsearch if desired.
if nargin < 3 || isempty(options)
    options = statset('stablefit');
else
    options = statset(statset('stablefit'),options);
end

% Maximize the log-likelihood with respect to the transformed parameters
% When possible, use 'fmincon' for speed
if (license('test', 'Optimization_Toolbox'))
    [parmhat,~,err,output] = fmincon(@(params)stable_nloglf(x,params),phi0, ...
        [],[],[],[],[],[],[],options);
else
    [parmhat,~,err,output] = fminsearch(@(params)stable_nloglf(x,params),phi0,options);
end

if (err == 0)
    % fminsearch may print its own output text; in any case give something
    % more statistical here, controllable via warning IDs.
    if output.funcCount >= options.MaxFunEvals
        warning(message('stats:addstable:EvalLimit'));
    else
        warning(message('stats:addstable:IterLimit'));
    end
elseif (err < 0)
    error(message('stats:addstable:NoSolution'));
end

% Transform the estimated parameters back
parmhat = varTrans(parmhat,'backward');

% Check if parameters are valid
try
    checkargs(parmhat(1),parmhat(2),parmhat(3),parmhat(4));
catch ME
    error(message('stats:addstable:FitError',ME.message));
end

if nargout > 1
    [~,acov] = stablelike(parmhat,x);
    probs = [alpha/2; 1-alpha/2];
    se = sqrt(diag(acov))';
    parmci = norminv([probs probs probs probs], [parmhat; parmhat], [se; se]);
end
end

function params = intMle(x)
%INTMLE Parameter estimation of stable data by McCulloch's quantile method.
% This method is only suitable for ALPHA>0.5.

% Reference:
%   [1] J. H. McCulloch (1986) "Simple Consistent Estimators of Stable
%       Distribution Parameters". Cummun. Stat. -Simula., 15(4).

if ~isvector(x)
    error(message('stats:addstable:VectorRequired'));
end

Xpcts = prctile(x,[95,75,50,25,5]);
nuAlpha = min(25,max((Xpcts(1) - Xpcts(5))/(Xpcts(2) - Xpcts(4)),2.439));
nuBeta = min(1,max((Xpcts(1) + Xpcts(5) - 2*Xpcts(3))/(Xpcts(1) - Xpcts(5)),-1));

% Input tables for Alpha and Beta estimation
nuA = [2.439 2.5 2.6 2.7 2.8 3.0 3.2 3.5 4.0 5.0 6.0 8.0 10 15 25];
nuB = [0 0.1 0.2 0.3 0.5 0.7 1];
[na, nb] = meshgrid( nuA , nuB );
alphaTbl=  [2.000 2.000 2.000 2.000 2.000 2.000 2.000;...
    1.916 1.924 1.924 1.924 1.924 1.924 1.924;...
    1.808 1.813 1.829 1.829 1.829 1.829 1.829;...
    1.729 1.730 1.737 1.745 1.745 1.745 1.745;...
    1.664 1.663 1.663 1.668 1.676 1.676 1.676;...
    1.563 1.560 1.553 1.548 1.547 1.547 1.547;...
    1.484 1.480 1.471 1.460 1.448 1.438 1.438;...
    1.391 1.386 1.378 1.364 1.337 1.318 1.318;...
    1.279 1.273 1.266 1.250 1.210 1.184 1.150;...
    1.128 1.121 1.114 1.101 1.067 1.027 0.973;...
    1.029 1.021 1.014 1.004 0.974 0.935 0.874;...
    0.896 0.892 0.887 0.883 0.855 0.823 0.769;...
    0.818 0.812 0.806 0.801 0.780 0.756 0.691;...
    0.698 0.695 0.692 0.689 0.676 0.656 0.595;...
    0.593 0.590 0.588 0.586 0.579 0.563 0.513]';
betaTbl=  [ 0.000 2.160 1.000 1.000 1.000 1.000 1.000;...
    0.000 1.592 3.390 1.000 1.000 1.000 1.000;...
    0.000 0.759 1.800 1.000 1.000 1.000 1.000;...
    0.000 0.482 1.048 1.694 1.000 1.000 1.000;...
    0.000 0.360 0.760 1.232 2.229 1.000 1.000;...
    0.000 0.253 0.518 0.823 1.575 1.000 1.000;...
    0.000 0.203 0.410 0.632 1.244 1.906 1.000;...
    0.000 0.165 0.332 0.499 0.943 1.560 1.000;...
    0.000 0.136 0.271 0.404 0.689 1.230 2.195;...
    0.000 0.109 0.216 0.323 0.539 0.827 1.917;...
    0.000 0.096 0.190 0.284 0.472 0.693 1.759;...
    0.000 0.082 0.163 0.243 0.412 0.601 1.596;...
    0.000 0.074 0.147 0.220 0.377 0.546 1.482;...
    0.000 0.064 0.128 0.191 0.330 0.478 1.362;...
    0.000 0.056 0.112 0.167 0.285 0.428 1.274]';

Alpha = interp2(na,nb,alphaTbl,nuAlpha,abs(nuBeta));

Beta = sign(nuBeta) * interp2(na,nb,betaTbl,nuAlpha,abs(nuBeta));

% Reset Alpha if necessary.
if Alpha > 2
    Alpha = 2;
    Beta = sign(nuBeta);
end

% Reset Beta if necessary.
if Beta > 1
    Beta = 1;
elseif Beta < -1
    Beta = -1;
end

% Input tables for Gam and Delta estimation
va = [2 1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 1.1 1.0 0.9 0.8 0.7 0.6 0.5];
vb = [0 0.25 0.5 0.75 1];
gamTbl = [  1.908 1.908 1.908 1.908 1.908;...
    1.914 1.915 1.916 1.918 1.921;...
    1.921 1.922 1.927 1.936 1.947;...
    1.927 1.930 1.943 1.961 1.987;...
    1.933 1.940 1.962 1.997 2.043;...
    1.939 1.952 1.988 2.045 2.116;...
    1.946 1.967 2.022 2.106 2.211;...
    1.955 1.984 2.067 2.188 2.333;...
    1.965 2.007 2.125 2.294 2.491;...
    1.980 2.040 2.205 2.435 2.696;...
    2.000 2.085 2.311 2.624 2.973;...
    2.040 2.149 2.461 2.886 3.356;...
    2.098 2.244 2.676 3.265 3.912;...
    2.189 2.392 3.004 3.844 4.775;...
    2.337 2.635 3.542 4.808 6.247;...
    2.588 3.037 4.534 6.636 9.144]';
deltaTbl = [0.000  0.000  0.000  0.000  0.000;...
    0.000 -0.017 -0.032 -0.049 -0.064;...
    0.000 -0.030 -0.061 -0.092 -0.123;...
    0.000 -0.043 -0.088 -0.132 -0.179;...
    0.000 -0.056 -0.111 -0.170 -0.232;...
    0.000 -0.066 -0.134 -0.206 -0.283;...
    0.000 -0.075 -0.154 -0.241 -0.335;...
    0.000 -0.084 -0.173 -0.276 -0.390;...
    0.000 -0.090 -0.192 -0.310 -0.447;...
    0.000 -0.095 -0.208 -0.346 -0.508;...
    0.000 -0.098 -0.223 -0.383 -0.576;...
    0.000 -0.099 -0.237 -0.424 -0.652;...
    0.000 -0.096 -0.250 -0.469 -0.742;...
    0.000 -0.089 -0.262 -0.520 -0.853;...
    0.000 -0.078 -0.272 -0.581 -0.997;...
    0.000 -0.061 -0.279 -0.659 -1.198]';
[ma, mb] = meshgrid( va, vb );
nuGam = interp2(ma,mb,gamTbl,Alpha,abs(Beta));
Gam = (Xpcts(2)-Xpcts(4))/nuGam;

% Reset the value of Gam if necessary
if Gam==0 || isnan(Gam) || isinf(Gam)
    s = std(x);
    if s>0
        Gam = s;
    else
        Gam = eps;
    end
end

nuKsi = sign(Beta) * interp2(ma,mb,deltaTbl,Alpha,abs(Beta));
Ksi = Xpcts(3) + Gam*nuKsi;

if Alpha == 1
    Delta = Ksi;
else
    Delta = Ksi - Beta*Gam*tan(pi*Alpha/2);
end
params(1) = Alpha;
params(2) = Beta;
params(3) = Gam;
params(4) = Delta;
end

function y = tabpdf(x,alpha,beta)
%TABPDF Density for standardized stable data by interpolating tabulated density table.
x = atan(x);

% Use cubic interpolation and enforce no extrapolation (this makes
% extrapolation NaN's)
% Keep the griddedInterpolant 'persistent' for speed with repeated function
% calls
persistent G;
if (isempty(G))
    s = load('private/StablePdfTable.mat');
    G = griddedInterpolant({s.b, s.a, s.xgd}, s.p, 'linear','none');
end

y = G({beta,alpha,x});
y = reshape(y,size(x));

% Assign 0 to NaN's caused by extrapolation.
y(isnan(y)) = 0;

end

function nll = neglog_pdf(x,alpha,beta)
%NEGLOG_PDF Negative log-likelihood for standardized data using interpolated densities.

nll = -log(max(realmin,tabpdf(x,alpha,beta))); % In case of small or negative densities
end


function nll = stable_nloglf(x,params)
%STABLE_NLOGLF Objective function for stable maximum likelihood.
n = numel(x);
params = varTrans(params,'backward');
alpha = params(1);
beta = params(2);
gam = params(3);
delta = params(4);
z = (x-delta)./gam;
nll = sum(neglog_pdf(z,alpha,beta)) + n*log(gam);
end

function phi = varTrans(theta,direct)
%VARTRANS Transform the variables to apply the optimization of the log-likelihood function.
%
% When DIRECTION = 'FORWARD', THETA is a vector of parameter
% values in the stable distribution and PHI is a vector of transformed
% parameter values.
% When DIRECTION = 'BACKWARD', THETA is a vector of transformed parameter
% and PHI is a vector of the parameter values used in stable distribution.
%
% The optimization algorithm fminsearch is for unconstrained variables. The
% variables THETA=(ALPHA,BETA,GAM,DELTA) in stable distribution are constrained
% such that 0<ALPHA<=2, -1<=BETA<=1, GAM>0 and -Inf<DELTA<Inf. To eliminate the
% bounds of THETA, we can do the following transformation:
%       ALPHA = a + (b - a)./(1 + exp(-PHI(1))), where a=0, b=2
%       BETA = a + (b - a)./(1 + exp(-PHI(2))), where a=-1, b=1
%       GAM = exp(PHI(3))
%       DELTA = PHI(4)
% then the parameters PHI are unbounded.

narginchk(2,Inf);
phi = ones(size(theta));

if strcmpi(direct,'forward')
    phi(1) = -log(2./theta(1) - 1);
    phi(2) = -log((2./(theta(2)+1) - 1));
    phi(3) = log(theta(3));
    phi(4) = theta(4);
else % direct = backward
    phi(1) = 2./(1 + exp(-theta(1)));
    phi(2) = -1 + 2./(1 + exp(-theta(2)));
    phi(3) = exp(theta(3));
    phi(4) = theta(4);
end
end



