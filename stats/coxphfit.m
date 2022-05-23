function [b,logL,H,stats] = coxphfit(X,y,varargin)
%COXPHFIT Cox proportional hazards regression.
%   B=COXPHFIT(X,Y) fits Cox's proportional hazards regression model to the
%   data in the vector or 2-column matrix Y, using the columns of the
%   matrix X as predictors. Y is typically a vector providing the event
%   time of right censored time-to-event data. X and Y must have the same
%   number of rows, and X should not contain a column of ones. The result B
%   is a vector of coefficient estimates. This model states that the hazard
%   rate for the distribution of Y can be written as h(t)*exp(X*B) where
%   h(t) is a common baseline hazard function.
%
%   For data with time-dependent covariates, Y is specified as a 2-column
%   matrix. Each individual appears in multiple rows, with one row for each
%   of the distinct values of the time-dependent variables. The
%   corresponding row of Y provides the (start,stop] time interval. This is
%   the counting process form as described in reference [1]. For example:
%
%          Age   Med            Time0 Time1
%     X = [ 20 false       Y = [ 0     5       Censor = [0
%           30  true             0    12                 1
%           30 false]           12    20]                1]
%
%   This represents two patients. The age 20 patient was not taking
%   medication for the interval from week 0 to 5. The age 30 patient was
%   taking medication for the interval from week 0 to 12, but not for the
%   interval from week 12 to 20. The value of the 'Censoring' parameter
%   (see below) is 0 if the event (death) was observed at the end of the
%   interval, and 1 if it was not. In this example, Med is considered to be
%   a time dependent covariate while Age is not.
%
%   [...] = COXPHFIT(X,Y,'PARAM1',VALUE1,'PARAM2',VALUE2,...) specifies
%   additional parameter name/value pairs chosen from the following:
%
%      Name          Value
%      'B0'          A vector containing initial values for the estimated
%                    coefficients B. Default is 0.01/SD where SD represents
%                    the standard deviation of each predictor in X.
%      'Baseline'    The X values at which the baseline hazard is to be
%                    computed.  Default is mean(X), so the hazard at X is
%                    h(t)*exp((X-mean(X))*B).  Enter 0 to compute the
%                    baseline relative to 0, so the hazard at X is
%                    h(t)*exp(X*B).
%      'Censoring'   A boolean array of the same size as Y that is 1 for
%                    observations that are right-censored and 0 for
%                    observations that are observed exactly.  Default is
%                    all observations observed exactly.
%      'Frequency'   An array of the same size as Y containing non-negative
%                    values. FREQUENCY typically contains integer frequencies
%                    for the corresponding observations, but may contain any 
%                    non-integer values which represent observation weights. 
%                    Default is 1 per row of X and Y.
%      'Options'     A structure specifying control parameters for the
%                    iterative algorithm used to estimate B.  This argument
%                    can be created by a call to STATSET.  For parameter
%                    names and default values, type STATSET('coxphfit').
%      'Strata'      An array of the same number of rows as Y containing the 
%                    stratification variables. The values of the stratification 
%                    variables must be real numbers. Observations are grouped
%                    into strata based on the different combinations of values
%                    in the stratification variables. Each stratum defines 
%                    a set of observations having a common baseline hazard.
%      'Ties'        A string specifying the method to handle tied failure 
%                    times. The methods are equivalent if there are no ties.
%                    Two methods are available:
%                      'breslow' Breslow's approximation to the partial  
%                                likelihood. It is the default.
%                      'efron'   Efron's approximation to the partial 
%                                likelihood.
%
%   [B,LOGL,H,STATS]=COXPHFIT(...) returns additional results.  LOGL
%   is the log likelihood. For unstratified model, H is a two-column matrix 
%   containing y values in column 1 and the estimated baseline cumulative 
%   hazard evaluated at those values in column 2. For stratified model, H
%   is a (2+K)-column matrix with the last K columns being the K 
%   stratification variables in the same order as in STRATA array. STATS is
%   a structure containing the following fields:
%       'beta'         coefficient estimates (same as B output)
%       'se'           standard errors of coefficient estimates
%       'z'            z statistics for B (B divided by standard error)
%       'p'            p-values for B
%       'covb'         estimated covariance matrix for B
%       'csres'        Cox-Snell residuals
%       'devres'       deviance residuals
%       'martres'      martingale residuals
%       'schres'       Schoenfeld residuals
%       'sschres'      scaled Schoenfeld residuals
%       'scores'       score residuals
%       'sscores'      scaled score residuals
%
%   For Cox-Snell, martingale and deviance residuals, the result is a column 
%   vector with one row per observation. For (scaled) Schoenfeld and (scaled)
%   score residuals, the result is a matrix of the same size as X. (scaled) 
%   Schoenfeld residuals of censored observations are defined to be NaNs. 
%   Residuals obtained from models with time-dependent covariates have the
%   same number of rows as Y.
%
%   Example:
%       % Generate Weibull data with A depending on predictor X
%       x = 4*rand(100,1); A = 50*exp(-0.5*x); B = 2;
%       y = wblrnd(A,B);
%    
%       % Fit Cox model
%       [b,logL,H,st] = coxphfit(x,y);
%    
%       % Show Cox estimate of baseline survivor and known Weibull function
%       figure;
%       stairs(H(:,1),exp(-H(:,2)))
%       xx = linspace(0,100);
%       line(xx,1-wblcdf(xx,50*exp(-0.5*mean(x)),B),'color','r')
%       title(sprintf('Baseline survivor function for X=%g',mean(x)));
%
%       % Assess model adequacy by residual plots
%       load(fullfile(matlabroot,'examples','stats','readmissiontimes.mat'))
%       [b1,logL1,H1,st1] = coxphfit([Age Sex Weight Smoker],ReadmissionTime,...
%                           'censoring',Censored,'ties','efron');
%
%       % Draw a deviance residual plot
%       Deviance = st1.devres;
%       SmokerDev = Deviance(Smoker==1);
%       NonsmokerDev = Deviance(Smoker==0);
%       figure;
%       plot((1:length(SmokerDev)),SmokerDev,'ro',...
%            (1:length(NonsmokerDev)),NonsmokerDev,'b*')
%       % Add a horizontal reference line with intercept 0
%       hline = refline([0,0]);
%       hline.Color = 'k';
%       legend('Smoker','Non-smoker','Location','Best')
%       xlabel('Observation number')
%       ylabel('Deviance residuals')
%
%       % Draw a martingale residual plot
%       Schres = st1.schres(:,4);      
%       figure;
%       plot(ReadmissionTime(Smoker==1),Schres(Smoker==1),'ro',...
%            ReadmissionTime(Smoker==0),Schres(Smoker==0),'b*')
%       hline = refline([0,0]);
%       hline.Color = 'k';
%       legend('Smoker','Non-smoker','Location','Best')
%       xlabel('Time')
%       ylabel('Schoenfeld residuals')
%       
%       % Refit the model by stratifying a covariate
%       [b2,logL2,H2,st2] = coxphfit([Age Sex Weight], ReadmissionTime,...
%                           'censoring',Censored,'ties','efron','strata',Smoker);
%
%   See also ECDF, STATSET, WBLFIT.
%
% References:
%   [1] Therneau, Terry M., and Patricia M. Grambsch, Modeling survival 
%       data: extending the Cox model. Springer Science & Business Media, 2000.

%   Copyright 2005-2018 The MathWorks, Inc.


if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

narginchk(2,inf);
% Check the required data arguments
if ndims(X)>2 || ~isreal(X)
    error(message('stats:coxphfit:BadX'));
end
if ~(isvector(y) || size(y,2)==2) || ~isreal(y)
    error(message('stats:coxphfit:BadY'));
end
if size(y,2)==2
    if any(y(:,1) >= y(:,2))
        error(message('stats:coxphfit:BadY2'));
    end
end

% Process the optional arguments
okargs =   {'baseline' 'censoring' 'frequency' {'init' 'b0'}  'strata' 'ties'     'options' };
defaults = {[]         []          []           []            []       'breslow'  []};
[baseX, cens, freq, init, strata, ties, options] = ...
                    internal.stats.parseArgs(okargs,defaults,varargin{:});

if ~isempty(cens) && (~isvector(cens) || ~all(ismember(cens,0:1)))
    error(message('stats:coxphfit:BadCensoring'));
end
if ~isempty(freq) && (~isvector(freq) || ~isreal(freq) || any(freq<0))
    error(message('stats:coxphfit:BadFrequency'));
end
if ~isempty(baseX) && ~(isnumeric(baseX) && (isscalar(baseX) || ...
                             (isvector(baseX) && length(baseX)==size(X,2))))
    error(message('stats:coxphfit:BadBaseline'));
elseif isscalar(baseX)
    baseX = repmat(baseX,1,size(X,2));
end
strataFlag = 0; % Default is no strata
if ~isempty(strata) 
    strataFlag = 1;
    if ~isreal(strata) 
        error(message('stats:coxphfit:BadStrata'));
    end
end
ties = lower(ties);
if isempty(ties)
    ties = 'breslow';
elseif ~ischar(ties)
    error(message('stats:coxphfit:BadTies'));
elseif ~ismember(ties,{'breslow', 'efron'})
    error(message('stats:coxphfit:BadTiesName',ties));
end

% Make sure the inputs agree in size, and remove NaNs
freq(freq==0) = NaN;   % easiest way to deal with zero-frequency observations
[badin,~,y,X,cens,freq,strata]=statremovenan(y,X,cens,freq,strata);
if badin>0
    whichbad = {'Y' 'X' '''censoring''' '''frequency''' '''strata'''};
    error(message('stats:coxphfit:InputSizeMismatch', whichbad{ badin }));
end

% Default values of cens and frequency
[n,p] = size(X);
if isempty(cens)
    cens = false(n,1);
end
if isempty(freq)
    freq = ones(n,1);
end

% Data pre-processing before fitting
if ~strataFlag
    [X,sorty,cens,freq,atrisk,tied,lastloc,tiedygrps,idx] = ...
                                        dataProcessing(X,y,cens,freq,ties);
else    
    ndeath = sum(~cens);
    nStrataVar = size(strata,2);
    
    % Find groups according to distinct combinations of values in strata
    [group,~,~,~,ngroups] = mgrp2idx(strata,n);
    groups = cell(1,ngroups);
    for gnum = 1:ngroups
        groups{gnum} = find(group==gnum);
    end
    
    % Pre-process data within each stratum
    grpdata = cell(ngroups,1); % Store preprocessed data in grpdata
    for j = 1:ngroups
        gX = X(groups{j},:);
        gy = y(groups{j},:);
        gcens = cens(groups{j},:);
        gfreq = freq(groups{j},:);
        gstrata = strata(groups{j},:);
        [gX,gsorty,gcens,gfreq,gatrisk,gtied,glastloc,gtiedgrps,gyidx] = ...
                                    dataProcessing(gX,gy,gcens,gfreq,ties);
        gX = Xcentering(gX,gfreq,baseX);
        % Save grouped data
        grpdata{j} = {gX,gfreq,gcens,gatrisk,gtied,gsorty,gtiedgrps,...
                      glastloc,gstrata,gyidx};
    end
end

X = Xcentering(X,freq,baseX);
% Try to diagnose some potential problems
try   
    if rank([ones(n,1), X]) < (p+1)
        if n>1 && any(all(diff(X,1)==0,1))
            warning(message('stats:coxphfit:DisallowedContantTerm'));
        else
            warning(message('stats:coxphfit:RankDeficient'));
        end
    end
catch ME
    if isequal(ME.identifier,'MATLAB:svd:matrixWithNaNInf')
        warning(message('stats:coxphfit:XorFreqWithInf'));
    end
end

% Get starting values of coefficients b
sumf = max(1, sum(freq));
if isempty(init)
    stdX = sqrt((freq'*X.^2)/sumf)';
    b0 = zeros(size(stdX),class(X));
    t = (stdX ~= 0);
    b0(t) = 0.01 ./ stdX(t);
elseif isvector(init) && numel(init)==p && isreal(p) && isnumeric(p)
    b0 = init(:);
else
    error(message('stats:coxphfit:BadInit', p));
end

% Perform the fit
unboundwng = false;
if p==0
    b = zeros(0,1,class(X));
else
    % The default options include turning statsfminbx's display off.  This
    % function gives its own warning/error messages, and the caller can
    % turn display on to get the text output from statsfminbx if desired.
    options = statset(statset('coxphfit'),options);
    options = optimset(options);
    dflts = struct('DerivativeCheck','off', 'HessMult',[], ...
        'HessPattern',ones(p), 'PrecondBandWidth',Inf, ...
        'TypicalX',ones(p,1), 'MaxPCGIter',1, 'TolPCG',0.1);

    % Maximize the log-likelihood
    if ~strataFlag
        funfcn = {'fungradhess' 'coxphfit' @(b)negloglike(b) [] []};
    else
        funfcn = {'fungradhess' 'coxphfit' @(b)negloglikeStrata(b) [] []}; 
    end
    lastb = [];
    try
        [b, ~, ~, err, output] = ...
            statsfminbx(funfcn, b0, [], [], options, dflts, 1);
        if (err == 0)
            % Check for convergence failure
            if output.funcCount >= options.MaxFunEvals
                warning(message('stats:coxphfit:EvalLimit'));
            else
                warning(message('stats:coxphfit:IterLimit'));
            end
            unboundwng = true;
        elseif (err < 0)
            error(message('stats:coxphfit:NoSolution'));
        end

    catch ME
        if isequal(ME.identifier,'stats:statsfminbx:BadFunctionOutput') || ...
           isequal(ME.identifier,'MATLAB:eig:matrixWithNaNInf')
            warning(message('stats:coxphfit:EvalWarning'));
             b = lastb;
             unboundwng = true;
        else
            m = message('stats:coxphfit:FitError');
            throw(addCause(MException(m.Identifier,'%s',getString(m)),ME));
        end
    end
end

% Compute log likelihood at the solution, incl. baseline hazard
if ~strataFlag
    [logL,~,ddl,H,mlflag,resid] = LogL(X,b,freq,cens,atrisk,lastloc,tied,...
                                    sorty,ties,tiedygrps,idx);

    % Try to diagnose a likelihood with no local maximum
    if ~unboundwng && mlflag
        warning(message('stats:coxphfit:FitWarning'));
    end

    if nargout>=4
        % Compute standard errors, etc.
        [covb,se,z,pval] = statscov(b,ddl);
        stats.covb = covb;
        stats.beta = b;
        stats.se = se;
        stats.z = z;
        stats.p = pval;
        
        stats.csres = resid.csres;
        stats.devres = resid.devres;
        stats.martres = resid.martres;
        stats.schres = resid.schres;
        stats.sschres = sum(~cens)*stats.schres*covb + repmat(b',size(cens));
        stats.scores = resid.scores;
        stats.sscores = stats.scores*covb;
    end
else 
    [logL,~,ddl,H,res] = LogLstrata(ngroups,grpdata,b,ties,p,nStrataVar,ndeath,groups);
    
    % Compute standard errors, etc.
    if nargout>=4
        [covb,se,z,pval] = statscov(b,ddl);
        stats.covb = covb;
        stats.beta = b;
        stats.se = se;
        stats.z = z;
        stats.p = pval;

        stats.csres = res.csres;
        stats.devres = res.devres;
        stats.martres = res.martres;
        stats.schres = res.schres;
        stats.sschres = sum(~cens)*stats.schres*covb + repmat(b',size(cens));
        stats.scores = res.scores;
        stats.sscores = stats.scores*covb;
    end   
end

% -------------------------------------------
function [L,dl,ddl] = negloglike(b)
    % Compute negative log likelihood
    [L,dl,ddl] = LogL(X,b,freq,cens,atrisk,lastloc,tied,sorty,ties,tiedygrps,idx);
    L = -L;
    dl = -dl;
    ddl = -ddl;
    lastb = b;
end

function [L,dl,ddl] = negloglikeStrata(b)
    % Compute negative log likelihood for stratified model
    [L,dl,ddl] = LogLstrata(ngroups,grpdata,b,ties,p,nStrataVar,ndeath,groups);
    L = -L;
    dl = -dl;
    ddl = -ddl;
    lastb = b;
end
end

function [X,sorty,cens,freq,atrisk,tied,lastloc,tiedygrps,idx] = dataProcessing(X,y,cens,freq,ties)
    % Sort by increasing time
    if size(y,2)==1
        [sorty,idx] = sort(y);
    else
        [~,idx] = sort(y(:,2));
        sorty = y(idx,:);
    end
    X = X(idx,:);
    cens = cens(idx);
    freq = freq(idx);
    
    % Determine the observations at risk at each time
    if size(y,2)==1
        [~,atrisk] = ismember(sorty,flipud(sorty), 'legacy');
        atrisk = length(sorty) + 1 - atrisk;
        % Determine the last location of distinct values in sorty 
        [~,lastloc] = ismember(sorty,sorty, 'legacy');
        tied = diff(sorty) == 0;
    else
        atrisk = grp2idx(sorty(:,2));
        [~,lastloc] = ismember(atrisk, atrisk, 'legacy');
        tied = diff(sorty(:,2)) == 0;
    end
    tied = [false;tied] | [tied;false];
    
    % Find observation indices of tied noncensored sorty values
    if strcmp(ties,'breslow')
        tiedygrps = [];
    else
        uy = unique(sorty(:,end));
        ngrp = length(uy);
        tiedygrps = cell(ngrp,1);
        for i = 1:ngrp
            tiedygrps{i} = find(sorty(:,end)==uy(i) & ~cens);
        end
    end
end

function X = Xcentering(X,freq,baseX)
    % Recenter X to make it better conditioned; does not affect coefficients
    n = size(X,1);
    sumf = max(1, sum(freq));
    if ~isempty(X)
        if isempty(baseX)
            baseX = (freq'*X) / sumf;
        end
        X = X - repmat(baseX,n,1);
    end
end

function cums = timeDependentRiskset(sorty,Q,atrisk)
    %Compute quantities at risk for Cox model with time varying covariates.
    % sorty is the 2-column sorted matrix whose first and second columns represent
    % the starting and ending points of intervals in the counting process syntax.
    % Q is a matrix representing the quantity that need to be calculated for the
    % risk sets. atrisk is a column vector of the index indicating
    % observations at risk.

    uystop = unique(sorty(:,2));
    len = length(uystop);
    cums = zeros(len,size(Q,2));
    for k = 1:len
        cums(k,:) = sum( Q((uystop(k)>sorty(:,1) & (uystop(k)<=sorty(:,2))),:),1 );
    end
    cums = cums(atrisk,:);
end

function H = cumhazest(obsfreq, risksum, sorty)            
    %Compute estimate of baseline cumulative hazard
    terms = obsfreq ./ risksum;
    cumbh = cumsum(terms);
    H = [sorty(:,end), cumbh];
    H = H(obsfreq>0,:);
    if ~isempty(H)
        H = [H(1,1), 0; H];
    end
end

function [covb,se,z,p] = statscov(b,ddl)
    % Compute cov, standard error, z score and p value
    
    covb = cast(-inv(ddl),class(b));
    varb = diag(covb);
    % In cases where we did not converge to a maximum, the 2nd derivative
    % matrix may not be negative definite.  It's likely a warning was
    % displayed and that varb has -Inf values.  We'll give infinite
    % variances for parameters affected by this.
    se = Inf(size(varb),class(b));
    se(varb>0) = sqrt(varb(varb>0));
    % se may be empty, we need a column vector
    se = se(:);  
    z = b ./ se;
    p = 2*normcdf(-abs(z));
end

function [L,dl,ddl,H,res] = LogLstrata(ngroups,grpdata,b,ties,p,nStrataVar,ndeath,groupidx)
    % Compute log-likelihood, first derivative and second derivative of each stratum.
    L = 0;
    dl = zeros(1,p);
    ddl = zeros(p,p);
    if nargout>=4
        % Used to store cumulative hazards of each stratum; the baseline
        % hazards are along the sorted y within each stratum. For each 
        % stratum, the baseline hazards starts from a zero. 
        len = 0;
        H = zeros(ngroups+ndeath,2+nStrataVar);   
    end

    for j = 1:ngroups 
        gX = grpdata{j}{1};
        gfreq = grpdata{j}{2};
        gcens = grpdata{j}{3};
        gatrisk = grpdata{j}{4};
        gtied = grpdata{j}{5};
        gsorty = grpdata{j}{6};
        gtiedgrps = grpdata{j}{7};
        glastloc = grpdata{j}{8};
        gstrata = grpdata{j}{9};
        gyidx = grpdata{j}{10};
        [gL,gdl,gddl,gH,~,gresid] = LogL(gX,b,gfreq,gcens,gatrisk,glastloc,...
                                        gtied,gsorty,ties,gtiedgrps,gyidx);
        if nargout>=4
            glen = size(gH,1);
            H(len+1:len+glen,:) = [gH repmat(gstrata(1,:),glen,1)];
            len = len+glen;
        end
        if nargout>=5
            res.csres(groupidx{j},:) = gresid.csres;
            res.devres(groupidx{j},:) = gresid.devres;
            res.martres(groupidx{j},:) = gresid.martres;
            res.scores(groupidx{j},:) = gresid.scores;
            res.schres(groupidx{j},:) = gresid.schres;          
        end
        L = L + gL;
        dl = dl + gdl;
        ddl = ddl + gddl;
    end
end

function [L,dl,ddl,H,mlflag,res]=LogL(X,b,freq,cens,atrisk,lastloc,tied,sorty,ties,tiedygrps,yidx)
    % Compute log likelihood L
    obsfreq = freq .* ~cens;
    Xb = X*b;
    r = exp(Xb);
    rfull = freq.*r;
    if size(sorty,2)==1
        risksum = cumsum(rfull, 'reverse');
        risksum = risksum(atrisk);
    else
        risksum = timeDependentRiskset(sorty,rfull,atrisk);
    end
    
    if strcmpi(ties,'breslow')
        L = obsfreq'*(Xb - log(risksum));
        
        % Compute first derivative dL/db
        [n,p] = size(X);
        Xr = X .* repmat(rfull,1,p);
        if size(sorty,2)==1
            Xrsum = cumsum(Xr, 'reverse');
            Xrsum = Xrsum(atrisk,:);
        else
            Xrsum = timeDependentRiskset(sorty,Xr,atrisk);
        end
        A = Xrsum ./ repmat(risksum,1,p);
        dl = obsfreq' * (X-A);
        if nargout>=5
            % Set the mlflag (monotone likelihood flag) to indicate if the
            % likelihood appears to be monotone, not at an optimum.  This
            % can happen if, at each of the sorted failure times, the
            % specified linear combination of the X's is larger than that
            % of all other observations at risk.
            if n>2
                mlflag = all(cens(1:end-1) | (diff(Xb)<0 & ~tied(1:end-1)));
            else
                mlflag = true;
            end
        end
        
        % Compute second derivative d2L/db2
        t1 = repmat(1:p,1,p);
        t2 = sort(t1);
        XXr = X(:,t1) .* X(:,t2) .* repmat(r.*freq,1,p^2);
        if size(sorty,2)==1
            XXrsum = cumsum(XXr, 1, 'reverse');
            XXrsum = XXrsum(atrisk,:,:) ./ repmat(risksum,[1,p^2]);
        else
            XXrsum = timeDependentRiskset(sorty,XXr,atrisk)./repmat(risksum,[1,p^2]);
        end
        ddl = reshape(-obsfreq'*XXrsum, [p,p]);
        ddl = ddl + A'*(A.*repmat(obsfreq,1,p));

        if nargout>=4
            % Compute estimate of baseline cumulative hazard
            H = cumhazest(obsfreq, risksum, sorty);
        end
        if nargout>=6
            % Baseline hazards are defined at distinct failure times: for
            % tied data set, they share the value of the last obs of each 
            % tied group; thus, to calculate residuals, set all obs of each
            % tied group to have that value(the last one of each tied group).
            cumbh = cumsum(obsfreq ./ risksum);
            cumbh = cumbh(lastloc);
            
            % Cox-Snell residuals
            csres = cumbh.*r;
            res.csres(yidx,:) = csres;
            
            % Martingale residuals
            if size(sorty,2)==1
                martres = (~cens) - csres;
            else
                tempres = internal.stats.coxmresmex(~cens',sorty(:,1)',sorty(:,2)',freq',r',0);
                martres = (~cens) - tempres';
            end
            res.martres(yidx,:) = martres;
            
            % Deviance residuals
            devres = sign(martres).*sqrt(-2*(martres + (~cens).*log((~cens)-martres)));
            % When martingale residual is zero, deviance residual is defined as zero
            devres(martres==0) = 0;
            res.devres(yidx,:) = devres;
            
            % Schoenfeld residuals
            schres = (~cens) .* (X-A);
            schrescopy = schres;
            % Schoenfeld residuals are undefined at censored times, set them to be NaNs 
            schres(cens==1,:) = NaN;
            res.schres(yidx,:) = schres;
            
            % Score residuals
            if size(X,1)>1
                diffH = [cumbh(1);cumbh(2:end)-cumbh(1:end-1)];
            else
                diffH = cumbh;
            end
            if size(sorty,2)==1
                res.scores(yidx,:) = schrescopy - X.*r.*cumbh + r.*cumsum(A.*diffH);
            else
                tempres = internal.stats.coxsresmex(~cens', X', sorty(:,1)', sorty(:,2)', freq', r' ,0);
                res.scores(yidx,:) = schrescopy - tempres';
            end
        end
    else        
        % Efron's method        
        [n,p] = size(X);
        fr = freq.*r;
        Xr = X .* repmat(fr,1,p);
        if size(sorty,2)==1
            Xrsum = cumsum(Xr, 'reverse');
            Xrsum = Xrsum(atrisk,:);
        else
            Xrsum = timeDependentRiskset(sorty,Xr,atrisk);
        end
        A = Xrsum ./ repmat(risksum,1,p);
        
        t1 = repmat(1:p,1,p);
        t2 = sort(t1);
        XXr = X(:,t1) .* X(:,t2) .* repmat(fr,1,p^2);
        if size(sorty,2)==1
            XXrsum = cumsum(XXr, 1, 'reverse');
            XXrsum = XXrsum(atrisk,:,:);
        else
            XXrsum = timeDependentRiskset(sorty,XXr,atrisk);
        end
        B = XXrsum ./ repmat(risksum,[1,p^2]);
        
        wtrisksum = risksum; % Weighted version of risksum
        wtfreq = freq; % Weighted version of freq
        
        % Average within tied groups
        for j = 1:length(tiedygrps)
            idxDeath = tiedygrps{j};
            len = numel(idxDeath);
            wtfreq(idxDeath(:)) = sum(freq(idxDeath(:)))./len;
            if len>1
                idxDeath = idxDeath(:);          
                avgrisksum = risksum(idxDeath(2:end)) - (1:len-1)'.*...
                             sum(fr(idxDeath(:)))./len;             
                wtrisksum(idxDeath(2:end)) = avgrisksum;
                A(idxDeath(2:end),:) = ( Xrsum(idxDeath(2:end),:) -...
                        (1:len-1)'.*sum(Xr(idxDeath(:),:))./len )./avgrisksum;
                B(idxDeath(2:end),:) = ( XXrsum(idxDeath(2:end),:) - ...
                        (1:len-1)'.*sum(XXr(idxDeath(:),:))./len )./avgrisksum;          
            end
        end
        
        % Compute log-likelihood, first derivative, and second derivative
        wtobsfreq = wtfreq.*(~cens);
        L = obsfreq'* Xb - wtobsfreq'*log(wtrisksum) ;
        dl = obsfreq'*X - wtobsfreq'*A;
        ddl = reshape(-wtobsfreq'*B, [p,p]) + A'*(A.*repmat(wtobsfreq,1,p));
        
        if nargout>=4
            % Compute estimate of baseline cumulative hazard
            H = cumhazest(obsfreq, risksum, sorty);
        end
        
        if nargout>=5
            % Set the mlflag (monotone likelihood flag) to indicate if the
            % likelihood appears to be monotone, not at an optimum.  This
            % can happen if, at each of the sorted failure times, the
            % specified linear combination of the X's is larger than that
            % of all other observations at risk.
            if n>2
                mlflag = all(cens(1:end-1) | (diff(Xb)<0 & ~tied(1:end-1)));
            else
                mlflag = true;
            end
        end
        
        if nargout>=6
            % Residuals of Efron's method
            % Cox-Snell residuals
            csterms = cumsum( obsfreq./risksum );
            csres = csterms(lastloc).*r;
            res.csres(yidx,:) = csres;

            % For martingale, Schoenfeld, and score residuals, reset values
            % at tied times
            terms = wtobsfreq ./ wtrisksum;
            cumbh = cumsum(terms);
            cumbhmart = cumbh(lastloc);
            hazard = wtfreq./wtrisksum;
            cumAh = cumsum(terms.*A);
            cumAhsco = cumAh(lastloc,:);
            for j = 1:length(tiedygrps)
                idxDeath = tiedygrps{j};
                len = numel(idxDeath);
                if len>1
                    % For martingale residuals: the baseline hazards are 
                    % averaged within each tied group, and the cumulative 
                    % baseline hazards have the same value as the last one 
                    % of that group
                    idxDeath = idxDeath(:);
                    frac = (1:len-1)'./len;
                    tiedbh = wtfreq(idxDeath(2:end)).*(1-frac)./wtrisksum(idxDeath(2:end));
                    tiedcumbh = cumsum([cumbh(idxDeath(1));tiedbh]);
                    cumbhmart(idxDeath) = repmat(tiedcumbh(end),len,1);
                    
                    % For Schoenfeld residuals
                    tiedA = ( Xrsum(idxDeath(2:end),:)-frac.*sum(Xr(idxDeath,:)) )./...
                                ( risksum(idxDeath(2:end)) - frac.*sum(fr(idxDeath(:))) );
                    tiedcumA = cumsum([A(idxDeath(1),:);tiedA]);
                    A(idxDeath,:) = repmat(tiedcumA(end,:),len,1)./len;
                    
                    % For score residuals
                    tiedcumAh = cumsum([cumAh(idxDeath(1),:);
                                        tiedA.*(1-frac).*hazard(idxDeath(2:end))]);
                    cumAhsco(idxDeath,:) = repmat(tiedcumAh(end,:),len,1);
                end
            end
            
            % Martingale residuals
            if size(sorty,2)==1
                resid = cumbhmart.*r;
                resid = (~cens) - resid;
            else
                resid = internal.stats.coxmresmex(~cens',sorty(:,1)',sorty(:,2)',freq',r',1);
                resid = (~cens) - resid';
            end
            res.martres(yidx,:) = resid;
            
            % Deviance residuals
            devres = sign(resid).*sqrt(-2*(resid + (~cens).*log((~cens)-resid)));
            devres(resid==0) = 0;
            res.devres(yidx,:) = devres;
            
            % Schoenfeld residuals
            schres = (~cens) .* (X-A);
            schrescopy = schres;
            schres(cens==1,:) = NaN;
            res.schres(yidx,:) = schres;
            
            % Score residuals 
            if size(sorty,2)==1
                res.scores(yidx,:) = schrescopy - X.*r.*cumbhmart + r.*cumAhsco;
            else
                tempres = internal.stats.coxsresmex(~cens', X', sorty(:,1)', sorty(:,2)', freq', r', 1);
                res.scores(yidx,:) = schrescopy - tempres';
            end
        end
    end

    if nargout>1 && any(isnan(dl(:)))
        % If we ran in to problems in computing the derivatives and got a
        % NaN result, this will cause problems in statsfminbx.  Best to set
        % the objective function to infinity, causing statsfminbx to back
        % up.
        L = -Inf;
    end
end

