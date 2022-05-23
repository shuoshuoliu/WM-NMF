function [settings, X] = cordexch(nfactors,nruns,model,varargin)
%CORDEXCH D-Optimal design of experiments - coordinate exchange algorithm.
%   [SETTINGS, X] = CORDEXCH(NFACTORS,NRUNS,MODEL) generates a D-optimal
%   design having NRUNS runs for NFACTORS factors.  SETTINGS is the
%   matrix of factor settings for the design, and X is the matrix of
%   term values (often called the design matrix).  MODEL is an optional
%   argument that controls the order of the regression model.
%   MODEL can be any of the following strings:
%
%     'linear'        constant and linear terms (the default)
%     'interaction'   includes constant, linear, and cross product terms.
%     'quadratic'     interactions plus squared terms.
%     'purequadratic' includes constant, linear and squared terms.
%
%   Alternatively MODEL can be a matrix of term definitions as
%   accepted by the X2FX function.
%
%   [SETTINGS, X] = CORDEXCH(...,'PARAM1',VALUE1,'PARAM2',VALUE2,...)
%   provides more control over the design generation through a set of
%   parameter/value pairs.  Valid parameters are the following:
%
%      Parameter     Value
%      'bounds'      Lower and upper bounds for each factor, specified
%                    as a 2-by-NFACTORS matrix.  Alternatively, this value
%                    can be a cell array containing NFACTORS elements, each
%                    element specifying the vector of allowable values for
%                    the corresponding factor.
%      'categorical' Indices of categorical predictors.
%      'display'     Either 'on' or 'off' to control display of
%                    iteration number. (default = 'on' unless 'UseParallel'
%                    is TRUE, in which case default = 'off').
%      'excludefun'  Function to exclude undesirable runs.
%      'init'        Initial design as an NRUNS-by-NFACTORS matrix
%                    (default is a randomly selected set of points).
%      'levels'      Vector of number of levels for each factor. Not used
%                    when 'bounds' is specified as a cell array.
%      'maxiter'     Maximum number of iterations (default = 10).
%      'tries'       Number of times to try to generate a design from a
%                    new starting point, using random points for each
%                    try except possibly the first (default 1). 
%      'options'     A structure that contains options specifying whether to
%                    compute multiple tries in parallel, and specifying how
%                    to use random numbers when generating the starting points
%                    for the tries. This argument can be created by a call to 
%                    STATSET. CORDEXCH uses the following fields:
%                        'UseParallel'
%                        'UseSubstreams'
%                        'Streams'
%                    For information on these fields see PARALLELSTATS.
%                    NOTE: If 'UseParallel' is TRUE and 'UseSubstreams' is FALSE,
%                    then the length of 'Streams' must equal the number of workers 
%                    used by CORDEXCH.  If a parallel pool is already open, this 
%                    will be the size of the parallel pool.  If a parallel pool 
%                    is not already open, then MATLAB may try to open a pool for 
%                    you (depending on your installation and preferences).
%                    To ensure more predictable results, it is best to use 
%                    the PARPOOL command and explicitly create a parallel pool 
%                    prior to invoking CORDEXCH with 'UseParallel' set to TRUE. 
%
%   The CORDEXCH function searches for a D-optimal design using a
%   coordinate exchange algorithm.  It creates a starting design, and then
%   iterates by changing each coordinate of each design point in an attempt
%   to reduce the variance of the coefficients that would be estimated
%   using this design.
%
%   If the 'excludefcn' function is F, it must support the syntax B=F(S) 
%   where S is a matrix of K-by-NFACTORS columns containing settings,
%   and B is a vector of K boolean values.  B(j) is true if the jth row
%   of S should be excluded.
%
%   Example:
%      % Design for two factors, quadratic model
%      sortrows(cordexch(2,9,'q'))
%
%      % Design for 2 of 3 factors making up a mixture, where factor
%      % values are up to 50%, and the two factors must not make up
%      % less than 15% or greater than 85% of the whole mixture
%      f = @(x) sum(x,2)>85 | sum(x,2)<15;
%      bnds = [0 0;50 50];
%      x=sortrows(cordexch(2,9,'q','bounds',bnds,'levels',101,'excl',f))
%      plot(x(:,1),x(:,2),'bo')
%
%   See also ROWEXCH, DAUGMENT, DCOVARY, X2FX, STATSET, PARALLELSTATS,
%   PARPOOL, PARFOR, GCP.

%   The undocumented parameters 'start' and 'covariates' are used when
%   this function is called by the daugment and dcovary functions to
%   implement D-optimal designs with fixed rows or columns.
   
%   Copyright 1993-2013 The MathWorks, Inc.

if nargin > 2
    model = convertStringsToChars(model);
end

if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if ~(isnumeric(nfactors) && isscalar(nfactors) && isreal(nfactors) && nfactors >= 0 ...
        && isequal(nfactors,floor(nfactors)) && isfinite(nfactors))
    error(message('stats:cordexch:BadNumFactors'));
end
if ~(isnumeric(nruns) && isscalar(nruns) && isreal(nruns) && nruns >= 0 ...
        && isequal(nruns,floor(nruns)) && isfinite(nruns))
    error(message('stats:cordexch:BadNumRuns'));
end

paropt = statset('cordexch');
[paropt,~,varargin] = internal.stats.parseArgs({'Options'}, {paropt}, varargin{:});

[useParallel, RNGscheme, poolsz] = ...
    internal.stats.parallel.processParallelAndStreamOptions(paropt, true);

varargin = [varargin {'UseParallel' useParallel}]; 

usePool = useParallel && poolsz>0;

[startdes,covariates,dodisp,xinit,maxiter,  ...
 tries,bnds,nlevels,excluder,categ,~] = ...
                  doptargcheck({},nfactors,nruns, varargin{:});
              
% Create design with starting rows, covariates, and zeros
nobs = size(startdes,1);
settings = zeros(nruns, nfactors);
if ~isempty(xinit)
   if ~isempty(excluder) && any(excluder(xinit))
       warning(message('stats:cordexch:InitExcluderViolated'));
   end
end

% Add initial set of rows never to change, if any
if ~isempty(startdes)
   settings = [startdes; settings];
end

% Add fixed covariates, if any
if ~isempty(covariates)
   settings = [settings covariates];
end

if nargin == 2 || isempty(model)
   model = 'linear';
end
if ~ischar(model)
   if size(settings,2) ~= size(model,2)
      error(message('stats:cordexch:InputSizeMismatch'));
   end
   modelorder = max(model(:));   % max order of any factor
else
   modelorder = 2;               % max of named models that we provide
end

if isempty(nlevels)
   nlevels = (modelorder+1) * ones(1,nfactors);
   nlevels(categ) = 2;
end
catlevels = nlevels(categ);
iscat = ismember(1:nfactors, categ);

% Convert from actual categorical factor levels to 1:maxlevel
if ~isempty(xinit)
    xinit = levels2numbers(xinit,bnds,categ,nlevels);
end

% Create an array nxij containing all possible values for each factor
% if the number of levels is the same, or a cell array nxijcell if
% different factors have different numbers of levels
nxij = [];
nxijcell = {};
if all(nlevels==nlevels(1))
   usecell = false;
   nxij = zeros(nlevels(1),nfactors);
   for j=1:nfactors
      if iscell(bnds)
          bndsj = bnds{j};
      else
          bndsj = bnds(:,j);
      end
      if iscat(j)
          nxij(:,j) = (1:nlevels(j))';
      elseif iscell(bnds) && length(bndsj)>=nlevels(j)
          nxij(:,j) = bndsj;
      else
          nxij(:,j) = linspace(bndsj(1),bndsj(end),nlevels(j))';
      end
   end
   rowlist = zeros(nlevels(1),1);
else
   usecell = true;
   nxijcell = cell(nfactors,1);
   for j=1:nfactors
      if iscell(bnds)
          bndsj = bnds{j};
      else
          bndsj = bnds(:,j);
      end
      if iscat(j)
          nxijcell{j} = (1:nlevels(j))';
      elseif iscell(bnds) && length(bndsj)>=nlevels(j)
          nxijcell{j} = bndsj;
      else
          nxijcell{j} = linspace(bndsj(1),bndsj(end),nlevels(j))';
      end
   end
end

% Change exclusion function if necessary to deal with factor levels as numbers
if ~isempty(categ) && ~isempty(excluder)
    excluder = @(x) excluder(numbers2levels(x,bnds,categ,nlevels));
end

% Create Iteration Counter Figure.
quiet = isequal(dodisp,'off');
if ~quiet
    if useParallel
        warning(message('stats:cordexch:parallelDisplay'));
        if usePool
            internal.stats.parallel.distributeToPool( ...
                'workerID', num2cell(1:poolsz) );
        end
    else
        [f, settry, setiter] = statdoptdisplay(getString(message('stats:cordexch:CoordExch')));
    end
end

% Generate designs, pick best one
setstart = settings;
warnedexcl = false;
warnedstart = false;

% Must perform first try separate from parfor loop to do error checking
% and to gather some information.
  
if isempty(xinit) 
    [xinit,warnedexcl] = internal.stats.parallel.smartForSliceout( ...
        1,@firstinit,useParallel,RNGscheme);
end
setstart(nobs+1:end,1:nfactors) = xinit;

% Generate the term values
[Xstart,model,termstart] = x2fx(setstart,model,categ,catlevels);
 
% First time only, do error checking and gather some information
totdf = size(Xstart,2);
if totdf>(size(Xstart,1)+nobs)
    if ~quiet && ~useParallel
        close(f);
    end
    error(message('stats:cordexch:TooFewRuns'));
end
modelorder = max(model,[],1);   % max order of each factor
if any(~iscat & nlevels<modelorder(1:nfactors)+1)
    if ~quiet && ~useParallel
        close(f);
    end
    error(message('stats:cordexch:TooFewLevels'));
end

temp = zeros(1,totdf);
temp(termstart) = 1;
termfactors = cumsum(temp);

factorterms = false(nfactors,totdf);
for ROW=1:nfactors
   factorterms(ROW,:) = ismember(termfactors,find(model(:,ROW)>0));
end
 
% Generate a new design from here by coordinate exchange
% (this nested function uses factorterms and some other variables)
[bestset,X,logdetX] = internal.stats.parallel.smartForSliceout( ...
    1, @firstgen1, useParallel, RNGscheme);

if tries>1
    % Setup for the return values from smartForReduce, which are:
    reductionArgument.fh = @reductionOperator;
    reductionArgument.iv = cell(3,1);
    reductionArgument.iv{1} = warnedexcl;
    reductionArgument.iv{2} = {logdetX, 1, bestset, X};
    
    reductionVar = internal.stats.parallel.smartForReduce( ...
        tries-1, @loopBody, useParallel, RNGscheme, reductionArgument);

    warnedexcl = reductionVar{1};
    bestset    = reductionVar{2}{3};
    X          = reductionVar{2}{4};
end

if warnedexcl
    warning(message('stats:cordexch:ConstraintsViolated'));
end

% Convert to actual categorical factor levels
settings = numbers2levels(bestset,bnds,categ,nlevels);

if ~quiet && ~useParallel
    close(f);
end

%   ----------------------------------
%          Nested functions
%   ----------------------------------

%   --------------------------------------------------------
    function [settings,X,logdetX] = gen1design(settings,X,S)
        
        % This is a nested function that shares some variables with its caller
        
        [~,R]=qr(X,0);
        
        % Adjust starting design if it is rank deficient, because the algorithm
        % will not proceed otherwise.
        if rank(R)<size(R,2)
            if ~warnedstart
                warning(message('stats:cordexch:BadStartingDesign'));
                warnedstart = true;
            end
            R = adjustr(R);
            wasbad = 1;
        else
            wasbad = 0;
        end
        logdetX = 2*sum(log(abs(diag(R))));
        
        iter = 0;
        madeswitch = 1;
        dcutoff = 1 + sqrt(eps);
        
        while madeswitch > 0 && iter < maxiter
            madeswitch = 0;
            iter = iter + 1;
            
            % Update iteration counter.
            if ~quiet && ~useParallel
                setiter(iter);
            end
            
            %Loop over rows of factor settings matrix.
            [~,collist] = sort(rand(S,1,nfactors));
            for row = (nobs+1):(nobs+nruns)
                fx = X(row,:);
                E = [];
                %Loop over columns of factor settings matrix.
                for col = collist
                    if usecell
                        newset = nxijcell{col};
                        xnew = repmat(settings(row,:),numel(newset),1);
                        xnew(:,col) = newset;
                    else
                        rowlist(:) = row;
                        xnew = settings(rowlist,:);
                        xnew(:,col) = nxij(:,col);
                    end
                    
                    if ~isempty(excluder)
                        excluderows = excluder(xnew);
                        xnew = xnew(~excluderows,:);
                        if isempty(xnew)
                            continue
                        end
                    end
                    
                    % Update affected terms only
                    t = model(:,col)>0;
                    fxnew = fx(ones(size(xnew,1),1),:);
                    fxnew(:,factorterms(col,:)) = x2fx(xnew,model(t,:),categ,catlevels);
                    
                    % Compute change in determinant.
                    if isempty(E)
                        E = fx/R;
                        dxold = E*E';
                    end
                    F = fxnew/R;
                    dxnew = sum(F.*F,2);
                    dxno  = F*E';
                    
                    d = (1 + dxnew).*(1 - dxold) + dxno.^2;
                    
                    % Find the maximum change in the determinant, switch if >1
                    [d,idx] = max(d);
                    if d>dcutoff || (iter==1 && ~iscat(col))
                        madeswitch = 1;
                        logdetX = log(d) + logdetX;
                        settings(row,col) = xnew(idx,col);
                        X(row,:) = fxnew(idx,:);
                        fx = X(row,:);
                        [~,R] = qr(X,0);
                        if wasbad
                            if rank(R)<size(R,2)
                                R = adjustr(R);
                            else
                                wasbad = 0;
                            end
                            logdetX = 2*sum(log(abs(diag(R))));
                        end
                        E = [];
                    end
                end
            end
        end
        
    end %-gen1design

%   --------------------------------------------
    function [xinit,warnedexcl] = firstinit(~,S)
        % Optionally, print first "try" number and worker to screen.
        if ~quiet
            if useParallel
                if usePool
                    labindx = internal.stats.parallel.workerGetValue('workerID');
                else
                    labindx = 1;
                end
                disp(getString(message('stats:cordexch:WorkerTry', labindx, 1)));
            else
                settry(1);
            end
        end
        
        if isempty(S)
            S = RandStream.getGlobalStream;
        end
        warnedexcl = false;
        [xinit,warnedexcl] = randstart(bnds, nlevels, nruns, nfactors, ...
            excluder, warnedexcl, categ, S);
    end %- firstinit

%   ---------------------------------------------
    function [bestset,X,logdetX] = firstgen1(~,S)
        if isempty(S)
            S = RandStream.getGlobalStream;
        end
        [bestset,X,logdetX] = gen1design(setstart,Xstart,S);
    end %-firstgen1

%   ----------------------------------------
    function reductionVar = loopBody(iter,S)
        % Print "try" number and worker to screen.
        if ~quiet
            % Bump iter for "try" display because an initial iteration
            % was performed prior to first invocation of loopBody.
            if useParallel
                if usePool
                    labindx = internal.stats.parallel.workerGetValue('workerID');
                else
                    labindx = 1;
                end
                disp(getString(message('stats:cordexch:WorkerTry', labindx, iter+1)));
            else
                settry(iter+1);
            end
        end
        
        if isempty(S)
            S = RandStream.getGlobalStream;
        end
        
        [xinitTemp,warnedexcl] = randstart(bnds, nlevels, nruns, nfactors, ...
            excluder, warnedexcl, categ, S);
        
        setstartTemp = setstart;  % must make a temp. var in parfor loop
        setstartTemp(nobs+1:end,1:nfactors) = xinitTemp;
        
        % Generate the term values ('model' will already be in matrix form)
        [XstartTemp,~,~] = x2fx(setstartTemp,model,categ,catlevels);
        
        % Generate a new design from here by coordinate exchange
        % (this nested function uses factorterms and some other variables)
        
        [A,B,C] = gen1design(setstartTemp,XstartTemp,S);
        
        reductionVar{1} = warnedexcl;
        reductionVar{2} = {C iter A B};
    end %-loopBody

end % of outer function
     
% --------------------------------------------
function R = adjustr(R)
%ADJUSTR Adjust R a little bit so it will be non-singular

diagr = abs(diag(R));
p = size(R,2);
smallval = sqrt(eps)*max(diagr);
t = (diagr < smallval);
if any(t)
   tind = (1:p+1:p^2);
   R(tind(t)) = smallval;
end
end

% --------------------------------------------------
function [xinit,warnedexcl] = randstart(bnds, nlevels, nruns, nfactors, ...
                                        excluder, warnedexcl, categ, s)
%RANDSTART Create random but feasible initial design
xinit = zeros(0,nfactors);
startiter = 1;
if isempty(excluder)
   blocksize = nruns;
else
   blocksize = max(1000,nruns);
end
catcols = ismember(1:nfactors, categ);
contcols = ~catcols;
ncatcols = sum(catcols);
ncontfactors = sum(contcols);
if iscell(bnds)  % only need extremes here, convert to matrix for convenience
   tempbnds = zeros(2,nfactors);
   for j=1:nfactors
       v = bnds{j};
       tempbnds(1,j) = v(1);
       tempbnds(2,j) = v(end);
   end
   bnds = tempbnds;
end
bnds(1,catcols) = 1;
bnds(2,catcols) = nlevels(catcols);

% Repeatedly attempt to generate points, discarding excluded ones
while(size(xinit,1)<nruns && startiter<=20)
    block = zeros(blocksize,nfactors);
    block(:,contcols) = repmat(bnds(1,contcols),blocksize,1) + ...
            repmat(diff(bnds(:,contcols),1,1),blocksize,1) ...
                .* rand(s,blocksize,ncontfactors);
    block(:,catcols) = ceil(rand(s,blocksize,ncatcols) ...
                            .* repmat(bnds(2,catcols),blocksize,1));
    if ~isempty(excluder)
       excluderows = excluder(block);
       if startiter==20 % last time only, save discards in case needed
           discards = block(excluderows,:);
       end
       block(excluderows,:) = [];
    end       
    xinit = [xinit; block];
    startiter = startiter + 1;
end
if size(xinit,1)<nruns
    % Use bad points if absolutely necessary, hoping they can be removed
    % during the coordinate exchange process.
    if ~warnedexcl
       % Record the occurrence of the bad initial design so that a warning 
       % can be issued before exiting CORDEXCH.
       warnedexcl = true;
    end
    N = nruns - size(xinit,1);
    xinit = [xinit; discards(1:N,:)];
else
    xinit = xinit(1:nruns,:);
end
end

% -------------------------------------------------------
function xinit = levels2numbers(xinit,bnds,categ,nlevels)
% Renumber levels from values in bnds to numbers in 1:nlevels

for j=1:length(categ)
    cj = categ(j);
    if iscell(bnds)
       setlist = bnds{cj};
    else
       setlist = linspace(bnds(1,cj),bnds(2,cj),nlevels(cj));
    end
    [foundit,catrow] = ismember(xinit(:,cj),setlist);
    if any(~foundit)
       error(message('stats:cordexch:BadInit', cj));
    end
    xinit(:,cj) = catrow;
end
end

% -------------------------------------------------------
function settings = numbers2levels(settings,bnds,categ,nlevels)
% Renumber levels from numbers in 1:nlevels to values in bnds

for j=1:length(categ)
    cj = categ(j);
    if iscell(bnds)
        setlist = bnds{cj};
    else
        setlist = linspace(bnds(1,cj),bnds(2,cj),nlevels(cj));
    end
    settings(:,cj) = setlist(settings(:,cj));
end
end

% -------------------------------------------------------
function reductionVal = reductionOperator(reductionVal, update)
    reductionVal{1} = reductionVal{1} | update{1};
    reductionVal{2} = internal.stats.parallel.pickLarger(reductionVal{2},update{2});
end %-reductionOperator

