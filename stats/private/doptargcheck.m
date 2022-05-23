function [varargout] = doptargcheck(okargs,nfactors,nruns,varargin)
%DOPTARGCHECK Check arguments common to D-optimal design functions
%   Utility function used by CORDEXCH and ROWEXCH.

%   The undocumented parameters 'start' and 'covariates' are used when
%   this function is called by the daugment and dcovary functions to
%   implement D-optimal designs with fixed rows or columns, and are
%   supported only by the caller CORDEXCH.
   
%   Copyright 2005-2014 The MathWorks, Inc.


%
% Define all arguments and their defaults
%
z = zeros(0,nfactors);
pnames = {'start'      'covariates'  'display'   'init'       'maxiter' ...
          'tries'      'bounds'      'levels'    'excludefun' 'categorical' 'useparallel'};
pdflts = {z            []            []          []           10            ...
          1            []            []          []           []            true};

% Take a subset if caller allows only some of them
if ~isempty(okargs)
   allowed = ismember(pnames,okargs);
   pnamesub = pnames(allowed);
   pdfltsub = pdflts(allowed);
else
   pnamesub = pnames;
   pdfltsub = pdflts;
   allowed = true(size(pnames));
end

% Get specified values of allowed parameters
varargout = cell(1,length(pnamesub));
[varargout{:}] = internal.stats.parseArgs(pnamesub, pdfltsub, varargin{:});

% Deal these out to separate variables for convenience in error checking
allvalues = pdflts;
allvalues(allowed) = varargout;
[startdes,   covariates, dodisp,  settings, maxiter,  ...
 tries,      bnds,       nlevels, excluder, categ, useParallel] = deal(allvalues{:});

%
% Perform parameter checks
%

% If there is no 'display' parameter in command line, and if the caller
% supports the 'display' parameter, set the default depending on whether
% parallel or serial computation.
if isempty(dodisp)
    dispix = find(strcmpi('display',pnamesub));
    if ~isempty(dispix)
        if useParallel
            dodisp = 'off';
        else
            dodisp = 'on';
        end
        varargout{dispix} = dodisp;
    end
elseif ~isequal(dodisp,'on') && ~isequal(dodisp,'off')
    m = message('stats:doptargchk:BadDisplay');
    throwAsCaller(MException(m.Identifier,'%s',getString(m)));
end


% Check number of iterations
if ~isnumeric(maxiter) || ~isscalar(maxiter) || maxiter<1 || ~isreal(maxiter)
    m = message('stats:doptargchk:BadMaxIter');
    throwAsCaller(MException(m.Identifier,'%s',getString(m)));
end

% Check number of tries
if ~isnumeric(tries) || ~isscalar(tries) || tries<1 || ...
   ~isreal(tries) || tries~=round(tries) || ~isfinite(tries)
    m = message('stats:doptargchk:BadTries');
    throwAsCaller(MException(m.Identifier,'%s',getString(m)));
end

% Check exclusion function, if any
if ~isempty(excluder)
   if ~isa(excluder,'function_handle') && ...
      ~(ischar(excluder) && exist(excluder,'file')>0)
      m = message('stats:doptargchk:BadExcludeFun');
      throwAsCaller(MException(m.Identifier,'%s',getString(m)));
   end
end

% Check dimensions of initial design, if any
nvars = nfactors + size(covariates,2);
if ~isempty(settings)
   if size(settings,2)~=nvars || ~isnumeric(settings) || size(settings,1)~=nruns
      m = message('stats:doptargchk:BadInit');
      throwAsCaller(MException(m.Identifier,'%s',getString(m)));
   end
end

% Check dimensions of the specified rows and columns, if any
if ~isempty(startdes) && (size(startdes,2)~=nfactors || ~isnumeric(startdes))
    m = message('stats:doptargchk:BadStart');
    throwAsCaller(MException(m.Identifier,'%s',getString(m)));
end
if ~isempty(covariates) && (size(covariates,1)~=nruns || ~isnumeric(covariates))
    m = message('stats:doptargchk:InputSizeMismatch');
    throwAsCaller(MException(m.Identifier,'%s',getString(m)));
end

% Check the specification of categorical factors
if ~isempty(categ)
   if ~isvector(categ) || ~all(ismember(categ,1:nfactors))
      m = message('stats:doptargchk:BadCategorical');
      throwAsCaller(MException(m.Identifier,'%s',getString(m)));
   end
end

% Check the bounds for each factor, if specified
boundsrow = find(strcmpi('bounds',pnamesub));
if isempty(bnds) && ~isempty(boundsrow)
   bnds = ones(2,nfactors);
   bnds(1,:) = -1;
   if ~isempty(categ)
      bnds(1,categ) = 1;
      if isempty(nlevels)
          lev = 2;
      elseif isscalar(nlevels)
          lev = nlevels;
      else
          lev = nlevels(categ);
      end
      bnds(2,categ) = lev;
   end
   varargout{boundsrow} = bnds;
end
if isempty(boundsrow)
   % ok
elseif iscell(bnds)
   if ~isvector(bnds) || length(bnds)~=nfactors
      m = message('stats:doptargchk:BadBoundsLength');
      throwAsCaller(MException(m.Identifier,'%s',getString(m)));
   end
   obsnlevels = cellfun('prodofsize',bnds);
   if any(obsnlevels<2)
      m = message('stats:doptargchk:BadBoundsValues');
      throwAsCaller(MException(m.Identifier,'%s',getString(m)));
   end
elseif isnumeric(bnds)
   if ~isequal(size(bnds),[2 nfactors])
      m = message('stats:doptargchk:BadBoundsMatrix');
      throwAsCaller(MException(m.Identifier,'%s',getString(m)));
   end
else
   m = message('stats:doptargchk:BadBoundsType');
   throwAsCaller(MException(m.Identifier,'%s',getString(m)));
end    

% Check the numbers of levels for each factor, if specified
levrow = find(strcmpi('levels',pnamesub));
if iscell(bnds)
   if (numel(nlevels)~=numel(obsnlevels) && ...
              ~isscalar(nlevels) && ~isscalar(obsnlevels))|| ...
           ~all(nlevels(:)==obsnlevels(:))
      if ~isempty(nlevels) 
         warning(message('stats:doptargchk:BadLevels'));
      end
   nlevels = obsnlevels;
   end
end
if ~isempty(nlevels) && ~isempty(levrow)
   if ~isscalar(nlevels) && ~isvector(nlevels)
      m = message('stats:doptargchk:BadLevelsSize');
      throwAsCaller(MException(m.Identifier,'%s',getString(m)));
   elseif ~isscalar(nlevels) && length(nlevels)~=nfactors
      m = message('stats:doptargchk:BadLevelsLength');
      throwAsCaller(MException(m.Identifier,'%s',getString(m)));
   elseif any(nlevels<2 | nlevels~=round(nlevels))
      m = message('stats:doptargchk:BadLevelsValues');
      throwAsCaller(MException(m.Identifier,'%s',getString(m)));
   end
   if isscalar(nlevels)
       nlevels = nlevels*ones(1,nfactors);
   end
end
if ~isempty(levrow)
   varargout{levrow} = nlevels;
end

