function [s,emsg,classname] = getdistributions(distname,douser,dostore,doforce)
%GETDISTRIBUTIONS Get structure defining the distributions supported by
%distributionFitter

%   Copyright 2003-2020 The MathWorks, Inc.


% If a struct was passed in, store this for later use
emsg = '';
if nargin>0 && isstruct(distname)
   dfgetset('alldistributions',distname);
   return
end

classname = '';
if nargin==1 && length(distname)>1
    switch(distname)
        case 'ev', distname = 'extreme value';
        case 'gev', distname = 'generalized extreme value';
        case 'gp', distname = 'generalized pareto';
        case 'hn', distname = 'half normal';
        case 'nbin', distname = 'negative binomial';
        case 'wbl', distname = 'weibull';
    end
    try
        dist = prob.ProbabilityDistributionRegistry.get(distname);
        if dist.toolboxcompliant
            classname = dist.classname;
            s = makespec(classname);
            return
        end
    catch ME
        if  ~strcmp(ME.identifier,'stats:internal:getParamVal:BadValueListChoices') ...
         && ~strcmp(ME.identifier,'stats:internal:getParamVal:BadValue')
            rethrow(ME)
        end
    end
end

% Get old value if already created and stored
if (nargin<3 || dostore) && (nargin<4 || ~doforce)
   s = dfgetset('alldistributions');
else
   s = '';
end

% If not created yet, create it now
if isempty(s)

   if nargin<2 || douser
      % Get user-defined distributions (won't be done if we already
      % had a distribution list created before this function was called)
      try
          s = stats.internal.dfit.getuserdists(s);
      catch ME
          if nargout>=2
              emsg = ME.message;
              causes = ME.cause;
              for j=1:length(causes)
                  emsg = sprintf('%s\n%s',emsg,causes{j}.message);
              end
          else
              rethrow(ME)
          end
      end
   end
   
   % Append new distributions, replacing old ones as needed
   s = appendnew(s);

   % Sort by name
   lowernames = lower(char(s.name));
   [~, ind] = sortrows(lowernames);
   s = s(ind);

   % Store it for next time
   if nargin<3 || dostore
       dfgetset('alldistributions',s);
   end
end

if nargin>0 && ~isempty(distname)
   % Return only the distribution(s) requested, not all of them
   allnames = {s.code};
   distnum = strncmpi(distname, allnames, length(distname));
   s = s(distnum);
end

% ------------------------------------
function a = appendnew(b)
% Add two new fields to the pre-R2013a spec structures
if ~isempty(b)
    b(end).logci = [];
    b(end).fittable = [];
    b(end).plim = [];
    oldcodes = {b.code};
else
    oldcodes = {};
end
    
newnames = prob.ProbabilityDistributionRegistry.list('fittable');
for j=1:length(newnames)
    dist = prob.ProbabilityDistributionRegistry.get(newnames{j});
    if ~dist.toolboxcompliant
        continue
    end
    classname = dist.classname;
    spec = makespec(classname);
    newcode = spec.code;
    loc = find(strcmp(newcode,oldcodes));
    if ~isempty(loc)
        b(loc) = [];
        oldcodes(loc) = [];
    end
        oldcodes = [oldcodes(:); newcode];
   
    b = [b(:); spec];
end
a = b;

% ------------ generalized Pareto functions are a special case
function [phat,pci] = localgpfit(x,theta,alpha,varargin)
%LOCALGPFIT Version of gpfit that handles a fixed threshold param

if any(x<=theta)
    error(message('stats:gpfit:BadDataThreshold'));
end

if nargout < 2
    phat = [gpfit(x-theta,alpha,varargin{:}) theta];
else
    [phat,pci] = gpfit(x-theta,alpha);
    phat = [phat theta];
    pci = [pci [theta; theta]];
end

function [nlogL,acov] = localgplike(params,data)
%LOCALGPLIKE Version of gplike that handles a fixed threshold param

theta = params(3);
params = params(1:2);
if nargout < 2
    nlogL = gplike(params,data-theta);
else
    [nlogL,acov] = gplike(params,data-theta);
    acov = [acov [0; 0]; [0 0 0]];
end

function [range,closed] = localgpsupport(params)
k = params(1);
theta = params(3);
if k<0
    sigma = params(2);
    range = sort([theta, theta-sigma/k]);
    closed = [true true];
else
    range = [theta Inf];
    closed = [false false];
end

function s = makespec(classname)
% Make old-style spec from new prob dist object
s = eval([classname '.getInfo()']);


% functions below need to be present when an old ProbDist object is loaded
function localbinofit(~)
function localbinosupport(~)
function localgevsupport(~)
function localbinolike(~)

