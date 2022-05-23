function [news,newrows]=dfgetuserdists(olds,userfun)
%GETUSERDISTS Get user-defined distributions for dfittool
%   [NEWS,NEWROWS]=GETUSERDISTS(OLDS) appends user-defined distribution
%   information to the existing distribution information in the structure
%   OLDS and returns the combined information in the structure NEWS.
%   NEWROWS is a vector of the indices that are new or have changed.
%
%   [...]=GETUERDISTS(OLDS,USERFUN) uses the function USERFUN in
%   place of the default @dfittooldists.


%   Copyright 2003-2011 The MathWorks, Inc.

% If no user-defined function given, use the default if it exists
if nargin<2
   if exist('dfittooldists','file') 
      userfun = @dfittooldists;
   else
      userfun = [];
   end
end

news = olds;
newrows = [];

if isempty(userfun)
   return
end

% First try running the user's function
if isa(userfun,'function_handle')
   userfunname = func2str(userfun);
else
   userfunname = userfun;
end
try
   s = feval(userfun);
catch ME
   m = message('stats:dfittool:BadUserDistributions',userfunname);
   ME2 = MException(m.Identifier,'%s',getString(m));
   throw(addCause(ME2,ME));
end

% Next make sure the result is a structure
if ~isstruct(s)
   error(message('stats:dfittool:StructureRequired', userfunname));
end

newfields = fieldnames(s);
numnewdists = length(s);
numolddists = length(olds);
requiredfields = {'name' 'pnames' 'cdffunc' 'pdffunc' 'invfunc' 'fitfunc'};

% Next make sure the result has all required fields
for j=1:length(requiredfields)
   if ~any(strcmp(requiredfields{j},newfields))
      error(message('stats:dfittool:MissingField', requiredfields{ j }, userfunname));
   end
end

% Make sure the field values are as expected
checked = cell(1,numnewdists);
for j=1:numnewdists
   sj = s(j);
   sj = checkfields(sj);
   checked{j} = sj;
end

% See if we are overwriting existing fields
if numolddists>0
   oldnames = {olds.name};
   oldcodes = {olds.code};
else
   oldnames = cell(0);
   oldcodes = cell(0);
end

for j=1:numnewdists
   % See if the proposed name or code is in use already
   sj = checked{j};
   newfields = fieldnames(sj);  % may need updating since previous assignment
   name = sj.name;
   code = sj.code;
   oldnrow = find(strcmp(name,oldnames));
   oldcrow = find(strcmp(code,oldcodes));

   if isempty(oldcrow)
      if ~isempty(oldnrow)
         newrow = oldnrow;          % replace distribution with same name
      else
         newrow = numolddists+1;    % new distribution
      end
   else
      % Trying to re-define an existing distribution
      if ~isempty(oldnrow) && ~isequal(oldcrow,oldnrow)
         error(message('stats:dfittool:DuplicateName', code));
      end
      newrow = oldcrow;
   end

   % Update fields in old structure.
   % Can't concatenate with [] if field names differ.
   for fieldnum = 1:length(newfields)
      fieldname = newfields{fieldnum};
      olds(newrow).(fieldname) = sj.(fieldname);
   end

   % Update arrays to guard against duplicates within the new structure
   if newrow>numolddists
      oldnames = [oldnames {name}];
      oldcodes = [oldcodes {code}];
      numolddists = numolddists+1;
   end
   newrows = [newrows newrow];
end

% Return updated structure as new structure
news = olds;
      
% ------------------------------------------
function sj = checkfields(sj)
%CHECKFIELDS Check that a distribution structure's fields are all valid

% Check required fields
testnames = {'name'};
for j=1:length(testnames)
   field = testnames{j};
   val = checkstring(sj,field,'',false);
   sj.(field) = val;
end

field = 'pnames';
val = checktext(sj,field,'',false,[]);
sj.(field) = val;
nparams = length(val);

testnames = {'pdffunc' 'cdffunc' 'invfunc' 'fitfunc'};
for j=1:length(testnames)
   field = testnames{j};
   val = checkfunc(sj,field,false);
   sj.(field) = val;
end

% Check optional fields and fill in defaults
sj = dfprobspecdefaults(sj);

testnames = {'likefunc' 'logcdffunc' 'loginvfunc' ...
             'randfunc' 'checkparam' 'cifunc' 'supportfunc'};
for j=1:length(testnames)
   field = testnames{j};
   val = checkfunc(sj,field,true);
   sj.(field) = val;
end

field = 'prequired';
val = checklogical(sj,field,false(1,nparams),nparams);
sj.(field) = val;

field = 'pdescription';
val = checktext(sj,field,{},true,nparams);
sj.(field) = val;

field = 'closedbound';
val = checklogical(sj,field,[false false],2);
sj.(field) = val;

field = 'support';
val = checksupport(sj);
sj.(field) = val;


% ------------------------------------------
function val = checklogical(s,field,default,nvals)
%CHECKLOGICAL Check that a field has a valid logical value

if nargin<4
   nvals = 1;
end

if ~isfield(s,field) || isempty(s.(field))
   val = default;
else
   val = s.(field);
end

if (numel(val) ~= nvals)
   error(message('stats:dfittool:WrongSize',field,nvals));
elseif ~(islogical(val) || isnumeric(val))
   error(message('stats:dfittool:NotLogical',field));
else
   val = (val ~= 0);
end


% ------------------------------------------
function val = checktext(s,field,default,optional,nvals)
%CHECKTEXT Check that a field has a value that is an array of strings

if nargin<5
   nvals = 1;
end

if ~isfield(s,field) || isempty(s.(field))
   if optional
      val = default;
      return
   else
      error(message('stats:dfittool:EmptyNotAllowed', field));
   end
end
val = s.(field);
if iscellstr(val)
   if ~isempty(nvals) && numel(val)~=nvals
      error(message('stats:dfittool:BadSize', field, nvals));
   end
elseif ischar(val)
   if ~isempty(nvals) && size(val,1)~=nvals
      error(message('stats:dfittool:BadSize', field, nvals));
   else
      val = cellstr(val);
   end
else
   error(message('stats:dfittool:NotCharacterCell', field));
end

% ------------------------------------------
function val = checkstring(s,field,default,optional)
%CHECKSTRING Check that a field has a valid string

if ~isfield(s,field) || isempty(s.(field))
   if optional
      val = default;
      return
   else
      error(message('stats:dfittool:EmptyNotAllowed', field));
   end
end
val = s.(field);
if ~ischar(val) || (~isequal(size(val), [1,length(val)]))
   error(message('stats:dfittool:NotString', field));
end


% ------------------------------------------
function val = checkfunc(s,field,optional)
%CHECKFUNC Check that a field has a valid function value

val = '';
if ~isfield(s,field) || isempty(s.(field))
   if ~optional
      error(message('stats:dfittool:EmptyNotAllowed', field));
   end
   return
end

val = s.(field);
if ~isa(val,'function_handle') || ~isscalar(val)
   error(message('stats:dfittool:NotFunctionHandle', field));
end

% ------------------------------------------
function val = checksupport(s)
%CHECKSUPPORT Check that the support field is valid

default = [-Inf Inf];
field = 'support';
if ~isfield(s,field) || isempty(s.(field))
   val = default;
   return;
end
val = s.(field);
if ~isnumeric(val) || numel(val)~=2 || (val(1)>=val(2)) || any(isnan(val))
   error(message('stats:dfittool:BadSupport'));
end
