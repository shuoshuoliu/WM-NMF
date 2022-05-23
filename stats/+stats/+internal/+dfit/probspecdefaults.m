function s = probspecdefaults(s)
%PROBSPECDEFAULTS Fill default entries into probability spec


%   Copyright 2010-2011 The MathWorks, Inc.

% Make sure logical fields have defaults rather than empties
testnames = {'hasconfbounds' 'iscontinuous' 'islocscale' 'uselogpp' ...
             'censoring'     'paramvec'     'optimopts'};
defaults  = {false           true           false        false      ...
             false           true           false};
for j=1:length(testnames)
   field = testnames{j};
   default = defaults{j};
   if ~isfield(s,field) || isempty(s.(field))
      val = default;
   else
      val = s.(field);
      if ~isscalar(val) || ~(islogical(val) || isnumeric(val))
          error(message('stats:dfittool:NotLogical', field));
      end
      val = (val~=0);
   end
   s.(field) = val;
end

% Check other optional fields and fill in defaults
testnames = {'code'};
defaults  = {lower(s.name)};
for j=1:length(testnames)
   field = testnames{j};
   default = defaults{j};
   if ~isfield(s,field) || isempty(s.(field))
      val = default;
   else
      val = s.(field);
      if ~isrow(val) || ~ischar(val)
          error(message('stats:dfittool:NotCharacter', field));
      end
   end
   s.(field) = val;
end
