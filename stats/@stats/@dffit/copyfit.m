function new=copyfit(original)
%COPYFIT Copy constructor


% Copyright 2003-2020 The MathWorks, Inc.

s = settings;
if (s.stats.DistributionFitter.LegacyDistributionFitter.ActiveValue == false) 
    new = copyfitApp(original);
    return;
end
    

% it may be coming in with a java bean wrapper
original=handle(original);

% Determine a new, unique name.
name = original.name;
taken = true;
i=1;

% keep from prepending multiple "copy x of"'s.
COPY = sprintf(' %s ',getString(message('stats:dfstrings:assignment_Copy')));
ind=strfind(name,COPY);
if ~isempty(ind)
    name=name(1:(ind(end)-1));
end

% search for first unique name
fitdb = dfswitchyard('getfitdb');
while taken
   newName = sprintf('%s%s%d',name,COPY,i);
   a = down(fitdb);
   taken = false;
   while(~isempty(a))
      if isequal(a.name,newName)
         taken = true;
         break
      end
      a = right(a);
   end
   i=i+1;
end

new = initdistfit(stats.dffit,newName);

% copy all fields from the old to the new, but
% skip any of the ones on the toskip list.
fields = fieldnames(new);
toskip = {'listeners' 'name' 'plot' 'linehandle' 'ColorMarkerLine' 'boundline'};
for i=1:length(fields)
   if ~ismember(fields{i},toskip)
      set(new,fields{i},get(original,fields{i}));
   end
end

fitdb = dfswitchyard('getfitdb');
connect(new, fitdb,'up');

% restore other properties
if original.plot
   new.plot=1;
end

