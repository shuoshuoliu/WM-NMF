function dfsetdistributions(dft,dists)
%DFSETDISTRIBUTIONS Set distribution information into the gui


%   Copyright 2003-2012 The MathWorks, Inc.

% Store for later use in M
dfgetset('alldistributions',dists);

% Set into the gui, placing nonparametric fit into sorted list
dft.clearFitTypes;

j = strcmp('kernel',{dists.code});
dists(j) = [];

nonparname = getString(message('stats:dfittool:NameNonparametric'));
allnames = [{dists.name}, {nonparname}];
[~,sortidx] = sort(allnames);
insertpos = find(sortidx == length(allnames));

for j=1:insertpos-1
   a = dists(j);
   preq = getpreq(a.prequired);   
   dft.addFitType('addparamfit',a.name, a.code, a.pnames, a.pdescription, ...
                  preq, a.support(1), a.support(2), ...
                  a.closedbound(1), a.closedbound(2), ...
                  a.censoring, ~a.iscontinuous);
end

dft.addFitType('addsmoothfit', nonparname, 'nonparametric',...
               {'a'}, {'a'}, {'true'}, -Inf, Inf, false, false, true, false);

for j=insertpos:length(dists)
   a = dists(j);
   preq = getpreq(a.prequired);   
   dft.addFitType('addparamfit',a.name, a.code, a.pnames, a.pdescription, ...
                  preq, a.support(1), a.support(2), ...
                  a.closedbound(1), a.closedbound(2), ...
                  a.censoring, ~a.iscontinuous);
end


% ---- Having trouble with booleans, so call java with text
function prequired=getpreq(boolvec)
prequired = cell(size(boolvec));
for j=1:length(boolvec)
   if boolvec(j)
      txt = 'true';
   else
      txt = 'false';
   end
   prequired{j} = txt;
end
