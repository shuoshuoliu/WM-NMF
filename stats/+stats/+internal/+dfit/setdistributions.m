function setdistributions(dft,dists)
%SETDISTRIBUTIONS Set distribution information into the gui

%   Copyright 2003-2020 The MathWorks, Inc.

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
   dft.addFitType('addparamfit',a.name, a.code, a.pnames, a.pdescription, ...
                  a.prequired, a.support(1), a.support(2), ...
                  a.closedbound(1), a.closedbound(2), ...
                  a.censoring, ~a.iscontinuous);
end

dft.addFitType('addsmoothfit', nonparname, 'nonparametric',...
               {'a'}, {'a'}, {'true'}, -Inf, Inf, false, false, true, false);

for j=insertpos:length(dists)
   a = dists(j);
   dft.addFitType('addparamfit',a.name, a.code, a.pnames, a.pdescription, ...
                  a.prequired, a.support(1), a.support(2), ...
                  a.closedbound(1), a.closedbound(2), ...
                  a.censoring, ~a.iscontinuous);
end
