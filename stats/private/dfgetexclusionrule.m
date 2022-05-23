function hExcl = dfgetexclusionrule(ename)
%GETEXCLUSIONRULE Get an exclusion rule by name


% Copyright 2003-2004 The MathWorks, Inc.

db = getoutlierdb;
hExcl = find(db,'name',ename);

