function schema
%


% Copyright 2003-2013 The MathWorks, Inc.
mlock
pk = findpackage('stats');

% Create a new class

c = schema.class(pk, 'outlierdb');

schema.prop(c, 'current', 'ustring');
p=schema.prop(c, 'listeners', 'MATLAB array');
p.AccessFlags.Serialize = 'off';
