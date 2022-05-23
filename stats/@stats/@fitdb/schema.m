function schema
%

% Copyright 2003-2013 The MathWorks, Inc.
mlock
pk = findpackage('stats');

% Create a new class

c = schema.class(pk, 'fitdb');

schema.prop(c, 'current', 'ustring');

p=schema.prop(c, 'listeners', 'MATLAB array');
p.AccessFlags.Serialize = 'off';

p=schema.prop(c, 'newOptions', 'handle');
p.AccessFlags.Serialize = 'off';

p=schema.prop(c, 'newModel', 'string');
p.AccessFlags.Serialize = 'off';

p=schema.prop(c, 'newCoeff', 'string vector');
p.AccessFlags.Serialize = 'off';

p=schema.prop(c, 'newProps', 'string vector');
p.AccessFlags.Serialize = 'off';
