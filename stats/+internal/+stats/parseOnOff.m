function ison = parseOnOff(onoff,paramname)
% Process argument that can be true, false, 0, 1, 'on', or 'off', and
% return true if it is 'on', 1 or true.

%   Copyright 2011-2015 The MathWorks, Inc.

ison = statslib.internal.parseOnOff(onoff,paramname);
