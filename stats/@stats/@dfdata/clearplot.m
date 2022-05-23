function clearplot(ds)
%CLEARPLOT Clear current plot data


%   Copyright 2003-2008 The MathWorks, Inc.

ds.plotx = [];
ds.ploty = [];
if ~isempty(ds.line) && ishghandle(ds.line)
   delete(ds.line);
end
