function clearplot(hFit)
%CLEARPLOT Clear plot information from fit object


%   Copyright 2003-2008 The MathWorks, Inc.

hFit.x = [];
hFit.y = [];
if ~isempty(hFit.linehandle) && ishghandle(hFit.linehandle)
   delete(hFit.linehandle);
end
