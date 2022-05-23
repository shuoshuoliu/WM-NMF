function minmax = xlim(ds)
%XLIM Return the X plot limits for this dataset


% Copyright 2003-2004 The MathWorks, Inc.

if isempty(ds.xlim)
   minmax = ds.datalim;   % use data limits if nothing else available
else
   minmax = ds.xlim;
end
