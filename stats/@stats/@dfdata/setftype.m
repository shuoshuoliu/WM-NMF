function setftype(ds,ftype)
%SETFTYPE Set function type for plotting of data set


%   Copyright 2003-2011 The MathWorks, Inc.

% Must have two inputs of proper type
if nargin<2 || ~ischar(ftype)
   error(message('stats:dfdata:setftype:FunctionTypeMissingOrNotChar'));
end

% Verify input value
oktypes = {'cdf' 'icdf' 'pdf' 'probplot' 'survivor' 'cumhazard'};
ftype = internal.stats.getParamVal(ftype,oktypes,'function type');

% If changing type, clear out information no longer correct
oldtype = ds.ftype;
if ~isequal(oldtype,ftype)
   clearplot(ds);
end

% Now we can set the new value
ds.ftype = ftype;
