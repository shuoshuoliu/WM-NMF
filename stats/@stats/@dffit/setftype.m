function setftype(hFit,ftype)
%SETFTYPE Method to set function type for plotting of fit


%   Copyright 2003-2011 The MathWorks, Inc.

% Must have two inputs of proper type
if nargin<2 || ~ischar(ftype)
   error(message('stats:dffit:setftype:FunctionTypeMissingOrNotChar'));
end

% Verify input value
oktypes = {'cdf' 'icdf' 'pdf' 'probplot' 'survivor' 'cumhazard'};
ftype = internal.stats.getParamVal(ftype,oktypes,'function type');

% If changing type, clear out information no longer correct
oldtype = hFit.ftype;
if ~isequal(oldtype,ftype)
   clearplot(hFit);
end

% Now we can set the new value
hFit.ftype = ftype;
