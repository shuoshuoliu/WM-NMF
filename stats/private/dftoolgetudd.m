function javaudd=dftoolgetudd(uddcmd,varargin);


%   Copyright 2003-2004 The MathWorks, Inc.

% unwrap any UDD objects
for i=1:length(varargin)
   if isa(varargin{i}, 'com.mathworks.jmi.bean.UDDObject')
      varargin{i}=handle(varargin{i});
   end
end

% wrap the return UDD object
if nargin == 1
   javaudd=java(eval(uddcmd));
else
   javaudd=java(feval(uddcmd,varargin{:}));
end
