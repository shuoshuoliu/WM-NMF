function [prototype] = numericPrototype(X, varNums)
%NUMERICPROTOTYPE Return a prototype numeric array based on the input data
%   PROTOTYPE = NUMERICPROTOTYPE(X) returns an empty single or double array
%   based on the numeric data in X. X may be an array or a table.
%
%   If X is a table, PROTOTYPE = NUMERICPROTOTYPE(X,VARNUMS) returns a
%   prototype based on the numeric data in TBL(:,VARNUMS).
%
%   The storage type (for example, gpuArray) of any numeric data in X will
%   be preserved.
%
%   Examples:
%   internal.stats.numericPrototype(logical(1)) => empty double
%   internal.stats.numericPrototype(single(1)) => empty single
%   internal.stats.numericPrototype(double(1)) => empty double
%   internal.stats.numericPrototype(gpuArray(double(1))) => empty double
%   internal.stats.numericPrototype(table(gpuArray(1), logical(1))) => empty gpuArray double
%   internal.stats.numericPrototype(table(gpuArray(1), single(1))) => empty gpuArray single
%
% See also internal.stats.dominantType.

%   Copyright 2020 The MathWorks, Inc.

if istable(X)
    if nargin > 1
        dataPrototypes = X([],varNums);
    else
        dataPrototypes = X([],:);
    end
    % dataPrototypes(:, vartype('numeric')) does not return numeric
    % gpuArray variables, we need to use the following workaround:
    numericDataPrototypes = dataPrototypes(:,varfun(@isnumeric,dataPrototypes,'OutputFormat','uniform'));
    
    % Get an empty array of the dominant type of numericDataPrototypes.
    % Single will dominate double.
    prototype = numericDataPrototypes{[],:};
    
    % Ensure we return a 0-by-0 array
    prototype = prototype([]);
else
    prototype = X([]);
end

% Promote the type of anything other than single (for example, logical,
% uint64) to double.  The storage type (for example, gpuArray) will be
% preserved.
if ~isfloat(prototype)
    prototype = double(prototype);
end
end

