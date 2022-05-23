function [err, commonSize, numElements] = statsizechkM(nparams,varargin)
%STATSIZECHK Check for compatible array sizes.
%   [ERR,COMMONSIZE,NUMELEMENTS] = STATSIZECHK(NPARAMS,A,B,...,M,N,...) or
%   [ERR,COMMONSIZE,NUMELEMENTS] = STATSIZECHK(NPARAMS,A,B,...,[M,N,...])
%   in effect computes size( A + B + ... + zeros(M,N,...) ), and catches
%   any size mismatches.  NPARAMS is the number of array input arguments.
%
%   This version is used when one or more inputs are not standard numeric
%   types.

%   Copyright 2015-2018 The MathWorks, Inc.

err = 0;
commonSize = [];
numElements = 0;

if ~isscalar(nparams) || ~isnumeric(nparams)
    error(message('stats:statsizechk:BadNParams'));
end

% Get sizes from parameters
paramSizes = cellfun(@size, varargin(1:nparams), 'UniformOutput', false);

% Get output size from trailing args (if any). If specified, other inputs
% must match this size (or be scalar). If not specified, work out the final
% size from the first non-scalar input.
outputSize = [];
if nargin > nparams+1
    outputSize = iExtractTrailingDims(varargin(nparams+1:end));
end

% Scan each input looking for inconsistencies
for ii=1:numel(paramSizes)
    % Ignore scalars
    if isequal(paramSizes{ii}, [1 1])
        continue;
    end
    % For non-scalars, compare to commonSize if set, or set it.
    if isempty(outputSize)
        outputSize = paramSizes{ii};
    else
        if ~isequal(paramSizes{ii}, outputSize)
            % Inconsistency.
            err = 1;
            return
        end
    end
end

% If outputSize is still not set then everything is scalar
if isempty(outputSize)
    commonSize = [1 1];
    numElements = 1;
else
    commonSize = outputSize;
    numElements = prod(commonSize);
end




function sz = iExtractTrailingDims(dimArgs)
% Helper to form a dimension vector from some trailing dimension arguments
%
% form can be M              -> MxM
%             M, N, P, ...   -> MxNxPx...
%             [M, N, P, ...] -> MxNxPx...
if numel(dimArgs)==1
    if isscalar(dimArgs{1})
        % Treat a scalar N as meaning NxN
        sz = dimArgs{1}*[1 1];
    else
        % One vector dimension
        sz = dimArgs{1};
    end
else
    % Multiple dims. Each should be a scalar.
    if any(~cellfun(@isscalar, dimArgs))
        error(message('MATLAB:NonScalarInput'));
    end
    % Concatenate the scalars to form the dimension vector
    sz = [dimArgs{:}];
end
% Strip any trailing ones
if numel(sz)>2
    lastNonUnity = find(sz~=1, 1, 'last');
    sz = sz(1:lastNonUnity);
end 