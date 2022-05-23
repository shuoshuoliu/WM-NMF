function plotargchk(varargin)
% Internal utility to check plotting arguments

%   Copyright 2012-2020 The MathWorks, Inc.


if nargin==0
    return
end

if mod(nargin,2)==1
    m = message('stats:internal:parseArgs:WrongNumberArgs');
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
end

nvPairNames = varargin(1:2:end);
[nvPairNames{:}] = convertStringsToChars(nvPairNames{:});
if ~iscellstr(nvPairNames)
    m = message('stats:internal:parseArgs:IllegalParamName');
    throwAsCaller(MException(m.Identifier, '%s', getString(m)));
end
