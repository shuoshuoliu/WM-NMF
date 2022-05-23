function [varargout] = checkCompatibleSizes(varargin)
% checkCompatibleSizes Check that all inputs have a size that is compatible.
%     For internal use only.

%   Copyright 2019 The MathWorks, Inc.

% Inputs need to match in size without generalized scalar expansion.
varargout = varargin;
[errorcode, varargout{:}] = distchck(length(varargin),varargin{:});
if errorcode > 0
    throwAsCaller(MException(message('stats:normcdf:InputSizeMismatch')));
end
end
