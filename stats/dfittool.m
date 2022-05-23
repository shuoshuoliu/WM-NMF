function varargout=dfittool(varargin)
%dfittool Open Distribution Fitter.
%   This function is obsolete. Use distributionFitter instead.
%
%   See also: distributionFitter.

if nargin > 0
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargout==0
    distributionFitter(varargin{:});
else
    [varargout{1:nargout}] = distributionFitter(varargin{:});
end

%   Copyright 2016 The MathWorks, Inc.
