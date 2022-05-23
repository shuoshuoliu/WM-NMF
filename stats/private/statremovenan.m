function [badin,wasnan,varargout]=statremovenan(varargin)
%STATREMOVENAN Remove NaN values from inputs
%
%  See also internal.stats.insertnan.

%   Copyright 1993-2019 The MathWorks, Inc.


[badin,wasnan,varargout{1:nargout-2}] = internal.stats.removenan(varargin{:});
