function h=vline(varargin)
%VLINE draws vertical lines
%   VLINE(X) draws one vertical line for each element in vector X
%   VLINE(AX,X) draws the lines to the axes specified in AX
%   VLINE(X,...) accepts HG param/value pairs for line objects
%   H=VLINE(...) returns a handle to each line created
%
%   Note:  Be sure to include the initial AX argument if there is any
%          chance that X could be a valid handle.

%   Copyright 2007-2014 The MathWorks, Inc.



if isempty(varargin) || (ishghandle(varargin{1}) && length(varargin)==1)
    error(message('stats:hline:NotEnoughArgs'));
end

if isscalar(varargin{1}) && ishghandle(varargin{1})
    ax=varargin{1};
    varargin=varargin(2:end);
else
    ax=gca;
end

x = varargin{1};
varargin=varargin(2:end);
hh = matlab.graphics.chart.primitive.ConstantLine('parent',ax,varargin{:}); 
hh.Value = x;
hh.changedependvar('x');
if nargout>0
    h=hh;
end
