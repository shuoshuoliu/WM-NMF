function h = refcurve(varargin)
%REFCURVE Add a reference curve (polynomial) to a plot.
%   REFCURVE(P) adds a line with the given polynomial coefficients 
%   to the current figure.
%   REFCURVE(AX,...) add a line into AX.
%   H = REFCURVE(P) returns the handle to the line object
%   in H.
%   REFCURVE with no arguments plots the function Y = 0.
%   Example: 
%       y = p(1)*x^d + p(2)*x^(d-1) + ... + p(d)*x + p(d+1)
%       Shows a polynomial of degree, d.
%       Note that p(1) goes with the highest order term.
%
%   See also POLYFIT, POLYVAL.   

%   B.A. Jones 2-3-95
%   Copyright 1993-2019 The MathWorks, Inc.


if isa(varargin{1},'matlab.graphics.axis.Axes') ||...
        isa(varargin{1},'matlab.ui.control.UIAxes')
    if nargin == 1
        p = zeros(2,1);
    else
        p = varargin{2};
    end
    cax = varargin{1};
    fig = ancestor(cax,'figure');
else
    if nargin == 0
        p = zeros(2,1);
    else
        p = varargin{1};
        cax = gca;
        fig = gcf;
    end
end


xlimits = get(cax,'Xlim');
ylimits = get(cax,'Ylim');

np = get(fig,'NextPlot');
set(fig,'NextPlot','add');

xdat = linspace(xlimits(1),xlimits(2),100);
ydat = polyval(p,xdat);
maxy = max(ydat);
miny = min(ydat);

if maxy > ylimits(2)
  if miny < ylimits(1)
     set(cax,'YLim',[miny maxy]);
  else
     set(cax,'YLim',[ylimits(1) maxy]);
  end
else
  if miny < ylimits(1)
     set(cax,'YLim',[miny ylimits(2)]);
  end
end

if nargout == 1
   h = line(cax,xdat,ydat);
   set(h,'LineStyle','-');
else
   hh = line(cax,xdat,ydat);
   set(hh,'LineStyle','-');
end

set(fig,'NextPlot',np);
