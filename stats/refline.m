function h = refline(varargin)
%REFLINE Add a reference line to a plot.
%   REFLINE(SLOPE,INTERCEPT) adds a line with the given SLOPE and
%   INTERCEPT to the current figure.
%
%   REFLINE(SLOPE) where SLOPE is a two element vector adds the line
%        y = SLOPE(2) + SLOPE(1)*x 
%   to the figure. (See POLYFIT.)
%
%   H = REFLINE(SLOPE,INTERCEPT) returns the handle to the line object
%   in H.
%    
%   H = REFLINE(AX,...) plots into AX instead of GCA.
%
%   REFLINE with no input arguments superimposes the least squares line on 
%   the plot based on points recognized by LSLINE.
%
%   See also POLYFIT, POLYVAL, LSLINE.   

%   Copyright 1993-2014 The MathWorks, Inc. 

[ax,args,nargin] = axescheck(varargin{:});

if nargin >2
    error(message('MATLAB:TooManyInputs'));
end
if isempty(ax)
    ax = gca;
end

if nargin == 0
   hh = lsline(ax);
   if nargout >0
       h = hh;
   end
   return;
end

if nargin == 1
   if max(size(args{:})) == 2
      slope = args{1}(1);
      intercept = args{1}(2);
   else
      slope = args{1};
      intercept = 0;
   end
end

if nargin == 2
    slope = args{1};
    intercept=args{2};
end

xlimits = get(ax,'Xlim');
ylimits = get(ax,'Ylim');

np = get(ancestor(ax,'Figure'),'NextPlot');
set(ancestor(ax,'Figure'),'NextPlot','add');

if all(isfinite(xlimits))
    xdat = xlimits;
else
    xdat = ax.DataSpace.XLim;
end
ydat = intercept + slope.*xdat;
maxy = max(ydat);
miny = min(ydat);

if maxy > ylimits(2)
  if miny < ylimits(1)
     set(ax,'YLim',[miny maxy]);
  else
     set(ax,'YLim',[ylimits(1) maxy]);
  end
else
  if miny < ylimits(1)
     set(ax,'YLim',[miny ylimits(2)]);
  end
end

if nargout == 1
   h = line(xdat,ydat,'Parent',ax);
   set(h,'LineStyle','-');
else
   hh = line(xdat,ydat,'Parent',ax);
   set(hh,'LineStyle','-');
end

set(ancestor(ax,'Figure'),'NextPlot',np);
