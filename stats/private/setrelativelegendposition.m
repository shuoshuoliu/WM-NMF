function setrelativelegendposition(relpos,fig,axh,legh)
%SETRELATIVELEGENDPOSITION Set (x,y) coords of legend center relative to axes


% Copyright 1993-2004 The MathWorks, Inc.

if nargin<2 || isempty(fig)
   fig = gcf;
end
if nargin<3 || isempty(axh)
   axh = get(fig,'CurrentAxes');
end
if nargin<4 || isempty(legh)
   legh = get(axh,'Legend');
end

if isempty(legh)
   return
end

if isnumeric(relpos) && numel(relpos)==2
   % Get the center of the legend in pixels
   oldu = get(legh,'units');
   legpos = get(legh,'position');
   legpos = hgconvertunits(fig,legpos,oldu,'pixels',fig);

   % Get the position of the axes in pixels
   oldu = get(axh,'units');
   axpos = get(axh,'position');
   axpos = hgconvertunits(fig,axpos,oldu,'pixels',fig);

   % Get the pixel coordinates of the relative position
   newctr = axpos(1:2) + relpos .* axpos(3:4);

   % Center the legend there
   legpos(1:2) = newctr - legpos(3:4)/2;
   legpos = hgconvertunits(fig,legpos,'pixels',oldu,fig);
   set(legh,'Position',legpos);
else
   set(legh,'Location',relpos);
end
