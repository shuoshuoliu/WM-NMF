function varargout = plotGroupedKSDensity(x,group,varargin)
%PLOTGROUPEDKSDENSITY plots grouped kernel smooth densities
%   plotGroupedKSDensity(X,GROUP) plots kernel smooth densities in each
%   group of the data X. X must be a vecotr. GROUP is a grouping variable
%   that can be accepted by GROUP2IDX. Missing data, NaN elments are
%   excluded for computation. No error checking.
%   
%   [XRANGE, PX] = plotGroupedKSDensity(X,GROUP,VARARGIN) computes the HG
%   axis 'xdata' and 'ydata' in each group and plot the lines. XRANGE
%   stores the xdata (a matrix) and PX stores the corresponding ydata (the
%   probablities).
%
%   Optional name/value pairs:
%     Parameter:       Value:
%      'AxisHandle'    A handle to the axis where the plot will be draw on.
%
%      'Color'         A color matrix that specify the color for each
%                      group. The number of rows must match the total
%                      number of groups, nG. The default is line(nG).
%
%      'LineWidth'     A vector that specify the line width property for
%                      each group. Vector length must be nG. The default is
%                      ones(1,nG).
%
%      'LineStyle'     A cell array of strings that specify the line style
%                      property for each group. It must be a vector of
%                      length nG. The default is all solid lines.
%
%      'Width'         Bandwidth of the kernel-smoothing window. The
%                      default is optimal for estimating normal densities,
%                      but you may want to choose a smaller value to reveal
%                      features such as multiple modes. It must be empty or
%                      a vector of length nG.
%
%      'Tag'           A string that is used to identify the plot. Default
%                      is 'groupedksplot'.
%
%      'Orientation'   Orientation of kernel density lines, specified as
%                      'vertical' or 'horizontal'. Default is vertical.
%
%      'AxisOn'        A logical value to determine whether the axis is
%                      'on' or not. If 'AxisOn' is false, a black
%                      horizontal line will be drawn as the X axis. The
%                      default is true.
% 
%      'Normalize'     A logical value to determine whether the density
%                      values are normalized or not. If 'Normalize' is
%                      true, all values are divided by the max range among
%                      all groups. The default is false.
%       
%   Example:
%   load fisheriris
%   [xrang, px] = internal.stats.plotGroupedKSDensity(meas(:,1),species);
%   internal.stats.plotGroupedKSDensity(meas(:,1),species,...
%            'color','bcr','LineWidth',[2 1 2],'LineStyle',{'-',':','-.'});
%   

%   Copyright 2012-2015 The MathWorks, Inc.

if isempty(group)
    group = ones(1,numel(x));
end
[grpID,gn] = grp2idx(group);
grp   = unique(grpID); % unique integer group labels
nGrp  = numel(gn); % total number of groups

x = x(:);
% Remove missing data
wasNaN = isnan(x) | isnan(grpID);
x(wasNaN) = [];
if isempty(x)
    return;
end
grpID(wasNaN) = [];

paramNames = {'AxisHandle','Color','LineWidth', 'LineStyle', 'Width', 'Tag', 'Orientation','AxisOn','Normalize'};
defaults   = {[], hsv(nGrp), ones(1,nGrp), [], [], 'groupedksplot','vertical',false,false};

[h,clr,lw,ls,ww,tag,or,isAxisOn,isNormalized] = internal.stats.parseArgs(paramNames,defaults,varargin{:});
if isempty(h)
   h = gca; 
end

cxmax = max(x) ;
cxmin = min(x) ;

if cxmax == cxmin
    [~, xrange]=ksdensity(x);
    xLim = [min(xrange),max(xrange)];
else
    dx = 0.1*range(x) ;
    xLim = [cxmin-dx, cxmax+dx];
    xrange = xLim(1):0.01*dx: xLim(2);
end
px = zeros(nGrp,size(xrange,2));

for i = 1:nGrp
    xg = x(grpID == grp(i));
    if ~isempty(xg)        
        if isempty(ww)
            px(i,:) = ksdensity(xg,xrange);
        else
            px(i,:) = ksdensity(xg,xrange,'Width',ww(i));
        end
    else
        px(i,:) = nan(1,size(xrange,2));
    end
end

if isNormalized
    px = px/max(range(px,2));
end

if isempty(clr)
    clr = lines(nGrp);
elseif ischar(clr) && isvector(clr)
    clr = clr(:);
end

% Now draw the kernel density line of each group
if strcmpi(or,'horizontal')
    hXLines = plot(h,px,xrange);
else
    hXLines = plot(h, xrange,px);
end
set(hXLines,'Tag',tag);
% Set the line properties accordingly
for i = 1:nGrp
    set(hXLines(i),'Color',clr(i,:),'LineWidth',lw(i));
    if ~isempty(ls)
        set(hXLines(i),'LineStyle',ls{i});
    end
end
% Draw a black horizontal line as the X axis because axis may be 'off'.
if ~isAxisOn
    if strcmpi(or,'vertical')
        line(xLim,[0 0],'Color','k');
    else
        line([0 0],xLim,'Color','k');
    end
end

if nargout > 0
    varargout{1} = xrange;
    varargout{2} = px;
end

end
