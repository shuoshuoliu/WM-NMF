function hg = plotGroupedHist(x,group,varargin)  
%PLOTGROUPEDHIST plots grouped histogram
%   plotGroupedHist(X,GROUP) plots histogram in each group of the data X. 
%   X must be a vector. GROUP is a grouping variable that can be accepted 
%   by GROUP2IDX. Missing data, NaN elments are excluded for computation. 
%   No error checking.
%
%   Optional name/value pairs:
%     Parameter:       Value:
%      'AxisHandle'    A handle to the axis where the plot will be draw on.
%      'Color'         A color matrix that specify the color for each
%                      group. The number of rows must match the total
%                      number of groups, nG. The default is line(nG).
%      'AxisOn'        A logical value to determine whether the axis is
%                      'on' or not. If 'AxisOn' is false, a black
%                      horizontal line will be drawn as the X axis. The
%                      default is true.
%      'NBins'         A positive integer value, a two-element vector, a
%                      1-by-1, 1-by-2, or a 2-by-1 cell array, specifying
%                      the number of bins for each group in the X and Y
%                      histograms.  All numbers should be positive integers
%                      greater than or equal to 2.
%      'PlotGroup'     A logical value indicating if grouped histograms or
%                      grouped kernel density plots are created when a
%                      grouping variable is given. The default is true. The
%                      histogram of the whole data set will be created if
%                      it is set to false.
%
%   H = PLOTGROUPEDHIST(...) returns an array of handles to the histograms
%   of each group.
%
%   Example:
%   load fisheriris
%   internal.stats.plotGroupedHist(meas(:,1),species,'color','bcr');
%   hg = internal.stats.plotGroupedHist(meas(:,1),species);

%   Copyright 2014-2015 The MathWorks, Inc.

hg = statslib.internal.plotGroupedHist(x,group,varargin{:});

end




