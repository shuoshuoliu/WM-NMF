function clrRGB = colorStringToRGB(clrString)
%COLORSTRINGTORGB converts a color string or cell strings to RGB valued
%matrix. If clrString is a numerical matrix, COLORSTRINGTORGB checks if it
%has 3 columns. Each row of the output corresponds to each element of the
%input.
%
%   Examples:
%   >> clrRGB = internal.stats.colorStringToRGB('rgb')
%   clrRGB =
%      1     0     0
%      0     1     0
%      0     0     1
%   >> clrRGB = internal.stats.colorStringToRGB({'r','blue','m'})
%   clrRGB =
%      1     0     0
%      0     0     1
%      1     0     1
%   >> clrRGB = internal.stats.colorStringToRGB({'red'})
%   clrRGB = 1     0     0

%   Copyright 2012-2015 The MathWorks, Inc.


clrRGB = statslib.internal.colorStringToRGB(clrString);

end
