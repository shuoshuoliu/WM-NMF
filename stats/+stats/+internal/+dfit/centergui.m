function centergui(fig)
%CENTERGUI Move a figure window to the center of the screen.

%   Copyright 2019-2020 The MathWorks, Inc.

% centergui was created from movegui, but has less flexibility; 
% "center" is the only option, and there is no drawnow.
% 
% centergui was created this because movegui calls drawnow, which can
% cause distributionFitter to get into states where events are triggered 
% before some objects have been completely contructed. If we fix 
% the general issue of events coming sooner than desired, we can delete 
% centergui and change calls to "movegui(f, 'center');

oldfunits = get(fig, 'Units');
set(fig, 'Units', 'pixels');

% save figure position before making adjustments
oldpos = fig.Position;

% estimated value of toolbar
toolBarEstimate = 24;

% we can't rely on outerposition to place the uifigure
% correctly.  use reasonable defaults and place using regular
% position. 

if isunix
    % reasonable defaults to calculate outer position in unix

    % border estimate for figure window
    borderEstimate = 0;
    % padding value to account backward compatibility
    paddingEstimate = 6;
    % width adjustment is border value plus padding value of window
    widthAdjustment = borderEstimate + paddingEstimate;
    % estimated value of titlebar
    titleBarEstimate = 24;
    % estimated value of menubar
    menuBarEstimate = 22;
else
    % reasonable defaults to calculate outer position in windows

    % border estimate for figure window
    borderEstimate = 8;
    % border value of both left and right side of window
    widthAdjustment = borderEstimate * 2;
    % estimated value of titlebar
    titleBarEstimate = 31;
    % estimated value of menubar
    menuBarEstimate = 22;
end

% estimate the outer position
heightAdjustment = titleBarEstimate + borderEstimate;

% check if the figure has a menubar
haveMenubar = ~isempty(findall(fig,'type','uimenu'));

% check if the figure has any toolbars 
numToolbars = length(findall(fig,'type','uitoolbar'));

if haveMenubar
    heightAdjustment = heightAdjustment + menuBarEstimate;
end

if numToolbars > 0
    heightAdjustment = heightAdjustment + toolBarEstimate * numToolbars;
end

oldpos(3) = oldpos(3) + widthAdjustment;
oldpos(4) = oldpos(4) + heightAdjustment;

fleft   = oldpos(1);
fbottom = oldpos(2);
fwidth  = oldpos(3);
fheight = oldpos(4);

old0units = get(0, 'Units');
set(0, 'Units', 'pixels');
screensize = get(0, 'ScreenSize');
monitors = get(0,'MonitorPositions');
set(0, 'Units', old0units);

% Determine which monitor contains at least one of the corners of the figure window
% We cycle through each monitor and check the four corners of the figure. Starting with bottom left, moving clockwise. 
% If any one of the corners is found to be within a particular monitor we break the search and that monitor is used as the reference screen size for further calculations. 
for k = 1:size(monitors,1)
    monitorPos = monitors(k,:);    
    if (((fleft > monitorPos(1)) && (fleft < monitorPos(1) + monitorPos(3)) && (fbottom > monitorPos(2)) && (fbottom < monitorPos(2) + monitorPos(4))) || ... % bottom left
        ((fleft > monitorPos(1)) && (fleft < monitorPos(1) + monitorPos(3)) && (fbottom + fheight > monitorPos(2)) && (fbottom + fheight < monitorPos(2) + monitorPos(4))) || ... % left top
        ((fleft + fwidth > monitorPos(1)) && (fleft + fwidth < monitorPos(1) + monitorPos(3)) && (fbottom + fheight > monitorPos(2)) && (fbottom + fheight < monitorPos(2) + monitorPos(4))) || ... % top right 
        ((fleft + fwidth > monitorPos(1)) && (fleft + fwidth < monitorPos(1) + monitorPos(3)) && (fbottom > monitorPos(2)) && (fbottom < monitorPos(2) + monitorPos(4)))) % bottom right
        screensize = monitorPos;
        break;
    end
end

sx = screensize(1);
sy = screensize(2);
swidth = screensize(3);
sheight = screensize(4);
% make sure the figure is not bigger than the screen size
fwidth = min(fwidth, swidth);
fheight = min(fheight, sheight);

% swidth - fwidth == remaining width
rwidth  = swidth-fwidth;

% sheight - fheight == remaining height
rheight = sheight-fheight;

newpos = [rwidth/2, rheight/2];
newpos = newpos + [sx + borderEstimate, sy + borderEstimate];

newpos(3:4) = [fwidth, fheight];

% remove width and height adjustments added above
newpos(3) = newpos(3) - widthAdjustment;
newpos(4) = newpos(4) - heightAdjustment;
set(fig, 'Position', newpos);

set(fig, 'Units', oldfunits);
