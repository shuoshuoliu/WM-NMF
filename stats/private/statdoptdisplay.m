function [f,settry,setiter] = statdoptdisplay(algtxt)
%STATDOPTDISPLAY Helper to display status of D-optimal algorithm
%   [F,SETTRY,SETITER]=STATDOPTDISPLAY(ALGTXT) sets up a message window
%   with the algorithm described by ALGTXT, and returns a figure handle F,
%   a function handle SETTRY for setting the "try" number, and a function
%   handle SETITER for setting the iteration number.

%   Copyright 2005-2011 The MathWorks, Inc.


% Set up figure to hold status display
textprops = {'FontName','fixedwidth','FontSize',12,'FontWeight','bold'};
screen = get(0,'ScreenSize');
f = figure('Units','Pixels', ...
    'Menubar','none',...
    'Position',[25 screen(4)-150 300 60],'HandleVisibility','callback',...
    'NumberTitle','off', ...
    'Name',getString(message('stats:statdoptdisplay:figName_DOptimalDesignGeneration')),...
    'WindowStyle', 'Normal', ...
    'DockControls', 'off');
ax = axes('Parent',f,'Visible','off');
t = text(0,0.8, ...
    sprintf('%s',getString(message('stats:statdoptdisplay:alignedText_Algorithm',algtxt))), ...
    textprops{:}, ...
    'Parent',ax);
a = get(t,'extent');
set(f,'Position', [25,screen(4)-150,300*max(1,a(3)), 90], 'Visible', 'On');

% Set up text strings to be filled in by nested function
t1 = text(0,0.5,'','Parent',ax,textprops{:});
t2 = text(0,0.2,'','Parent',ax,textprops{:});

% Start things off
settryfunc(1);
setiterfunc(1);

% Define return values
settry = @settryfunc;
setiter = @setiterfunc;

function settryfunc(ntry)
    % Nested function sets try number
    if ishghandle(t1)
        set(t1,'String',sprintf('%s',getString(message('stats:statdoptdisplay:alignedText_Try',ntry))));
        drawnow;
    end
    end

function setiterfunc(niter)
    % Nested function sets iteration number
    if ishghandle(t2)
        set(t2,'String',sprintf('%s',getString(message('stats:statdoptdisplay:alignedText_Iteration',niter))));
        drawnow;
    end
    end

end
