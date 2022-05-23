function dfdupfigure(dffig)
%DFDUPFIGURE Make a duplicate, editable copy of the distribution fitter figure


%   Copyright 2003-2014 The MathWorks, Inc.

% Copy the regular axes, not the legend axes
f = figure;
ax = findall(dffig,'Type','axes','Tag','main');
copyobj(ax,f);
newax = findall(f,'Type','axes','Tag','main');

% Adjust layout in new figure, but don't add axis controls
dfadjustlayout(f,'off');


% Remove any context menus and callbacks associated with the old figure
set(findall(newax,'Type','line'),...
    'DeleteFcn','','UIContextMenu',[],'ButtonDownFcn','');

% Make a new legend based on the original, if any
legh = get(ax,'Legend');

% Recover information as to which items in the plot are to
% appear in the legend, and with what associated text.
contents = dfgetset('legendContents');
if isempty(contents) || isempty(contents{1})
    return;
end
h0 = contents{1};
txt = contents{2};

if ~isempty(h0)
    c0 = get(ax,'Child');
    c1 = get(newax,'Child');
    h1 = h0;
    remove = false(size(h1));
    for j=1:length(h0)
        k = find(c0==h0(j));
        if isempty(k)
            remove(j) = true;
        else
            % Convert to lineseries
            h1(j) = internal.stats.line2chartline(c1(k(1)));
        end
    end
    h1(remove) = [];
    txt(remove) = [];
    legpos = getrelativelegendposition(dffig,ax,legh);
    if ~isempty(legpos)
        leginfo = dfgetlegendinfo(legh);
        newlegh = legend(newax,h1,txt,leginfo{:});
        setrelativelegendposition(legpos,f,newax,newlegh);
        set(newlegh,'Interpreter','none');
    end
end
