function dupfigure(dffig)
%DFDUPFIGURE Make a duplicate, editable copy of the distribution fitter figure


%   Copyright 2003-2020 The MathWorks, Inc.

% Copy the regular axes, not the legend axes
f = uifigure;
ax = findall(dffig,'Type','axes','Tag','dfMainAxes');
copyobj(ax,f);
newax = findall(f,'Type','axes','Tag','dfMainAxes');

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
    legpos = stats.internal.dfit.getrelativelegendposition(dffig,ax,legh);
    if ~isempty(legpos)
        leginfo = stats.internal.dfit.getlegendinfo(legh);
        newlegh = legend(newax,h1,txt,leginfo{:});
        stats.internal.dfit.setrelativelegendposition(legpos,f,newax,newlegh);
        set(newlegh,'Interpreter','none');
    end
end
