function h = meansgraph(gmeans, glower, gupper, gnames, mname)
%MEANSGRAPH Creates an interactive graph of several means
%   Used by MULTCOMPARE.  Not meant to be called directly.

%   Copyright 1993-2011 The MathWorks, Inc. 


% Check for call-back invocation
if (nargin == 1 && ~isempty(gcbo))
   [~,f] = gcbo;
   doclick(f, gmeans);
   return
end

% Make sure arguments are correct
if (nargin < 3)
   error(message('stats:meansgraph:TooFewInputs'));
end
if (min(size(gmeans))>1 || min(size(glower))>1 || min(size(gupper))>1)
   error(message('stats:meansgraph:VectorRequired'))
end
n = length(gmeans);
if (length(glower)~=n || length(gupper)~=n)
   error(message('stats:meansgraph:InputSizeMismatchLength'));
end

% Make sure gnames is a cell array of a single column.
if (nargin < 4)
   gnames = cellstr([repmat('Group ',n,1) num2str((1:n)')]);
elseif iscell(gnames)
   if (length(gnames)~=n || min(size(gnames))>1)
      error(message('stats:meansgraph:InputSizeMismatchGroupNames'));
   end
else
   if (size(gnames,1) ~= n)
      error(message('stats:meansgraph:InputSizeMismatchRows'));
   end
   gnames = cellstr(gnames);
end

if (nargin < 5), mname = 'means'; end

% Make comparison lines around the first group
x = [0.5 n+0.5];
x = [x NaN x];
y = [repmat(glower(1),1,2) NaN repmat(gupper(1),1,2)];
clf
plot(y,x,':', 'UserData','Comparison lines', 'Color',repmat(0.8,1,3));
title(getString(message('stats:meansgraph:ClickGroup')), 'fontweight','bold');
set(gca, 'YLim', [0.5, n+0.5]);

% Create initial plot with all groups
hold on
yci = zeros(2,1);
xci = zeros(2,1);
for j=1:n
   yci(1) = glower(j);
   yci(2) = gupper(j);
   xci(:) = n+1-j;
   ftxt = sprintf('meansgraph(%d)', j);
   plot(gmeans(j),n+1-j,'bo', 'UserData',j,  'ButtonDownFcn',ftxt);
   plot(yci,      xci,  'b-', 'UserData',-j, 'ButtonDownFcn',ftxt);
end

set(gca,'YTick',1:n, 'YTickLabel',gnames(n:-1:1,:), 'UserData',n, 'Box','on');
set(gcf,'UserData',mname);

% Color groups that are not significantly different
doclick(gcf,1);
hold off
if (nargout>0), h = gcf; end

function doclick(hfig,grp)
%DOCLICK Processes a click on group grp

% Get info from graph, return if there's an unexpected problem
mname = get(hfig,'UserData');
ax = get(hfig, 'CurrentAxes');
ngroups = get(ax, 'UserData');
gnames = get(ax, 'YTickLabel');
if (length(ngroups) ~= 1), return; end
if ((~isnumeric(ngroups)) || ngroups <= 0 || ngroups ~= round(ngroups))
   return
end

% Get bounds for this group
h = findobj(ax, 'UserData', -grp);
if (isempty(h)), return; end            % unexpected
x = get(h, 'XData');
U = max(x);
L = min(x);

% Change the comparison lines to use these values
h = findobj(ax, 'UserData', 'Comparison lines');
if (~isempty(h))
   x = get(h, 'XData');
   x(1:2) = L;
   x(4:5) = U;
   set(h, 'XData', x);
end

% Loop over all groups, adjusting colors
ndiff = 0;
for j=1:ngroups
   h = findobj(ax, 'UserData',-j, 'Type','line');
   if (isempty(h)), continue; end
   x = get(h, 'XData');
   if (j == grp)
      clr = 'b';
   elseif (max(x) < L || min(x) > U)
      ndiff = ndiff + 1;
      clr = 'r';
      diffgrp = j;
   else
      clr = repmat(0.8,1,3);
   end
   set(h, 'Color', clr);
   h = findobj(ax, 'UserData',j, 'Type','line');
   if (~isempty(h)), set(h, 'Color', clr); end
end

% Update text
grpname = deblank(gnames{ngroups + 1 - grp});
if (ndiff == 1)
   t = sprintf('%s',getString(message('stats:meansgraph:GroupsDifferent',...
               mname, grpname, deblank(gnames{ngroups + 1 - diffgrp}))));
else
   if (~isnan(str2double(grpname)))
      grpname = sprintf('%s',getString(message('stats:meansgraph:Group',grpname)));
   end
   if (ndiff == 0)
      t = sprintf('%s',getString(message('stats:meansgraph:NoGroupsSignificantlyDifferent', ...
                  mname, grpname)));
   else
      t = sprintf('%s',getString(message('stats:meansgraph:GroupsSignificantlyDifferent', ...
                  ndiff, mname, grpname)));
   end
end
xlab = get(ax, 'XLabel');
if (isempty(xlab))
   xlab = text('');
   set(ax, 'XLabel', xlab);
end
set(xlab, 'String', t);
