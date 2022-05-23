function disp(a)
%PIECEWISEDISTRIBUTION/DISP Display a piecewisedistribution object.
%   DISP(A) prints a text representation of the piecewisedistribution
%   object A, without printing the object name.  In all other ways it's
%   the same as leaving the semicolon off an expression.
%
%   See also PIECEWISEDISTRIBUTION.

%   Copyright 2006-2011 The MathWorks, Inc.



isLoose = strcmp(get(0,'FormatSpacing'),'loose');

if (isLoose)
   fprintf('\n');
end

% Create some text in preparation for display
nseg = nsegments(a);
QI = [[-Inf; a.Q], [a.Q; Inf]];  % intervals on quantile scale
PI = [[0;a.P], [a.P;1]];          % intervals on probability scale

C1 = cell(nseg,1);
C2 = cell(nseg,1);
for j=1:nseg
    C1{j} = sprintf('%g < x < %g',QI(j,:));
    C2{j} = sprintf('(%g < p < %g)',PI(j,:));
end

% Measure the text
max1 = max(cellfun(@length,C1));
max2 = max(cellfun(@length,C2));

% Display with uniform widths
if nseg==1
    fprintf('%s',getString(message('stats:piecewisedistribution:sprintf_Piecewise1Segment')));
else
    fprintf('%s',getString(message('stats:piecewisedistribution:sprintf_PiecewiseNSegments',nseg)));
end

for j=1:nseg
    fprintf('   %*s  %*s: %s\n',...
            max1,C1{j},max2,C2{j}, a.distribution(j).description);
end

if (isLoose)
   fprintf('\n');
end
