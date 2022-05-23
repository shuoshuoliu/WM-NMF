function [hAxis,hLines] = svmplotdata(x,group,theAxis)
% SVMPLOTDATA draw a plot for 2-D data in SVM functions

% Copyright 2004-2014 The MathWorks, Inc.



holdState = ishold;
if nargin == 2
    class1 = 'r+';
    class2 = 'g*';
else
    axes(theAxis);
    hold on;
    class1 = 'm+';
    class2 = 'c*';
end
Xp = x(group==1,:);
h1 = plot(Xp(:,1),Xp(:,2),class1);
hAxis = get(h1,'parent');
hold on
Xn =  x(group==-1,:);
h2 = plot(Xn(:,1),Xn(:,2),class2);
if isempty(hAxis)
    hAxis = get(h2,'parent');
end
%axis equal
drawnow
% reset hold state if it was off
if ~holdState
    hold off
end

% Concatenating an empty handle with a non-empty handle gives a single
% non-empty handle. Use cell array to ensure that two handles be always
% returned.
hLines = {h1,h2};
