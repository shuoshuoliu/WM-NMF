function name = getfitname()
%GETFITNAME returns a new distributionFitter fit name.

% Copyright 2003-2020 The MathWorks, Inc.

count=dfgetset('fitcount');
if isempty(count)
    count = 1;
end
taken = 1;
while taken
    name=getString(message('stats:dfstrings:name_FitNumber', count));
    if isempty(find(stats.internal.dfit.getfitdb,'name',name))
        taken = 0;
    else
        count=count+1;
    end
end
dfgetset('fitcount',count+1);

