function updateplot(hFit,newlim)
%UPDATEPLOT Update the plot of this fit


%   Copyright 2003-2020 The MathWorks, Inc.

if isequal(hFit.fittype, 'smooth')
    if nargin==1
        updatesmoothplot(hFit);
    else
        updatesmoothplot(hFit,newlim);
    end
else
    if nargin==1
        updateparamplot(hFit);
    else
        updateparamplot(hFit,newlim);
    end
end
s = settings;
if (s.stats.DistributionFitter.LegacyDistributionFitter.ActiveValue == false) 
    stats.internal.dfit.updatelegend(dfgetset('dffig'));
else
    dfswitchyard('dfupdatelegend', dfgetset('dffig'));
end
