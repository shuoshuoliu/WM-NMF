function evec = getexcluded(ds,outlier)
%GETEXCLUDED Get an exclusion vector for this dataset/outlier combination


%   Copyright 2003-2020 The MathWorks, Inc.

ydata = ds.y;
evec = logical(false(size(ydata)));

if ~isequal(outlier,'(none)') && ~isempty(outlier)
   % For convenience, accept either an outlier name or a handle
   if ischar(outlier)
       s = settings;
       if (s.stats.DistributionFitter.LegacyDistributionFitter.ActiveValue == false)
           hExcl = find(stats.internal.dfit.getoutlierdb, 'name', outlier);
       else
           hExcl = dfswitchyard('dfgetexclusionrule',outlier);
       end
   else
       hExcl = outlier;
   end

   % Exclude points too low
   [bnd,strict] = getlowerbound(hExcl);
   if isfinite(bnd)
      if strict
         evec = evec | (ydata < bnd);
      else
         evec = evec | (ydata <= bnd);
      end
   end

   % Exclude points too high
   [bnd,strict] = getupperbound(hExcl);
   if isfinite(bnd)
      if strict
         evec = evec | (ydata > bnd);
      else
         evec = evec | (ydata >= bnd);
      end
   end
end
