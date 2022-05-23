function boundwarn(h)
%BOUNDWARN Helper function for dfittool to warn about conf bounds

%   Copyright 2003-2011 The MathWorks, Inc.


if isa(h,'stats.dfdata')
   % Warn if the function type does not allow data set bounds
   ftype = dfgetset('ftype');
   if isequal(ftype,'pdf')
      pname = getString(message('stats:dfstrings:errormsg_Densityhistogram'));
   elseif isequal(ftype,'probplot')
      pname = getString(message('stats:dfstrings:errormsg_Probability'));
   elseif isequal(ftype,'icdf')
      pname = getString(message('stats:dfstrings:dropdown_Quantile'));
   else
      return;
   end
   wmsg = getString(message('stats:dfstrings:cellstr_NoDataConf', ...
                    pname,h.name));
   warndlg(wmsg,getString(message('stats:dfstrings:dlg_DistributionFittingWarning')),'modal');

else % isa(h,'stats.dffit')
   % Warn if the function type does not allow fit bounds
   if ~dfgetset('dobounds')
      ftype = dfgetset('ftype');
      if isequal(ftype,'pdf')
         pname = getString(message('stats:dfstrings:errormsg_Density'));
      elseif isequal(ftype,'probplot')
         pname = getString(message('stats:dfstrings:errormsg_Probability'));
      else
         return;
      end
      wmsg = getString(message('stats:dfstrings:cellstr_NoFitConf', ...
                       pname,h.name));
      warndlg(wmsg,getString(message('stats:dfstrings:dlg_DistributionFittingWarning')),'modal');
   end
end
