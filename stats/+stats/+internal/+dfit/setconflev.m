function changed = setconflev(dffig,clev)
%SETCONFLEV Set confidence level for curve fitting


%   Copyright 2003-2011 The MathWorks, Inc.

% Get new value
oldlev = dfgetset('conflev');
if isempty(clev)
   ctxt = inputdlg({getString(message('stats:dfstrings:dlg_ConfidenceLevelinPercent'))},...
                   getString(message('stats:dfstrings:dlg_SetConfidenceLevel')),1,{num2str(100*oldlev)},...
                   'on');
   if isempty(ctxt)
      clev = oldlev;
   else
      ctxt = ctxt{1};
      clev = str2double(ctxt);
      if ~isfinite(clev) || ~isreal(clev) || clev<=0 || clev>=100
         errordlg(getString(message('stats:dfstrings:dlg_RevertBadConfLev',...
                                    ctxt,sprintf('%g',100*oldlev))),...
                 getString(message('stats:dfstrings:dlg_Error')),'modal');
         clev = oldlev;
      else
         clev = clev/100;
      end
   end
end
changed = (oldlev~=clev);
if changed
   dfgetset('conflev',clev);
   
   % Update any existing data sets and fits
   dsdb = stats.internal.dfit.getdsdb;
   ds = down(dsdb);
   while(~isempty(ds))
      ds.confLev = clev;
      ds = right(ds);
   end
   fitdb = stats.internal.dfit.getfitdb;
   ft = down(fitdb);
   while(~isempty(ft))
      ft.confLev = clev;
      ft = right(ft);
   end
   stats.internal.dfit.updateylim;
end

% Check the appropriate menu item
h = [findall(dffig,'Type','uimenu','Tag','conflev90');
     findall(dffig,'Type','uimenu','Tag','conflev95');
     findall(dffig,'Type','uimenu','Tag','conflev99');
     findall(dffig,'Type','uimenu','Tag','conflevOther')];
set(h,'Checked','off');
verysmall = sqrt(eps);
if abs(clev-.90)<verysmall
   set(h(1),'Checked','on');
elseif abs(clev-.95)<verysmall
   set(h(2),'Checked','on');
elseif abs(clev-.99)<verysmall
   set(h(3),'Checked','on');
else
   set(h(4),'Checked','on');
end

