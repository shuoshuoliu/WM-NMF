function [p,stats] = vartestn(x,group,varargin)
%VARTESTN Test for equal variances across multiple groups.
%   VARTESTN(X) performs Bartlett's test for equal variances for the
%   columns of the matrix X.  This is a test of the null hypothesis
%   that the columns of X come from normal distributions with the same
%   variance, against the alternative that they come from normal
%   distributions with different variances.  The result is a display
%   of a box plot of the groups, and a summary table of statistics.
%
%   VARTESTN(X,GROUP) requires a column vector X, and a GROUP argument that
%   is a categorical variable, vector, string array, or cell array of
%   strings with one row for each element of X.  X values corresponding to
%   the same value of GROUP are placed in the same group.  The function
%   tests for equal variances across groups.
%
%   VARTESTN treats NaNs as missing values, and ignores them.
%
%   P = VARTESTN(...) returns the p-value, i.e., the probability of
%   observing the given result, or one more extreme, by chance if the null
%   hypothesis of equal variances is true.  Small values of P cast doubt
%   on the validity of the null hypothesis.
%
%   [P,STATS] = VARTESTN(...) returns a structure with the following
%   fields:
%      'chistat' -- the value of the test statistic
%      'df'      -- the degrees of freedom of the test
%
%  [...] = VARTESTN(...,'PARAM1',val1,'PARAM2',val2,...) specifies one or
%  more of the following name/value pairs:
%
%       Parameter       Value
%       'display'       'on' to display a boxplot and table, or 'off' to
%                       omit these displays. Default 'on'.
%       'testtype'      One of the following strings to control the type
%                       of test to perform:
%          'Bartlett'          Bartlett's test (default)
%          'LeveneQuadratic'   Levene's test computed by performing anova
%                              on the squared deviations of the data values
%                              from their group means
%          'LeveneAbsolute'    Levene's test computed by performing anova
%                              on the absolute deviations of the data
%                              values from their group means
%          'BrownForsythe'     Brown-Forsythe test computed by performing
%                              anova on the absolute deviations of the data
%                              values from the group medians
%          'OBrien'            O'Brien's modification of Levene's test with
%                              W=0.5.
%
%   The classical 'Bartlett' test is sensitive to the assumption that the
%   distribution in each group is normal. The other test types are more
%   robust to non-normal distributions, especially ones prone to outliers.
%   For these tests, the STATS output structure has a field named 'fstat'
%   containing the test statistic, and 'df1' and 'df2' containing its
%   numerator and denominator degrees of freedom.
%
%   Example:  Does the variance of mileage measurements differ
%             significantly from one model year to another?
%      load carsmall
%      vartestn(MPG,Model_Year)
%
%   See also VARTEST, VARTEST2, ANOVA1.

%   Copyright 2005-2020 The MathWorks, Inc.


narginchk(1,Inf);
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if (nargin < 2), group = []; end

% Process remaining arguments
displayopt = '';
testtype = '';

if nargin>=3    
    % See if the 2nd argument appears to be the start of name/value pairs,
    % with the group argument having been omitted
    if (ischar(group) && size(group,1)==1 ...
           && ~isempty(strncmpi(group,{'displayopt','testtype'},max(1,length(group)))))
        varargin = [{group} varargin];
        group = [];
    end
 
   if isequal(varargin{1},'on')||isequal(varargin{1},'off')||isempty(varargin{1})
          % Old syntax
          %   VARTESTN(X,GROUP,DISPLAYOPT,TESTTYPE)
          displayopt = varargin{1};
          if nargin>=4 
              testtype = varargin{2};
          end 
            
    else
        okargs =   {'displayopt' 'testtype'};
        defaults = { ''       ''};
        [displayopt, testtype] = ...
                         internal.stats.parseArgs(okargs,defaults,varargin{:});
    end
end

if isempty(displayopt)
    displayopt = 'on';
elseif ~isequal(displayopt,'on') && ~isequal(displayopt,'off')
    error(message('stats:vartestn:BadDisplayOpt'));
end
dodisplay = isequal(displayopt,'on');

if isempty(testtype)
    dorobust = false;
else
    ttypes = {'robust';'classical';'LeveneAbsolute';'OBrien';'BrownForsythe';'LeveneQuadratic';'Bartlett'};
    [~,i] = internal.stats.getParamVal(testtype,ttypes,'''testtype''');
    dorobust = i;
end

% Error if no data
if isempty(x)
   error(message('stats:vartestn:NoData'));
end

% Convert group to cell array from character array, make it a column
if (ischar(group) && ~isempty(group)), group = cellstr(group); end
if (size(group, 1) == 1), group = group'; end

% If x is a matrix, convert to vector form.
[n,m] = size(x);
if isempty(group)
   x = x(:);
   group = reshape(repmat((1:m), n, 1), n*m, 1);
elseif iscolumn(x)
   if numel(group) ~= n
       error(message('stats:vartestn:InputSizeMismatch'));
   end
else
   error(message('stats:vartestn:BadGroup'));
end

% Get rid of NaN
t = isnan(x) | ismissing(group);
if any(t)
   x(t) = [];
   group(t,:) = [];
end

% Compute group summary statistics
[igroup,gname] = grp2idx(group);
[gmean,gsem,gcount] = grpstats(x,igroup);
df = gcount-1;                  % group degrees of freedom
gvar = gcount .* gsem.^2;       % group variances
sumdf = sum(df);
if sumdf>0
   vp = sum(df.*gvar) / sumdf;  % pooled variance
else
   vp = NaN;
end
k = length(df);

if dorobust==1 || dorobust==6 || dorobust==3
   % dorobust=1 or 6: Robust Levene's test using squares
   % dorobust=3: Levene's test using absolute values  
   
   % Remove single-point groups and center each group
   spgroups = find(gcount<2);
   t = ~ismember(igroup,spgroups);
   xc = x(t) - gmean(igroup(t));
   ngroups = length(gcount) - length(spgroups);

   % Now do the anova and extract results from the anova table
   if dorobust==1 || dorobust==6
       xtransf = xc.^2;
       testname = getString(message('stats:vartestn:Levene'));
   else
       xtransf = abs(xc);
       testname = getString(message('stats:vartestn:LeveneAbs'));
   end
   if ngroups>1
       [p,atab] = anova1(xtransf, igroup(t),'off');
       B = atab{2,5};                             % F statistic
       Bdf = [atab{2,3}, atab{3,3}];              % both d.f.
   else
       p = NaN;
       B = NaN;
       Bdf = [0, length(xc)-ngroups];
   end
   Bdftable = sprintf('%d, %d',Bdf(1),Bdf(2));   % as text for display
   statname = 'fstat';

   
elseif dorobust==4
   % O'Brien's test (W=0.5)
   spgroups = find(gcount<2);
   t = ~ismember(igroup,spgroups);
   xc = x(t) - gmean(igroup(t));
   xcs = xc.^2;
   W = 0.5;
   xcw = ((W+gcount(igroup(t))-2).*gcount(igroup(t)).*xcs-W.*(gcount(igroup(t))-1).*gvar(igroup(t)))./...
            ((gcount(igroup(t))-1).*(gcount(igroup(t))-2)) ;

   ngroups = length(gcount) - length(spgroups);

   if ngroups>1
      [p,atab] = anova1(xcw, igroup(t),'off');
      B = atab{2,5};                             % F statistic
      Bdf = [atab{2,3}, atab{3,3}];              % both d.f.
   else
      p = NaN;
      B = NaN;
      Bdf = [0, length(xcw)-ngroups];
   end
   Bdftable = sprintf('%d, %d',Bdf(1),Bdf(2));   % as text for display
   testname = getString(message('stats:vartestn:OBrien'));
   statname = 'fstat';


elseif dorobust==5 
   % Brown-Forsythe's test: using the absolute deviations from the group medians
   spgroups = find(gcount<2);
   t = ~ismember(igroup,spgroups);
   gmedian = grpstats(x,igroup,@median); % group medians
   xcbf = x(t) - gmedian(igroup(t));
   ngroups = length(gcount) - length(spgroups);

   if ngroups>1
      [p,atab] = anova1(abs(xcbf), igroup(t),'off');
      B = atab{2,5};                             % F statistic
      Bdf = [atab{2,3}, atab{3,3}];              % both d.f.
   else
      p = NaN;
      B = NaN;
      Bdf = [0, length(xcbf)-ngroups];
   end
   Bdftable = sprintf('%d, %d',Bdf(1),Bdf(2));   % as text for display
   testname = getString(message('stats:vartestn:BrownForsythe'));
   statname = 'fstat';
    
   
else  % dorobust = 2 or 7
   % Classical Bartlett's test
   Bdf = max(0, sum(df>0)-1);
   t = df>0;
   if Bdf>0 && sumdf>0
      B = log(vp) * sum(df) - sum(df(t).*log(gvar(t)));
      C = 1 + (sum(1./df(t)) - 1/sum(df))/(3*Bdf);
      B = B/C;
   else
      B = NaN;
   end
    
   p = chi2pval(B,Bdf);
   Bdftable = Bdf;
   testname = getString(message('stats:vartestn:Bartlett'));
   statname = 'chisqstat';
end

if dodisplay
   Table = cell(k+6,4);
   Table(1,:) = { getString(message('stats:vartestn:Group')) ...
                  getString(message('stats:vartestn:Count')) ...
                  getString(message('stats:vartestn:Mean')) ...
                  getString(message('stats:vartestn:StdDev')) };
   Table(2:k+1,1) = gname;
   Table(2:k+1,2) = num2cell(gcount);
   Table(2:k+1,3) = num2cell(gmean);
   Table(2:k+1,4) = num2cell(sqrt(gvar));

   Table{k+2,1} = getString(message('stats:vartestn:Pooled'));
   Table{k+2,2} = sum(gcount);
   Table{k+2,3} = sum(gcount.*gmean) / sum(gcount);
   Table{k+2,4} = sqrt(vp);

   Table{k+3,1} = ' ';
   Table{k+4,1} = testname;
   Table{k+4,2} = B;
   Table{k+5,1} = getString(message('stats:vartestn:Df'));
   Table{k+5,2} = Bdftable;
   Table{k+6,1} = getString(message('stats:vartestn:Pval'));
   Table{k+6,2} = p;
   
   tblfig = statdisptable(Table, ...
                          getString(message('stats:vartestn:figMainTitle')),...
                          getString(message('stats:vartestn:displayHeader')));
   set(tblfig,'tag','table');
end

% Create output stats structure if requested, used by MULTCOMPARE
if (nargout > 1)
   stats = struct(statname,B,'df',Bdf);
end

if ~dodisplay
   return;
end

% Make a new figure, and draw into it explicitly
f1 = figure('pos',get(gcf,'pos') + [0,-200,0,0],'tag','boxplot');
ax = axes('Parent',f1);
boxplot(ax,x,group);
