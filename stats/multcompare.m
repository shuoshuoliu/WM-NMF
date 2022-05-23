function [comparison,means,h,gnames] = multcompare(stats,varargin)
%MULTCOMPARE Perform a multiple comparison of means or other estimates
%   MULTCOMPARE performs a multiple comparison using one-way anova or
%   anocova results to determine which estimates (such as means,
%   slopes, or intercepts) are significantly different.
%
%   COMPARISON = MULTCOMPARE(STATS) performs a multiple comparison
%   using a STATS structure that is obtained as output from any of
%   the following functions:  anova1, anova2, anovan, aoctool,
%   kruskalwallis, friedman.  The return value COMPARISON is a matrix
%   with one row per comparison and six columns.  Columns 1-2 are the
%   indices of the two samples being compared.  Columns 3-5 are a lower
%   bound, estimate, and upper bound for their difference. Column 6 is the 
%   p-value for each individual comparison. 
%
%   COMPARISON = MULTCOMPARE(STATS, 'PARAM1',val1, 'PARAM2',val2,...)
%   specifies one or more of the following name/value pairs:
%   
%     'alpha'       Specifies the confidence level as 100*(1-ALPHA)%
%                   (default 0.05).
%     'display'     Either 'on' (the default) to display a graph of the
%                   estimates with comparison intervals around them, or
%                   'off' to omit the graph.
%     'ctype'       The type of critical value to use.  Choices are
%                   'tukey-kramer' (default), 'dunn-sidak', 'bonferroni',
%                   'scheffe'.  Enter two or more choices separated by
%                   spaces, to use the minimum of those critical values.
%     'dimension'   A vector specifying the dimension or dimensions over
%                   which the population marginal means are to be
%                   calculated.  Used only if STATS comes from anovan.
%                   The default is to compute over the first dimension
%                   associated with a categorical (non-continuous) factor.
%                   The value [1 3], for example, computes the population
%                   marginal mean for each combination of the first and
%                   third predictor values.
%     'estimate'    Estimate to compare.  Choices depend on the source of
%                   the stats structure:
%         anova1:  ignored, compare group means
%         anova2:  'column' (default) or 'row' means
%         anovan:  ignored, compare population marginal means
%         aoctool:  'slope', 'intercept', or 'pmm' (default is 'slope'
%                   for separate-slopes models, 'intercept' otherwise)
%         kruskalwallis:  ignored, compare average ranks of columns
%         friedman:  ignored, compare average ranks of columns
%
%   [COMPARISON,MEANS,H,GNAMES] = MULTCOMPARE(...) returns additional
%   outputs.  MEANS is a matrix with columns equal to the estimates
%   and their standard errors.  H is a handle to the figure containing
%   the graph.  GNAMES is a cell array with one row for each group,
%   containing the names of the groups.
%
%   The intervals shown in the graph are computed so that to a very close
%   approximation, two estimates being compared are significantly different
%   if their intervals are disjoint, and are not significantly different if
%   their intervals overlap.  (This is exact for multiple comparison
%   of means from anova1, if all means are based on the same sample size.)
%   You can click on any estimate to see which means are significantly
%   different from it.
%
%   Two additional CTYPE choices are available.  The 'hsd' option stands
%   for "honestly significant differences" and is the same as the
%   'tukey-kramer' option.  The 'lsd' option stands for "least significant
%   difference" and uses plain t-tests; it provides no protection against
%   the multiple comparison problem unless it follows a preliminary overall
%   test such as an F test.
%
%   MULTCOMPARE does not support multiple comparisons using anovan output
%   for a model that includes random or nested effects.  The calculations
%   for a random effects model produce a warning that all effects are
%   treated as fixed.  Nested models are not accepted.
%
%   Example:  Perform 1-way anova, and display group means with their names
%
%      load carsmall
%      [p,t,st] = anova1(MPG,Origin,'off');
%      [c,m,h,nms] = multcompare(st,'display','off');
%      [nms num2cell(m)]   
%
%   See also ANOVA1, ANOVA2, ANOVAN, AOCTOOL, FRIEDMAN, KRUSKALWALLIS.

%   Also supports older calling sequence:
%      [...] = MULTCOMPARE(STATS,ALPHA,DISPLAYOPT,CTYPE,ESTIMATE,DIM)
%
%   Reference: Y. Hochberg and A.C. Tamhane, "Multiple Comparison
%   Procedures," Wiley, New York, 1987.
%
%   The Tukey-Kramer critical value is the default.  This is known
%   to be the best choice for one-way anova comparisons, but the
%   conjecture that this is best for other comparisons is
%   unproven.  The Bonferroni and Scheffe critical values are
%   conservative in all cases.

%   Copyright 1993-2012 The MathWorks, Inc.
%   Copyright 1993-2013 The MathWorks, Inc.
%   Copyright 1993-2014 The MathWorks, Inc.


% Check inputs and assign default values
if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

narginchk(1,Inf);

if isempty(varargin) || (~isempty(varargin{1}) && ischar(varargin{1}))
   okargs =   {'alpha' 'displayopt' 'ctype' 'estimate' 'dimension'};
   defaults = {0.05    'on'         ''      ''         []};
   [alpha,displayopt,ctype,estimate,dim] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});
else
   % Old-style calling sequence with fixed argument positions
   if (nargin >= 2) && ~isempty(varargin{1})
      alpha = varargin{1};
   else
      alpha = 0.05;
   end
   if (nargin >= 3) && ~isempty(varargin{2})
      displayopt = varargin{2};
   else
      displayopt = 'on';
   end
   if (nargin>=4) && ~isempty(varargin{3})
      ctype = varargin{3};
   else
      ctype = '';
   end
   if (nargin>=5) && ~isempty(varargin{4})
      estimate = varargin{4};
   else
      estimate = '';
   end
   if (nargin>=6) && ~isempty(varargin{5})
      dim = varargin{5};
   else
      dim = [];
   end
end

if ~isstruct(stats) || ~isfield(stats,'source')
   error(message('stats:multcompare:BadStats'));
end
if (length(alpha)~=1 || ~isfinite(alpha) || alpha<=0 || alpha>=1)
   error(message('stats:multcompare:BadAlpha'));
end
if ~(isequal(displayopt, 'on') || isequal(displayopt, 'off'))
   error(message('stats:multcompare:BadDisplayOpt'));
end
dodisp = isequal(displayopt, 'on');

source = stats.source;

switch(source)
 case 'anova1'
   mname = getString(message('stats:multcompare:mnameStringAnova1'));
   gmeans = stats.means(:);
   gnames = stats.gnames;
   n = stats.n(:);
   df = stats.df;
   s = stats.s;
   ng = sum(n>0);
   if (df < 1)
      error(message('stats:multcompare:NotEnoughDataANOVA'));
   end
   
   gcov = diag((s^2)./n);
   [mn,se] = tValue(gmeans,gcov);
   t = mn./se;
   
   % Get Tukey-Kramer critical value
   if (isempty(ctype)), ctype = 'tukey-kramer'; end
   [crit,pval] = getcrit(ctype, alpha, df, ng, t);

 case 'anova2'
   docols = true;
   if (~isempty(estimate))
      estimate = internal.stats.getParamVal(estimate,{'row' 'column'},'ESTIMATE');
      docols = isequal(estimate,'column');
   end
   if (docols)
      gmeans = stats.colmeans(:);
      n = stats.coln(:);
      mname = getString(message('stats:multcompare:mnameStringAnova2Col'));
   else
      gmeans = stats.rowmeans(:);
      n = stats.rown(:);
      mname = getString(message('stats:multcompare:mnameStringAnova2Row'));
   end
   ng = length(gmeans);
   sigma = sqrt(stats.sigmasq);
   gnames = strjust(num2str((1:ng)'), 'left');
   df = stats.df;
   if (df < 1)
      error(message('stats:multcompare:NotEnoughDataANOVA'));
   end
   
   gcov = ((sigma^2)/n) * eye(ng);
   [mn,se] = tValue(gmeans,gcov);
   t = mn./se;
   
   % Get Tukey-Kramer critical value
   if (isempty(ctype)), ctype = 'tukey-kramer'; end
   [crit,pval] = getcrit(ctype, alpha, df, ng, t); 

   % This whole activity is a little strange if the model includes
   % interactions, especially if they are important.
   if (stats.inter && dodisp)     % model included an interaction term
      if (stats.pval < alpha)
         theString = getString(message('stats:multcompare:NoteSignifInteraction'));
      else
         theString = getString(message('stats:multcompare:NoteInsignifInteraction'));
      end
      theCell = textwrap({theString},80);
      fprintf('%s\n',theCell{:});
   end

 case 'anovan'
   mname = getString(message('stats:multcompare:mnameStringAnovan'));
   
   % We do not handle nested models
   try
       vnested = stats.vnested;
   catch me %#ok<NASGU>
       vnested = [];
   end
   if any(vnested(:))
       error(message('stats:multcompare:NoNesting'));
   end


   % Our calculations treat all effects as fixed
   if ~isempty(stats.ems)
      warning(message('stats:multcompare:IgnoringRandomEffects'))
   end
   
   nvars = length(stats.nlevels);
   P0 = stats.nullproject;
   
   % Make sure DIM is a scalar or vector of factor numbers.
   if isempty(dim)
       dim = find(stats.nlevels>1,1);
   end
   dim = dim(:);
   if isempty(dim) || any(dim<1 | dim>nvars | dim~=round(dim))
      error(message('stats:multcompare:BadDim', nvars));
   end
   dim = sort(dim);
   dim(diff(dim)==0) = [];
   if any(stats.nlevels(dim)<2)
      error(message('stats:multcompare:DimSpecifiesZeroDFFactor'));
   end
   
   % Create all combinations of the specified factors
   try
       continuous = stats.continuous;
   catch me %#ok<NASGU>
       continuous = zeros(nvars,1);
   end
   ffdesign = fullfact(stats.nlevels(dim));
   nrows = size(ffdesign,1);
   
   % Create a design matrix for these combinations
   dums = cell(nvars, 1);
   for j=1:length(dim)
      dj = dim(j);
      dums{dj} = idummy(ffdesign(:,j),3);  % dummy variables for each factor
   end
   
   % Create a vector of average values for remaining factors
   for j=1:nvars
      if isempty(dums{j});
          if continuous(j)
             dums{j} = stats.vmeans(j) * ones(nrows,1);
          else
             nlev = stats.nlevels(j);
             dums{j} = (1/nlev) * ones(nrows,nlev);
          end
      end
   end

   % Fill in x columns for each term
   termcols = stats.termcols(:);
   termstart = cumsum([1; termcols]);
   terms = [zeros(1,nvars); stats.terms];
   ncols = sum(termcols);
   x = zeros(size(ffdesign,1), ncols);
   x(:,1) = 1;
   for j=1:length(termcols)
      tm = terms(j,:);
      t0 = termstart(j);
      t1 = termstart(j) + termcols(j) - 1;
      if all(tm==0)
         x(:,t0:t1) = 1;
      else
         x0 = [];
         for k=nvars:-1:1
            if tm(k)
               x0 = termcross(x0,dums{k});
            end
         end
         x(:,t0:t1) = x0;
      end
   end

   % Compute estimates and their standard errors
   gmeans = x * stats.coeffs;
   xproj = (x*P0)';
   tmp = stats.Rtr \ xproj;
   if (stats.dfe == 0)
      mse = NaN;
   else
      mse = max(stats.mse,0);
   end
   gcov = mse * (tmp' * tmp);
   
   % Find non-estimable means and set them to NaN
   Xrows = stats.rowbasis';           % row basis of original X matrix
   bb = Xrows \ (x');                 % fit rows of x to row basis
   xres = Xrows * bb - x';            % get residuals
   xres = sum(abs(xres));             % sum of absolute residuals
   cutoff = sqrt(eps(class(xres))) * size(xres,2); % cutoff for large residuals
   gmeans(xres > cutoff) = NaN;       % not in row space of original X
   
   [mn,se] = tValue(gmeans,gcov);
   t = mn./se;
   
   % Get Tukey-Kramer critical value
   if (isempty(ctype)), ctype = 'tukey-kramer'; end
   [crit,pval] = getcrit(ctype, alpha, stats.dfe, length(gmeans),t);

   % Get names for each group
   ngroups = size(ffdesign,1);
   gnames = cell(ngroups,1);
   allnames = stats.grpnames;
   varnames = stats.varnames;
   for j=1:ngroups
      v1 = dim(1);
      vals1 = allnames{v1};
      nm = sprintf('%s=%s',varnames{v1},vals1{ffdesign(j,1)});
      for i=2:size(ffdesign,2)
         v2 = dim(i);
         vals2 = allnames{v2};
         nm = sprintf('%s,%s=%s',nm,varnames{v2},vals2{ffdesign(j,i)});
      end
      gnames{j} = nm;
   end

 case 'aoctool'
   model = stats.model;
   if (model==1 || model==3)
      error(message('stats:multcompare:NoMultipleParameters'));
   end
   gnames = stats.gnames;
   n = stats.n(:);
   ng = length(n);
   df = stats.df;
   if (df < 1)
      error(message('stats:multcompare:NotEnoughDataAOC'));
   end

   % Get either slope or intercept estimates and covariances
   if (isempty(estimate))
      if (model == 5)
         estimate = 'slope';
      else
         estimate = 'intercept';
      end
   else
      estimate = internal.stats.getParamVal(estimate,{'slope' 'intercept' 'pmm'},'ESTIMATE');
   end
   switch(estimate)
    case 'slope'
      if (~isfield(stats, 'slopes'))
         error(message('stats:multcompare:BadStatsNoSlope'));
      end
      gmeans = stats.slopes;
      gcov = stats.slopecov;
      mname = getString(message('stats:multcompare:mnameStringAoctoolSlope'));
    case 'intercept'
      if (~isfield(stats, 'intercepts'))
         error(message('stats:multcompare:BadStatsNoIntercept'));
      end
      gmeans = stats.intercepts;
      gcov = stats.intercov;
      mname = getString(message('stats:multcompare:mnameStringAoctoolIntercept'));
    case 'pmm'
      gmeans = stats.pmm;
      gcov = stats.pmmcov;
      mname = getString(message('stats:multcompare:mnameStringAoctoolPmm'));
   end

   if (any(any(isinf(gcov))))
      error(message('stats:multcompare:InfiniteVariance', mname));
   end
   
   [mn,se] = tValue(gmeans,gcov);
   t = mn./se;
   
   % Get Tukey-Kramer critical value
   if (isempty(ctype)), ctype = 'tukey-kramer'; end
   [crit,pval] = getcrit(ctype, alpha, df, ng, t);

 case 'kruskalwallis'
   gmeans = stats.meanranks(:);
   gnames = stats.gnames;
   n = stats.n(:);
   sumt = stats.sumt;
   ng = length(n);
   N = sum(n);
   mname = getString(message('stats:multcompare:mnameStringKruskalwallis'));

   gcov = diag(((N*(N+1)/12) - (sumt/(12*(N-1)))) ./ n);
   [mn,se] = tValue(gmeans,gcov);
   t = mn./se;
   
   % Get critical value; H&T recommend the Tukey-Kramer value
   if (isempty(ctype)), ctype = 'tukey-kramer'; end
   [crit,pval] = getcrit(ctype, alpha, Inf, ng, t);
   
   % Note that the intervals in M can be used for testing but not
   % for simultaneous confidence intervals.  See H&T, p. 249.
   if (dodisp)
      disp(getString(message('stats:multcompare:NoteNotSimul')));
   end

 case 'friedman'
   gmeans = stats.meanranks(:);
   n = stats.n;
   ng = length(gmeans);
   sigma = stats.sigma;
   mname = getString(message('stats:multcompare:mnameStringFriedman'));
   gnames = strjust(num2str((1:ng)'), 'left');

   gcov = ((sigma^2) / n) * eye(ng);
   [mn,se] = tValue(gmeans,gcov);
   t = mn./se;
   
   % Get critical value; H&T recommend the Tukey-Kramer value
   if (isempty(ctype)), ctype = 'tukey-kramer'; end
   [crit,pval] = getcrit(ctype, alpha, Inf, ng, t);

   % Note that the intervals in M can be used for testing but not
   % for simultaneous confidence intervals.  See H&T, p. 249.
   if (dodisp)
      disp(getString(message('stats:multcompare:NoteNotSimul')));
   end

 otherwise
   error(message('stats:multcompare:BadStats'));
end


% Create output matrix showing tests for all pairwise comparisons
% and graph that approximates these tests.
[M,MM,hh] = makeM(gmeans, gcov, mn, se, crit, gnames, mname, dodisp, pval);

comparison = M;
if (nargout>1), means = MM; end
if (nargout>2), h = hh; end

% -----------------------------------------------
function [crit,pval] = getcrit(ctype, alpha, df, ng, t)
% Get the minimum of the specified critical values
crit = Inf;
[onetype,ctype] = strtok(ctype);

while(~isempty(onetype))
   if (length(onetype) == 1)
      switch onetype
       case 't', onetype = 'tukey-kramer';
       case 'd', onetype = 'dunn-sidak';
       case 'b', onetype = 'bonferroni';
       case 's', onetype = 'scheffe';
       case 'h', onetype = 'tukey-kramer';
       case 'l', onetype = 'lsd';
      end
   end
   if (isequal(onetype, 'hsd')), onetype = 'tukey-kramer'; end
   
   switch onetype
    case 'tukey-kramer' % or hsd
     crit1 = internal.stats.stdrinv(1-alpha, df, ng) / sqrt(2);
     
     % The T-K algorithm is inaccurate for small alpha, so compute
     % an upper bound for it and make sure it's in range.
     ub = getcrit('dunn-sidak', alpha, df, ng, t);
     if (crit1 > ub), crit1 = ub; end
     pval = 0*t;
     for j=1:numel(t)
         pval(j) = 1 - stdrcdf(sqrt(2)*abs(t(j)),df,ng);
     end

    case 'dunn-sidak'
     kstar = nchoosek(ng, 2);
     alf = 1-(1-alpha).^(1/kstar);
     if (isinf(df))
        crit1 = norminv(1-alf/2);
     else
        crit1 = tinv(1-alf/2, df);
     end
     pval = 1 - (1-2*tcdf(-abs(t),df)).^kstar;

    case 'bonferroni'
     kstar = nchoosek(ng, 2);
     if (isinf(df))
        crit1 = norminv(1 - alpha / (2*kstar));
     else
        crit1 = tinv(1 - alpha / (2*kstar), df);
     end
     pval = 2*kstar*tcdf(-abs(t),df);

    case 'lsd'
     if (isinf(df))
        crit1 = norminv(1 - alpha / 2);
     else
        crit1 = tinv(1 - alpha / 2, df);
     end
     pval = 2*tcdf(-abs(t),df);

    case 'scheffe'
     if (isinf(df))
        tmp = chi2inv(1-alpha, ng-1) / (ng-1);
     else
        tmp = finv(1-alpha, ng-1, df);
     end
     crit1 = sqrt((ng-1) * tmp);
     pval = fcdf((t.^2)/(ng-1),ng-1,df,'upper');
     
    otherwise
     error(message('stats:multcompare:BadCType', ctype));
   end

   pval(pval>1) = 1;
   if (~isnan(crit1)), crit = min(crit, crit1); end
   [onetype,ctype] = strtok(ctype);
end

% -----------------------------------------------
function [M,MM,hh] = makeM(gmeans, gcov, mn, se, crit, gnames, mname, dodisp, pval)
% Create matrix to test differences, matrix of means, graph to display test

ng = length(gmeans);
MM = zeros(ng,2);
MM(:,1) = gmeans(:);
MM(:,2) = sqrt(diag(gcov));
MM(isnan(MM(:,1)),2) = NaN;

M = nchoosek(1:ng, 2);      % all pairs of group numbers
M(1,5) = 0;                 % expand M to proper size
g1 = M(:,1);
g2 = M(:,2);
M(:,4) = mn;
i12 = sub2ind(size(gcov), g1, g2);
delta = crit * se;
M(:,3) = M(:,4) - delta;
M(:,5) = M(:,4) + delta;
M(:,6) = pval;

% If requested, make a graph that approximates these tests
if (dodisp)
   % Find W values according to H&T (3.32, p. 98)
   d = zeros(ng, ng);
   d(i12) = se;
   sum1 = sum(sum(d));
   d = d + d';
   sum2 = sum(d);
   if (ng > 2)
      w = ((ng-1) * sum2 - sum1) ./ ((ng-1)*(ng-2));
   else
      w = sum1 * ones(2, 1) / 2;
   end
   halfwidth = crit * w(:);
   hh = meansgraph(gmeans, gmeans-halfwidth, gmeans+halfwidth, ...
                   gnames, mname);
   set(hh, 'Name', getString(message('stats:multcompare:MultcompareFigureTitleString', mname)));
else
   hh = [];
end

% -----------------------------------------------
function [mn,se]= tValue(gmeans,gcov)
% Make sure NaN groups don't affect other results
t = isnan(gmeans);
if any(t)
    gcov(t,:) = 0;
    gcov(:,t) = 0;
end
ng = length(gmeans);
M = nchoosek(1:ng, 2);      % all pairs of group numbers
g1 = M(:,1);
g2 = M(:,2);
mn = gmeans(g1) - gmeans(g2);
i12 = sub2ind(size(gcov), g1, g2);
gvar = diag(gcov);
se = sqrt(gvar(g1) + gvar(g2) - 2 * gcov(i12));              





