function txt = getresults(hFit)
%GETRESULTS Return text summarizing the results of the fit


%   Copyright 2003-2011 The MathWorks, Inc.

% Get text describing support of this fit
spt = hFit.support;
distspec = hFit.distspec;

% First see if inequalities are strict
if isempty(distspec)
   closedbound = [0 0];
else
   closedbound = distspec.closedbound;
end
ops = {'<' '<='};
op1 = ops{closedbound(1)+1};
op2 = ops{closedbound(2)+1};

if isequal(spt,'unbounded')
   spt = '-Inf < y < Inf';
elseif isequal(spt,'positive')
   spt = sprintf('0 %s y %s Inf', op1, op2);
else
   spt = sprintf('%g %s y %s %g', spt(1), op1, op2, spt(2));
end


% Get a cell array with lines that depend on the fit type
if isequal(hFit.fittype, 'smooth')
   t = getsmoothresults(hFit,spt);
else
   t = getparamresults(hFit,spt);
end

% Create a text string from this information
txt = '';
for j=1:length(t)
   a = t{j};
   for k=1:size(a,1)
      txt = sprintf('%s\n%s',txt,a(k,:));
   end
end

% Remove leading newline
nl = sprintf('\n');
if isequal(nl, txt(1:length(nl)))
   txt(1:length(nl)) = [];
end


% ------------------------------------
function t = getsmoothresults(hFit,spt)
%GETSMOOTHRESULTS Get results for smooth fit
% Fill cell array with some summary information
t = cell(1,4);
kernel = hFit.kernel;
switch(kernel)
  case 'normal'
    kernel = getString(message('stats:dfstrings:kernel_Normal'));
  case 'box'
    kernel = getString(message('stats:dfstrings:kernel_Box'));
  case 'triangle'
    kernel = getString(message('stats:dfstrings:kernel_Triangle'));
  case 'epanechnikov'
    kernel = getString(message('stats:dfstrings:kernel_Epanechnikov'));
  otherwise
    % unexpected, just use the kernel name as stored
end

t{1} = getString(message('stats:dfstrings:Table_Kernel', kernel));
t{2} = getString(message('stats:dfstrings:Table_Bandwidth', sprintf('%g',hFit.bandwidth)));
t{3} = getString(message('stats:dfstrings:Table_Domain', spt));


% ------------------------------------
function t = getparamresults(hFit,spt)
%GETPARAMRESULTS Get results for parametric fit
% Fill cell array with some summary information
t = cell(1,7);
dist = hFit.distspec;
t{1} = getString(message('stats:dfstrings:Table2_Distribution', dist.name));
if isempty(hFit.loglik)
   t{2} = '';
else
   t{2} = getString(message('stats:dfstrings:Table2_LogLikelihood', sprintf('%g',hFit.loglik)));
end
t{3} = getString(message('stats:dfstrings:Table2_Domain_1', spt));
try
if ~isfield(dist,'statfunc') || isempty(dist.statfunc)
        distmean = mean(hFit.probdist);
        distvar = var(hFit.probdist);
else
   pcell = num2cell(hFit.params);
   [distmean,distvar] = feval(dist.statfunc,pcell{:});
    end
   t{4} = sprintf('%s\n%s\n',...
             getString(message('stats:dfstrings:Table2_Mean', sprintf('%g',distmean))),...
             getString(message('stats:dfstrings:Table2_Variance', sprintf('%g',distvar))));
catch
    t{4} = ' ';
end

% Add the parameter names, estimates, and standard errors
if isempty(hFit.pcov)
   t{5} = maketable(getString(message('stats:dfstrings:Table1_Parameter')), ...
                    dist.pnames, ...
                    {getString(message('stats:dfstrings:Table1_Estimate'))}, ...
                    hFit.params(:));
else
   t{5} = maketable(getString(message('stats:dfstrings:Table1_Parameter')), ...
                    dist.pnames, ...
                    {   getString(message('stats:dfstrings:Table1_Estimate')) ...
                        getString(message('stats:dfstrings:Table1_StdErr'))}, ...
                    [hFit.params(:) sqrt(diag(hFit.pcov))]);
end

% Add the covariance matrix, if any
if ~isempty(hFit.pcov)
   t{6} = getString(message('stats:dfstrings:sprintf_EstimatedCovariance'));
   t{7} = maketable(' ', dist.pnames, dist.pnames, hFit.pcov);
end


% -----------------
function t = maketable(rcname, rnames, cnames, data)
%MAKETABLE Make display of a table with row and column names

r = size(data,1);
c = size(data,2);

% Provide empty column names if necessary
if isempty(rcname) && ~isempty(cnames)
   rcname = ' ';
elseif isempty(cnames) && ~isempty(rcname)
   cnames = repmat({' '},1,c);
end

% Blanks to go between columns
blanks = repmat(' ', r+~(isempty(cnames) | isempty(rcname)), 2);

if isempty(cnames)
   % Create body of table without column names
   t = num2str(data(:,1),'%10g');
   for j=2:c
      t = [t, blanks, num2str(data(:,j),'%10g')];
   end
else
   % Create body of table with column names
   t = char(cnames{1}, num2str(data(:,1),'%10g'));
   for j=2:c
      t = [t, blanks, char(cnames{j}, num2str(data(:,j),'%10g'))];
   end
end

if ~isempty(rcname)
   rnames = char(rcname, rnames{:});
   t = [rnames, blanks, t];
end
