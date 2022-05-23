function hout = addedvarplot(varargin)
%ADDEDVARPLOT Create added-variable plot for stepwise regression
%   ADDEDVARPLOT(X,Y,VNUM,INMODEL) produces an added-variable plot for the
%   response Y and the predictor in column VNUM of X.  This plot illustrates
%   the incremental effect of this predictor in a regression model in which
%   the columns listed in INMODEL are used as predictors.  X is an N-by-P
%   matrix of predictor values.  Y is vector of N response values.  VNUM is a
%   scalar index specifying the column of X to use in the plot.  INMODEL is
%   either a vector of column numbers or a logical vector of P elements,
%   specifying the columns of X to use in the base model.  By default, all
%   elements of INMODEL are false (the model has no predictors).  ADDEDVARPLOT
%   automatically includes a constant term in the model.
%
%   ADDEDVARPLOT(X,Y,VNUM,INMODEL,STATS) uses the structure STATS containing
%   fitted model results created by the STEPWISEFIT function.  If STATS is
%   omitted, this function computes it.
%
%   ADDEDVARPLOT(AX,...) plots into the axes with handle AX.
%
%   An added-variable plot contains data and fitted lines.  Suppose X1 is
%   column VNUM of X.  The data curve plots Y versus X1 after removing the
%   effects of the other predictors specified by INMODEL.  The solid line is
%   a least squares fit to the data curve, and its slope is the coefficient
%   that X1 would have if it were included in the model.  The dotted lines
%   are 95% confidence bounds for the fitted line, and they can be used to
%   judge the significance of X1.
%
%   If VNUM also appears in INMODEL, the plot that ADDEDVARPLOT produces is
%   sometimes known as a partial regression leverage plot.
%
%   Example:  Perform a stepwise regression on the Hald data, and create
%             an added-variable plot for the predictor in column 2.
%      load hald
%      [b,se,p,in,stats] = stepwisefit(ingredients,heat);
%      addedvarplot(ingredients,heat,2,in,stats)

%   LinearModel.plotAdded may call this function with a vector of VNUM
%   values. It will never call with an "out" term. It supplies the value
%   DOCONST=true because it supplies a constant term in X.

%   Copyright 1993-2020 The MathWorks, Inc.

[ax, varargin, nargin] = axescheck(varargin{:});
narginchk(3,Inf);

% If axes is found then add it as Name-Value pair.
if ~isempty(ax)
    varargin = [varargin, 'Parent', {ax}];
end

% Extract input arguments.
x = varargin{1};
y = varargin{2};
vnum = varargin{3};
varargin = varargin(4:end);


if(nargin > 3)
    in = varargin{1};
    varargin = varargin(2:end);
end
if(nargin > 4)
    stats = varargin{1};
    varargin = varargin(2:end);
end
if(nargin > 5)
    f = varargin{1};
    varargin = varargin(2:end);
end
if(nargin > 6)
    doconst = varargin{1};
    varargin = varargin(2:end);
end

if nargin > 7
    [varargin{:}] = convertStringsToChars(varargin{:});
end



P = size(x, 2);

% Check for valid inputs
if ~isvector(y)
   error(message('stats:addedvarplot:VectorRequired'));
end

if (nargin < 4)
    in = false(1,P);
elseif islogical(in)
   if length(in)~=P
      error(message('stats:addedvarplot:InModelBadSize'));
   end
else
   if any(~ismember(in,1:P))
      error(message('stats:addedvarplot:InModelBadValue'));
   end
   in = ismember((1:P),in);
end

if isempty(vnum) || ~isvector(vnum)
    error(message('stats:addedvarplot:EmptyVarSelection'));
elseif ~(all(in(vnum)) || ~any(in(vnum)))
    error(message('stats:addedvarplot:BadVarSelection'));
end

if nargin<7
    doconst = true;
end

% Perform fit if fit results are not done; otherwise retrieve some results
newFit = nargin<5;
if newFit % never true when called from plotAdded
    [B,~,~,in,stats] = stepwisefit(x,y, 'maxiter', 0, 'display', 'off',...
        'inmodel',in);
    if any(stats.wasnan)
        x = x(~stats.wasnan,:);
        y = y(~stats.wasnan,:);
    end
else
    if ~isstruct(stats) || ~isfield(stats,'source') || ~isequal(stats.source, 'stepwisefit')
      error(message('stats:addedvarplot:BadStats'));
    end
    B = stats.B;
end
if isfield(stats,'mse')
    mse = stats.mse;
else
    mse = [];
end
if isfield(stats,'wts')
    wts = stats.wts;
else
    wts = ones(size(y));
end

% Argument 6 processed below

N = length(y);
sumw = sum(wts);
ymean = (wts'*y)/sumw;
bk = B(vnum);
alpha = 0.05;

if all(in(vnum))
   % Create a partial regression leverage plot for a column that is in.
   % First adjust y for the remaining predictors.
   r = y - x(:,in)*B(in);
   if doconst
       r = r - (wts'*r)/sumw;
   end
   
   % Adjust the chosen x column or columns
   xnotk = [ones(N,doconst) x(:,in & (~ismember(1:P,vnum)))];
   xk = x(:,vnum);
   bxnotk = qrfit(xnotk,xk);
   xr = xk - xnotk*bxnotk;
   
   if isscalar(vnum)
       % Simple for a single term: adjusted X effect plus residual
       yr = xr*B(vnum) + r;
       xmean = (wts'*x(:,vnum))/sumw;
       ttl = getString(message('stats:addedvarplot:PartialRegrXTitle',vnum));
   else
       % For multiple terms, find the linear combination of the terms that
       % best fits y, and convert to a unit direction vector; then subtract
       % the effect of the other terms
       xr = xr*B(vnum);
       yr = xr + r;
       xr = xr/norm(B(vnum));
       xmean = (wts'*(xk*(B(vnum)/norm(B(vnum)))))/sumw;
       ttl = getString(message('stats:addedvarplot:PartialRegrTitle'));
   end
else
   % Created added-variable plot for an X column that is now out
   if stats.dfe==0
      error(message('stats:addedvarplot:NotEnoughData'))
   end
   varlist = find(~in);
   outnum = varlist==vnum;
   xr = stats.xr(:,outnum);
   yr = stats.yr;
   ttl = getString(message('stats:addedvarplot:AddedVarTitle',vnum));
   xmean = (wts'*x(:,vnum))/sumw;
end


% Create informative title
in2 = in;
in2(vnum) = 0;
runstart = find([in2(1), in2(2:end)&~in2(1:end-1)]);
runend   = find([in2(1:end-1)&~in2(2:end), in2(end)]);
txt = '';
for j=1:length(runstart)
   if runstart(j)==runend(j)
      txt = sprintf('%s,X%d',txt,runstart(j));
   elseif runstart(j)==runend(j)-1
      txt = sprintf('%s,X%d,X%d',txt,runstart(j),runend(j));
   else
      txt = sprintf('%s,X%d-X%d',txt,runstart(j),runend(j));
   end
end
if ~isempty(txt)     % add to title after removing extra comma
   ttl = {ttl; getString(message('stats:addedvarplot:AdjustedFor',txt(2:end)))};
end

% Create plot
if nargin>=6 && ~isempty(f) && isempty(ax)
    ax = get(f,'CurrentAxes');
    if isempty(ax)
        ax = axes('Parent', f);
        set(ax,'Position',[.13 .11 .78 0.78]);
    end
    varargin = [{'Parent';ax};varargin(:)];
end

if ~isscalar(bk)
    bknorm = norm(bk);
    bk = bknorm;
end

% Get data points for the plot
xplot = xmean + xr;
yplot = ymean + yr;

% Fit to points, but use dfe from original fit
sw = sqrt(wts);
X = [sw, sw.*xplot];
b = X\(sw.*yplot);
n = length(yplot);
yplotfit = b(1) + b(2)*xplot;
if isempty(mse)
    mse = (wts'*(yplot-yplotfit).^2)/sumw * (n/(n-2));
end
[~,R] = qr(X,0);

% Compute confidence intervals for line at a grid of x values
xconf = linspace(min(xplot),max(xplot))';
E = [ones(size(xconf)),xconf]/R;
dy = -tinv(alpha/2,stats.dfe) * sqrt(sum(E.^2,2) * mse);

% Compute fitted line and bounds on this grid
yfit = b(1) + b(2)*xconf;
upper = yfit+dy;
lower = yfit-dy;

% Restore NaNs to get back to original row numbers
if any(stats.wasnan)
    [xplot,yplot] = statinsertnan(stats.wasnan,xplot,yplot);
end

h = plot(xplot,yplot,'x',varargin{:});
ax = ancestor(h,'axes');

washold = ishold(ax);
if ~washold
    hold(ax,'on');
end
h = [h; plot(ax,xconf,yfit,'r-',...
                [xconf;NaN;xconf],[upper;NaN;lower],'r:')];
if ~washold
    hold(ax,'off');
end

set(h(1),'Tag','data');
set(h(2),'Tag','fit');

title(ttl,'Parent', ax);
xlabel(getString(message('stats:addedvarplot:AdjustedX',vnum(1))), 'Parent', ax);
ylabel(getString(message('stats:addedvarplot:AdjustedY')), 'Parent', ax);

legend(ax,getString(message('stats:addedvarplot:AdjustedData')),...
       getString(message('stats:addedvarplot:FitEquation',sprintf('y=%g*x',bk))),...
       getString(message('stats:addedvarplot:ConfBounds',sprintf('%g',100*(1-alpha)))),...
       'Location','Best');

set(ax,'UserData',{'addedvarplot'}); % flag for the gname function

if nargout>0
    hout = h;
end

function b = qrfit(X,y)
[n,ncolX] = size(X);
[Q,R,perm] = qr(X,0);
if ~isempty(R)
    p = sum(abs(diag(R)) > max(n,ncolX)*eps(R(1)));
    if p < ncolX
        R = R(1:p,1:p);
        Q = Q(:,1:p);
        perm = perm(1:p);
    end
end

% Compute the LS coefficients, filling in zeros in elements corresponding
% to rows of X that were thrown out.
b = zeros(ncolX,size(y,2),"like",internal.stats.dominantType(X,y));
b(perm,:) = R \ (Q'*y);
