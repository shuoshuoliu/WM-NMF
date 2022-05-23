function X = lhsdesign(n,p,varargin)
%LHSDESIGN Generate a latin hypercube sample.
%   X=LHSDESIGN(N,P) generates a latin hypercube sample X containing N
%   values on each of P variables.  For each column, the N values are
%   randomly distributed with one from each interval (0,1/N), (1/N,2/N),
%   ..., (1-1/N,1), and they are randomly permuted.
%
%   X=LHSDESIGN(...,'PARAM1',val1,'PARAM2',val2,...) specifies parameter
%   name/value pairs to control the sample generation.  Valid parameters
%   are the following:
%
%      Parameter    Value
%      'smooth'     'on' (the default) to produce points as above, or
%                   'off' to produces points at the midpoints of
%                   the above intervals:  .5/N, 1.5/N, ..., 1-.5/N.
%      'iterations' The maximum number of iterations to perform in an
%                   attempt to improve the design (default=5)
%      'criterion'  The criterion to use to measure design improvement,
%                   chosen from 'maximin' (the default) to maximize the
%                   minimum distance between points, 'correlation' to
%                   reduce correlation, or 'none' to do no iteration.
%                   Smoothing is turned off with the 'correlation'
%                   criterion.
%
%   Latin hypercube designs are useful when you need a sample that is
%   random but that is guaranteed to be relatively uniformly distributed
%   over each dimension.
%
%   Example:  The following commands show that the output from lhsdesign
%             looks uniformly distributed in two dimensions, but too
%             uniform (non-random) in each single dimension.  Repeat the
%             same commands with x=rand(100,2) to see the difference.
%
%      x = lhsdesign(100,2);
%      subplot(2,2,1); plot(x(:,1), x(:,2), 'o');
%      subplot(2,2,2); hist(x(:,2));
%      subplot(2,2,3); hist(x(:,1));
%
%   See also LHSNORM, UNIFRND.

%   Copyright 1993-2020 The MathWorks, Inc.


% Check input arguments
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

okargs = {'iterations' 'criterion' 'smooth'};
defaults = {NaN 'maximin' 'on'};
[maxiter,crit,dosmooth] = internal.stats.parseArgs(okargs,defaults,varargin{:});

if isempty(maxiter)
   maxiter = NaN;
elseif ~isnumeric(maxiter) || ~isscalar(maxiter) || maxiter<0
   error(message('stats:lhsdesign:ScalarRequired'));
end
if isnan(maxiter), maxiter = 5; end

okcrit = {'none' 'maximin' 'correlation'};
crit = internal.stats.getParamVal(crit,okcrit,'Criterion');

if isempty(dosmooth)
   dosmooth = 'on';
elseif (~isequal(dosmooth,'on')) & (~isequal(dosmooth,'off'))
   error(message('stats:lhsdesign:BadSmooth'));
end

% Start with a plain lhs sample over a grid
X = getsample(n,p,dosmooth);

% Create designs, save best one
if isequal(crit,'none') || size(X,1)<2
    maxiter = 0;
end
switch(crit)
 case 'maximin'
   bestscore = score(X,crit);
   for j=2:maxiter
      x = getsample(n,p,dosmooth);
      
      newscore = score(x,crit);
      if newscore > bestscore
         X = x;
         bestscore = newscore;
      end
   end
 case 'correlation'
   bestscore = score(X,crit);
   for iter=2:maxiter
      % Forward ranked Gram-Schmidt step:
      for j=2:p
         for k=1:j-1
            z = takeout(X(:,j),X(:,k));
            X(:,k) = (rank(z) - 0.5) / n;
         end
      end
      % Backward ranked Gram-Schmidt step:
      for j=p-1:-1:1
         for k=p:-1:j+1
            z = takeout(X(:,j),X(:,k));
            X(:,k) = (rank(z) - 0.5) / n;
         end
      end
   
      % Check for convergence
      newscore = score(X,crit);
      if newscore <= bestscore
         break;
      else
         bestscore = newscore;
      end
   end
end

% ---------------------
function x = getsample(n,p,dosmooth)
x = rand(n,p);
for i=1:p
   x(:,i) = rank(x(:,i));
end
   if isequal(dosmooth,'on')
      x = x - rand(size(x));
   else
      x = x - 0.5;
   end
   x = x / n;
   
% ---------------------
function s = score(x,crit)
% compute score function, larger is better

if size(x,1)<2
    s = 0;       % score is meaningless with just one point
    return
end

switch(crit)
 case 'correlation'
   % Minimize the sum of between-column squared correlations
   c = corrcoef(x);
   s = -sum(sum(triu(c,1).^2));

 case 'maximin'
   % Maximize the minimum point-to-point difference
   [~,dist] = knnsearch(x,x,'k',2);
   s = min(dist(:,2));
 
end

% ------------------------
function z=takeout(x,y)

% Remove from y its projection onto x, ignoring constant terms
xc = x - mean(x);
yc = y - mean(y);
b = (xc-mean(xc))\(yc-mean(yc));
z = y - b*xc;

% -----------------------
function r=rank(x)

% Similar to tiedrank, but no adjustment for ties here
[~, rowidx] = sort(x);
r(rowidx) = 1:length(x);
r = r(:);
