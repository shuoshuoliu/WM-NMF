function [p, h, stats] = ranksum(x,y,varargin)
%RANKSUM Wilcoxon rank sum test for equal medians.
%   P = RANKSUM(X,Y) performs a two-sided rank sum test of the hypothesis
%   that two independent samples, in the vectors X and Y, come from
%   distributions with equal medians, and returns the p-value from the
%   test.  P is the probability of observing the given result, or one more
%   extreme, by chance if the null hypothesis ("medians are equal") is
%   true.  Small values of P cast doubt on the validity of the null
%   hypothesis.  The two sets of data are assumed to come from continuous
%   distributions that are identical except possibly for a location shift,
%   but are otherwise arbitrary.  X and Y can be different lengths.
%   RANKSUM treats NaNs in X or Y as missing values, and removes them.
%   The two-sided p-value is computed by doubling the most significant
%   one-sided value.
%
%   The Wilcoxon rank sum test is equivalent to the Mann-Whitney U test.
%
%   [P,H] = RANKSUM(...) returns the result of the hypothesis test,
%   performed at the 0.05 significance level, in H.  H=0 indicates that
%   the null hypothesis ("medians are equal") cannot be rejected at the 5%
%   level. H=1 indicates that the null hypothesis can be rejected at the
%   5% level.
%
%   [P,H] = RANKSUM(...,'alpha',ALPHA) returns the result of the hypothesis
%   test performed at the significance level ALPHA.
%
%   [P,H] = RANKSUM(...,'method',M) computes the p-value exactly if M is
%   'exact', or uses a normal approximation if M is 'approximate'.  If you
%   omit this argument, RANKSUM uses the exact method for small samples and
%   the approximate method for larger samples.
%
%   [P,H] = RANKSUM(...,'tail',TAIL) performs the test against the
%   alternative hypothesis specified by TAIL:
%       'both'  -- "medians are not equal" (two-tailed test, default)
%       'right' -- "median of X is greater than median of Y" (right-tailed test)
%       'left'  -- "median of X is less than median of Y" (left-tailed test)
%   TAIL must be a single string.
%
%   [P,H,STATS] = RANKSUM(...) returns STATS, a structure with one or two
%   fields.  The field 'ranksum' contains the value of the rank sum
%   statistic for X.  For the 'approximate' method, the field 'zval'
%   contains the value of the normal (Z) statistic.
%
%   See also SIGNTEST, SIGNRANK, KRUSKALWALLIS, TTEST2.

%   References:
%      [1] Hollander, M. and D. A. Wolfe.  Nonparametric Statistical
%          Methods. Wiley, 1973.
%      [2] Gibbons, J.D.  Nonparametric Statistical Inference,
%          2nd ed.  M. Dekker, 1985.


% Check that x and y are vectors
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if ~isvector(x) || ~isvector(y)
   error(message('stats:ranksum:InvalidData'));
end

% Remove missing data
x = x( ~isnan(x) );
y = y( ~isnan(y) );
if isempty(x) || isempty(y)
	error(message('stats:signrank:NotEnoughData'));
end

% Determine value for 'alpha' and parse inputs
alpha = 0.05;   % default
if nargin>2 && isnumeric(varargin{1})
   % Grandfathered syntax:  ranksum(x,y,alpha)
   alpha = varargin{1};
   varargin(1) = [];
end
%
oknames = {'alpha' 'method' 'tail'};
dflts   = {alpha   ''   'both'};
[alpha,method,tail] = internal.stats.parseArgs(oknames,dflts,varargin{:});

% Check value of 'alpha'
if ~isscalar(alpha) || ~isnumeric(alpha) || isnan(alpha) || (alpha <= 0) || (alpha >= 1)
   error(message('stats:ranksum:BadAlpha'));
end

% Check value of 'tail'
tail = internal.stats.getParamVal(tail, {'both'  'right'  'left'}, '''tail''');

% Determine and check value for 'method'
nx = numel(x);
ny = numel(y);
ns = min(nx,ny);
if isempty(method)
   if (ns < 10)  &&  ((nx+ny) < 20)
      method = 'exact';
   else
      method = 'approximate';
   end
elseif strcmpi(method, 'oldexact')
	method = 'oldexact';
else   % method not recognized, throw error
   method = internal.stats.getParamVal(method,{'exact' 'approximate'},'''method''');
end

% Determine computational 'technique'
switch method
	case 'approximate'
		technique = 'normal_approximation';
	case 'oldexact'
		technique = 'full_enumeration';
	case 'exact'
		if (nx+ny) < 10
			technique = 'full_enumeration';
		else
			technique = 'network_algorithm';
		end
end

%      %      %      %      %      %      %      %      %      %
% Calculations for Rank Sum Test

x = x(:);   % ensure columns
y = y(:);
if nx <= ny
   smsample = x;
   lgsample = y;
   same_order = true;
else
   smsample = y;
   lgsample = x;
   same_order = false;
end

% Compute the rank sum statistic based on the smaller sample
[ranks, tieadj] = tiedrank([smsample; lgsample]);
srank = ranks(1:ns);
w = sum(srank);


switch technique
	case 'full_enumeration'
		allpos = nchoosek(ranks,ns);   % enumerate all possibilities
		sumranks = sum(allpos,2);
		np = size(sumranks, 1);
		
		switch tail
			case 'both'
				plo = sum( sumranks <= w) / np ;
				phi = sum( sumranks >= w) / np ;
				p_tail = min(plo,phi);
				p = min(2*p_tail, 1);   % 2-sided, p>1 means middle is double-counted
				
			case 'right'
				switch same_order
					case true
						p = sum( sumranks >= w) / np ;
					case false
						p = sum( sumranks <= w) / np;
				end
				
			case 'left'
				switch same_order
					case true
						p = sum( sumranks <= w) / np ;
					case false
						p = sum( sumranks >= w) / np;
				end
				
		end
		
		%     %     %     %     %     %     %      %      %      %
		
	case 'network_algorithm'
		[p_net, pvals] = exactprob(smsample, lgsample, w);
		
		if any(isnan(p_net)) || any(isnan(pvals))
			warning(message('stats:ranksum:NanResult'));
			p = NaN;
			
		else
			switch tail
				case 'both'   % two-tailed test
					p = min(2*p_net, 1);   % p>1 means the middle is double-counted
					
				case 'right'   % right-tail test
					switch same_order
						case true
							p =  pvals(2) + pvals(3);
						case false
							p = pvals(2) + pvals(1);
					end
					
				case 'left'   % left-tail test
					switch same_order
						case true
							p =  pvals(2) + pvals(1);
						case false
							p = pvals(2) + pvals(3);
					end
					
			end   % conditional on 'tail'
		
		end
		
		%     %     %     %     %     %     %      %      %      %
				
	case 'normal_approximation'
		wmean = ns*(nx + ny + 1)/2;
		tiescor = 2 * tieadj / ((nx+ny) * (nx+ny-1));
		wvar  = nx*ny*((nx + ny + 1) - tiescor)/12;
		wc = w - wmean;

		% compute z-value, including continuity correction
		switch tail
			case 'both'
				z = (wc - 0.5 * sign(wc))/sqrt(wvar);
				if ~same_order
					z = -z;
				end
				p = 2*normcdf(-abs(z));
				
			case 'right'
				if same_order
					z = (wc - 0.5)/sqrt(wvar);
				else
					z = -(wc + 0.5)/sqrt(wvar);
				end
				
				p = normcdf(-z);
				
			case 'left'
				if same_order
					z = (wc + 0.5)/sqrt(wvar);
				else
					z = -(wc - 0.5)/sqrt(wvar);
				end
				
				p = normcdf(z);
		end


		if (nargout > 2)   % handle additional outputs
			stats.zval = z;
		end
		
end   % conditional on 'technique'



% Handle additional outputs
if nargout > 1,
   h = (p<=alpha);

   if (nargout > 2)
	   if same_order
		   stats.ranksum = w;
	   else
		   stats.ranksum = sum( ranks(ns+1:end) );
	   end
   end
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [p1, pvals] = exactprob(x,y,w)
%EXACTPROB Exact P-values for Wilcoxon Mann Whitney nonparametric test
%   [P1,PVALS]=EXACTPROB(X,Y,W) computes the p-value P for the test
%   statistic W in a Wilcoxon-Mann-Whitney nonparametric test of the
%   hypothesis that X and Y come from distributions with equal medians.

% Create a contingency table with element (i,j) indicating how many
% times u(j) appears in sample i
u = unique([x(:); y(:)]);
t = zeros(2,length(u));
t(1,:) = histc(x,u)';
t(2,:) = histc(y,u)';

% Compute weights for wmw test
colsum = sum(t,1);
tmp = cumsum(colsum);
wts = [0 tmp(1:end-1)] + .5*(1+diff([0 tmp]));

% Compute p-value using network algorithm for contingency tables
[p1, pvals] = statctexact(t,wts,w);