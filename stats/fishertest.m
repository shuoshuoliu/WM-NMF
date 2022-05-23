function [h,p,stats] = fishertest(x,varargin)
%FISHERTEST Fisher's exact test
%   H = FISHERTEST(X) performs Fisher's exact test on a 2-by-2 contingency
%   table X in matrix or table form. This is a test of the hypothesis that
%   there are no non-random associations between the two 2-level
%   categorical variables in X. FISHERTEST returns the result of the test
%   in H. H = 0 indicates that the null hypothesis ('no association')
%   cannot be rejected at the 5% significance level. H = 1 indicates that
%   the null hypothesis can be rejected at the 5% level. X must contain
%   only nonnegative integers. Function CROSSTAB can be used to generate
%   the contingency table from samples of two categorical variables.
%   Fisher's exact test is not suitable when all integers in X are very
%   large. Chi-square test is suggested in this case. 
%
%   [H,P] = FISHERTEST(X) returns the p-value. i.e., the probability of
%   observing the given result, or one more extreme, by chance if the null
%   hypothesis is true. Small values of P cast doubt on the validity of the
%   null hypothesis.
%
%   [H,P,STATS] = FISHERTEST(X) returns a structure with the following fields: 
%      'OddsRatio'           - the odds ratio 
%      'ConfidenceInterval'  - the asymptotic confidence interval for the
%                              odds ratio. If any of the four entries in
%                              the contingency table X is zero, the
%                              confidence interval will not be computed,
%                              and [-Inf Inf] will be displayed.
%
%   [...] = FISHERTEST(X,'PARAM1',VALUE1,'PARAM2',VALUE2) specifies
%   additional parameter name/value pairs chosen from the following:
%      Name       Value 
%      'alpha'  - A value ALPHA between 0 and 1 specifying the significance
%                 level as (100*ALPHA)%. Default is 0.05 for 5% significance.
%      'tail'   - A string specifying the alternative hypothesis:
%          'both'    odds ratio not equal to 1, indicating association
%                    between two variables (two-tailed test, default)
%          'right'   odds ratio greater than 1 (right-tailed test) 
%          'left'    odds ratio is less than 1 (left-tailed test)
%
%
%   Example:
%       x = [3 1;1 3] 
%       [h, p, stats] = fishertest(x);
%
%   See also CROSSTAB, CHI2GOF.

%   Copyright 2014-2015 The MathWorks, Inc.



if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if isa(x,'table');
    x = table2array(x);
elseif isa(x,'dataset')
    x = double(x);
end
if ~ismatrix(x) || ~all(size(x)==2)
    error(message('stats:fishertest:BadInputX'));
end
if any(x(:)<0) || any(isnan(x(:))) || any(isinf(x(:))) || ~internal.stats.isIntegerVals(x)
    error(message('stats:fishertest:NonnegativeInteger'));
end
if all(x(:)>=1e7)
    error(message('stats:fishertest:LargeEntry'));
end

alpha = 0.05;
tail = 'both';
okargs =   {'alpha' 'tail'};
defaults = {alpha   tail};
[alpha, tail] = internal.stats.parseArgs(okargs,defaults,varargin{:});

if  ~isscalar(alpha) || ~isnumeric(alpha) || isnan(alpha) || (alpha <= 0) || (alpha >= 1)
    error(message('stats:fishertest:BadAlpha'));
end
tail = internal.stats.getParamVal(tail,{'both' 'right'  'left'},'''tail''');

r1 = sum(x(1,:));
r2 = sum(x(2,:));
c1 = sum(x(:,1));
c2 = sum(x(:,2));
n = sum(x(:));

try
    if strcmp(tail,'left')       
        p = hygecdf(x(1,1),n,r1,c1);       
    else
        if min(r1,c1)<= min(r2,c2)
            x11 = (0 : min(r1,c1))';
        else
            x22 = (0 : min(r2,c2))';
            x12 = c2 - x22;
            x11 = r1 - x12;
        end       
        switch tail
            case 'both'
                p1 = hygepdf(x(1,1),n,r1,c1);
                p2 = hygepdf(x11,n,r1,c1);
                p = sum(p2(p2 < p1+10*eps(p1)));
            case 'right'
                xr = x11(x11 >= x(1,1));
                p = sum(hygepdf(xr,n,r1,c1));
        end
    end   
catch e
    if any(strcmpi(e.identifier,{'MATLAB:nomem','MATLAB:pmaxsize'}))
        error(message('stats:fishertest:LargeEntry'));
    end
end

h = (p<=alpha);

if nargout > 2
    or = x(1,1)*x(2,2)/x(1,2)/x(2,1);
    if any(x(:)==0)
        CI = [-Inf Inf];
    else
        se = sqrt(1/x(1,1)+1/x(1,2)+1/x(2,1)+1/x(2,2));
        LB = or*exp(-norminv(1-alpha/2)*se);
        UB = or*exp(norminv(1-alpha/2)*se);
        CI = [LB UB];
    end
    stats = struct('OddsRatio', or, 'ConfidenceInterval', CI);
end

end





