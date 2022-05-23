function [p,dw]=dwtest(r,X,varargin)
%DWTEST Durbin-Watson test for autocorrelation in linear regression.
%    [P,DW] = DWTEST(R,X) performs a Durbin-Watson test on the vector R of
%    residuals from a linear regression, where X is the design matrix from
%    that linear regression.  P is the computed p-value for the test, and
%    DW is the Durbin-Watson statistic.  The Durbin-Watson test is used to
%    test if the residuals are uncorrelated, against the alternative that
%    there is autocorrelation among them.
%
%   [...] = DWTEST(R,X,'PARAM1',val1,'PARAM2',val2,...) specifies one or
%   more of the following name/value pairs:
%
%       Parameter       Value
%       'method'        'exact' to calculate an exact p-value using the PAN
%                       algorithm (default if the sample size is less than
%                       400), or 'approximate' to calculate the p-value
%                       using a normal approximation (default if the sample
%                       size is 400 or larger).
%       'tail'          A string specifying the alternative hypothesis:
%          'both'   "serial correlation is not 0" (two-tailed test, default)
%          'right'  "serial correlation is greater than 0" (right-tailed test)
%          'left'   "serial correlation is less than 0" (left-tailed test)
%
%   Example:
%      % Fit a straight line to the census data and note the 
%      % autocorrelation in the residuals
%      load census
%      n = length(cdate);
%      X = [ones(n,1), cdate];
%      [b,bint,r1] = regress(pop,X);
%      p1 = dwtest(r1,X)
%      plot(cdate,r1,'b-', cdate,zeros(n,1),'k:')
% 
%      % Adding a squared term reduces the autocorrelation but it is still
%      % significantly different from zero
%      X = [ones(n,1), cdate, cdate.^2];
%      [b,bint,r2] = regress(pop,X);
%      p2 = dwtest(r2,X)
%      line(cdate,r2,'color','r')
%     
%   See also REGRESS.

%   Reference:
%   J. Durbin & G.S. Watson (1950), Testing for Serial Correlation in Least
%   Squares Regression I. Biometrika (37), 409-428.
%
%   R.W. Farebrother (1980), Pan's Procedure for the Tail Probabilities of
%   the Durbin-Watson Statistic. Applied Statistics, 29, 224-227.
%
%   Copyright 1993-2012 The MathWorks, Inc. 


if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

[n,p] = size(X);   
option = '';
alternative= '';

if nargin>=3
    if ismember(lower(varargin{1}),{'exact','approximate'})||isempty(varargin{1})
        % Old syntax
        option = varargin{1};
        if nargin>=4
            alternative = varargin{2};
        end
        
    elseif nargin==3
        error(message('stats:pvaluedw:BadMethod'));
        
    else
        okargs =   {'method' 'tail'};
        defaults = { ''       ''};
        [option, alternative] = internal.stats.parseArgs(okargs,defaults,varargin{:});
    end
end

if isempty(option)
    if n<400
        option = 'exact';
    else
        option = 'approximate';
    end;
else
    option = internal.stats.getParamVal(option,{'exact','approximate'},'''method''');
end

if isempty(alternative)
    alternative = 'both';
else
    alternative = internal.stats.getParamVal(alternative,{'left','both','right'},'''tail''');
end

dw = sum(diff(r).^2)/sum(r.^2); % durbin-watson statistic

% This calls the function of Pan algorithm/normal approximation
% Recall that the distribution of dw depends on the design matrix
% in the regression.
pdw = pvaluedw(dw,X,option);

% p-value depends on the alternative hypothesis.
switch lower(alternative)
    case 'both'  
        p = 2*min(pdw, 1-pdw); 
    case 'left'
        p = 1-pdw;
    case 'right'
        p = pdw;
end


        
        
