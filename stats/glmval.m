function [yhat,dylo,dyhi] = glmval(beta,x,link,varargin)
%GLMVAL Predict values for a generalized linear model.
%   YHAT = GLMVAL(B,X,LINK) computes predicted values for the generalized
%   linear model with link function LINK and predictor values X.  GLMVAL
%   automatically includes a constant term in the model (do not enter a column
%   of ones directly into X).  B is a vector of coefficient estimates as
%   returned by the GLMFIT function.  LINK can be any of the link function
%   specifications acceptable to the GLMFIT function.
%
%   [YHAT,DYLO,DYHI] = GLMVAL(B,X,LINK,STATS) also computes 95% confidence
%   bounds on the predicted Y values.  STATS is the stats structure
%   returned by GLMFIT.  DYLO and DYHI define a lower confidence bound of
%   YHAT-DYLO and an upper confidence bound of YHAT+DYHI.  Confidence bounds
%   are non-simultaneous and they apply to the fitted curve, not to a new
%   observation.
%
%   [...] = GLMVAL(...,'PARAM1',val1,'PARAM2',val2,...) allows you to
%   specify optional parameter name/value pairs to control the predicted
%   values.  Parameters are:
%
%      'Alpha'         Specifies the confidence level as 100(1-ALPHA)%.
%                      The default is 0.05 for 95% confidence.
%
%      'BinomialSize'  The size parameter (N) for a binomial model.  This
%                      may be a scalar, or a vector with one value for each
%                      row of X. Default is 1.
%
%      'Constant'      Either 'on' (the default) if the model fit included
%                      a constant term, or 'off' if not.  The coefficient of
%                      the constant term should be in the first element of B.
%
%      'Offset'        A vector to use as an additional predictor variable, but
%                      with a coefficient value fixed at 1.0. Default is to
%                      use no offset vector.
%
%      'Simultaneous'  Either true to compute simultaneous confidence
%                      intervals, or false (default) to compute
%                      non-simultaneous intervals.
%
%   Example:  Display the fitted probabilities from a probit regression
%   model for y on x.  Each y(i) is the number of successes in n(i) trials.
%
%       x = [2100 2300 2500 2700 2900 3100 3300 3500 3700 3900 4100 4300]';
%       n = [48 42 31 34 31 21 23 23 21 16 17 21]';
%       y = [1 2 0 3 8 8 14 17 19 15 17 21]';
%       b = glmfit(x, [y n], 'binomial', 'link', 'probit');
%       yfit = glmval(b, x, 'probit', 'binomialsize', n);
%       plot(x, y./n, 'o', x, yfit./n, '-')
%
%   See also GLMFIT.

%   References:
%      [1] Dobson, A.J. (2001) An Introduction to Generalized Linear
%          Models, 2nd edition, Chapman&Hall/CRC Press.
%      [2] McCullagh, P., and J.A. Nelder (1990) Generalized Linear
%          Models, 2nd edition, Chapman&Hall/CRC Press.
%      [3] Collett, D. (2002) Modelling Binary Data, 2nd edition,
%          Chapman&Hall/CRC Press.

%   Copyright 1993-2020 The MathWorks, Inc.


if nargin > 2
    link = convertStringsToChars(link);
end

if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 3
   error(message('stats:glmval:TooFewInputs'));
end

% Get STATS if it's there.
optnargin = nargin - 3;
if optnargin > 0 && ~ischar(varargin{1}) % not a parameter name, assume it's STATS
    stats = varargin{1};
    varargin(1) = [];
    optnargin = optnargin - 1;
else
    stats = [];
end

% Process optional name/value pairs.
if optnargin > 0 && ischar(varargin{1}) % assume it's a parameter name, not CONFLEV
    paramNames = {'confidence' 'size' 'offset' 'constant' 'simultaneous' 'alpha' 'binomialsize'};
    paramDflts = {        .95      1        0        'on' false          []      []};
    [clev,N,offset,const,simul,alpha,binomN,supplied] = ...
             internal.stats.parseArgs(paramNames, paramDflts, varargin{:});
    if supplied.confidence && supplied.alpha
        error(message('stats:glmval:ArgCombination','Confidence','Alpha'));
    end
    if supplied.alpha
        clev = 1-alpha;
    end
    if supplied.size && supplied.binomialsize
        error(message('stats:glmval:ArgCombination','Size','BinomialSize'));
    end
    if supplied.binomialsize
        N = binomN;
    end

else % the old syntax glmval(beta,x,link,stats,clev,N,offset,const)
    clev = .95;
    N = 1;
    offset = 0;
    const = 'on';
    simul = false;
    if optnargin > 0 && ~isempty(varargin{1}), clev = varargin{1}; end
    if optnargin > 1 && ~isempty(varargin{2}), N = varargin{2}; end
    if optnargin > 2 && ~isempty(varargin{3}), offset = varargin{3}; end
    if optnargin > 3 && ~isempty(varargin{4}), const = varargin{4}; end
end
isconst = isequal(const,'on');

% Instantiate functions for one of the canned links, or validate a
% user-defined link specification.
[~,~,ilinkFun] = stattestlink(link,underlyingType(x));

% Should X be changed to a column vector?
if isrow(x)
   if (length(beta)==2 && isconst) || (isscalar(beta) && ~isconst)
      x = x(:);
      if isvector(N), N = N(:); end
      if isvector(offset), offset = offset(:); end
   end
end

if ~isscalar(N) && (~iscolumn(N) || numel(N)~=size(x,1))
    error(message('stats:glmval:WrongArgSize','BinomialSize'));
end
if ~isscalar(offset) && (~iscolumn(offset) || numel(offset)~=size(x,1))
    error(message('stats:glmval:WrongArgSize','Offset'));
end

% Add constant column to X matrix, compute linear combination, and
% use the inverse link function to get a predicted value
if isconst, x = [ones(size(x,1),1,"like",x) x]; end
eta = x*beta + offset;
yhat = N .* ilinkFun(eta);

% Return bounds if requested
if nargout > 1
    if isempty(stats)
        error(message('stats:glmval:NeedSTATS'));
    end
    if ~isnan(stats.s) % dfe > 0 or estdisp == 'off'
        se = stats.se(:);
        cc = stats.coeffcorr;
        V = (se * se') .* cc;
        % R = cholcov(V); <-- this version fails if V is singular
        % vxb = sum((R * x').^2,1);
        vxb = sum((x*V).*x,2)';
        if stats.estdisp
            dfe = stats.dfe;
        else
            dfe = Inf;
        end
        if simul
            dfr = length(beta);
            if dfe==Inf
                crit = sqrt(chi2inv(clev,dfr));
            else
                crit = sqrt(dfr*finv(clev,dfr,dfe));
            end
        else
            crit = tinv((1+clev)/2, dfe);
        end

        dxb = crit * sqrt(vxb(:));
        dyhilo = [N.*ilinkFun(eta-dxb) N.*ilinkFun(eta+dxb)];
        dylo = yhat - min(dyhilo,[],2);
        dyhi = max(dyhilo,[],2) - yhat;
    else
        dylo = NaN(size(yhat),"like",yhat);
        dyhi = NaN(size(yhat),"like",yhat);
    end
end
