function p = betacdf(x,a,b,varargin)
%BETACDF Beta cumulative distribution function.
%   P = BETACDF(X,A,B) returns the beta cumulative distribution
%   function with parameters A and B at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   BETAINC does the computational work.
%
%   P = BETACDF(X,A,B,'upper') returns the upper tail probability of the beta 
%   distribution function with parameters A and B at the values in X.
%
%   See also BETAFIT, BETAINV, BETALIKE, BETAPDF, BETARND, BETASTAT, CDF,
%            BETAINC.
   
%   Reference:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.5.

%   Copyright 1993-2017 The MathWorks, Inc. 


if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin<3, 
   error(message('stats:betacdf:TooFewInputs')); 
end

[errorcode, x, a, b] = distchck(3,x,a,b);

if errorcode > 0
   error(message('stats:betacdf:InputSizeMismatch'));
end

% Weed out any out of range parameters or data.
okAB = (0 < a & a < Inf) & (0 < b & b < Inf);
k = (okAB & (0 <= x & x <= 1));
allOK = all(k(:));

% Fill in NaNs for out of range cases, fill in edges cases when X is outside 0 or 1.
if ~allOK
    p = NaN(size(k),'like',internal.stats.dominantType(x,a,b)); % Single if x, a, or b is
    if nargin>3 && strcmpi(varargin{end},'upper')
        p(okAB & x <= 0) = 1;
        p(okAB & x >= 1) = 0;
    else
        p(okAB & x < 0) = 0;
        p(okAB & x > 1) = 1;
    end
    
    % Remove the bad/edge cases, leaving the easy cases.  If there's
    % nothing remaining, return.
    if any(k(:))
        if numel(x) > 1, x = x(k); end
        if numel(a) > 1, a = a(k); end
        if numel(b) > 1, b = b(k); end
    else
        return;
    end
end

if nargin>3
    if ~strcmpi(varargin{end},'upper')
        error(message('stats:cdf:UpperTailProblem'));
    else
        varargin{end}=lower(varargin{end});
    end
end
pk = betainc(x,a,b,varargin{:});

% Broadcast the values to the correct place if need be.
if allOK
    p = pk;
else
    p(k) = pk;
end
