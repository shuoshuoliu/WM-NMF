function [varargout] = unifcdf(x,varargin)
%UNIFCDF Uniform (continuous) cumulative distribution function (cdf).
%   P = UNIFCDF(X,A,B) returns the cdf for the uniform distribution
%   on the interval [A,B] at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   By default, A = 0 and B = 1.
%
%   P = UNIFCDF(X,A,B,'upper') returns the upper tail probability
%   for the uniform distribution on the interval [A,B] at the values in X.
%
%   See also UNIFINV, UNIFIT, UNIFPDF, UNIFRND, UNIFSTAT, CDF.

%   Copyright 1993-2015 The MathWorks, Inc.


if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 1,
    error(message('stats:unifcdf:TooFewInputs'));
end

if nargin>1 && strcmpi(varargin{end},'upper')
    %Compute upper tail and remove 'upper' flag
    uflag=true;
    varargin(end) = [];
elseif nargin>1 && ischar(varargin{end})&& ~strcmpi(varargin{end},'upper')
    error(message('stats:cdf:UpperTailProblem'));
else
    uflag=false;
end
[varargout{1:max(1,nargout)}] = localunifcdf(x,uflag,varargin{:});
end

function p = localunifcdf(x,uflag,a,b)
if nargin <3
    a = 0;
    b = 1;
end

[errorcode,x,a,b] = distchck(3,x,a,b);

if errorcode > 0
    error(message('stats:unifcdf:InputSizeMismatch'));
end

% Initialize P to zero.
p = zeros(size(x),'like',internal.stats.dominantType(x,a,b));

k = find(x > a & x < b & a < b);
if  uflag == true
    %Compute upper tail
    p(x <= a & a < b) = 1;
    p(x >= b & a < b) = 0;
    if any(k)
        p(k) = (b(k)- x(k)) ./ (b(k) - a(k));
    end
else
    % (1) x <= a and a < b.
    p(x <= a & a < b) = 0;
    
    % (2) x >= b and a < b.
    p(x >= b & a < b) = 1;
    
    % (3) a < x < b.
    if any(k)
        p(k) = (x(k) - a(k)) ./ (b(k) - a(k));
    end
end

% (4) a >= b then set p to NaN.
p(a >= b) = NaN;

% (5) If x or a or b is NaN, set p to NaN.
p(isnan(x) | isnan(a) | isnan(b)) = NaN;

end
