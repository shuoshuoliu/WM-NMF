function [varargout] = raylcdf(x,varargin)
%RAYLCDF  Rayleigh cumulative distribution function.
%   P = RAYLCDF(X,B) returns the Rayleigh cumulative distribution 
%   function with parameter B at the values in X.
%
%   The size of P is the common size of X and B. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   P = RAYLCDF(X,B,'upper') returns the upper tail probability of 
%   the Rayleigh distribution with parameter B at the values in X.
%
%   See also RAYLFIT, RAYLINV, RAYLPDF, RAYLRND, RAYLSTAT, CDF.

%   Reference:
%      [1]  Evans, Merran, Hastings, Nicholas and Peacock, Brian,
%      "Statistical Distributions, Second Edition", Wiley
%      1993 p. 134-136.

%   Copyright 1993-2017 The MathWorks, Inc.


if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin < 1
    error(message('stats:raylcdf:TooFewInputs')); 
end
if nargin>1 && isequal(varargin{end},'upper')
    %Compute upper tail
    uflag=true;
    varargin(end) = [];
elseif nargin>1 && ischar(varargin{end})&& ~strcmpi(varargin{end},'upper')
    error(message('stats:cdf:UpperTailProblem'));
else
    uflag=false;
end
[varargout{1:max(1,nargout)}] = localraylcdf(x,uflag,varargin{:});

function p = localraylcdf(x,uflag,b)
if nargin<3
    b = 1;
end

[errorcode,x,b] = distchck(2,x,b);

if errorcode > 0
    error(message('stats:raylcdf:InputSizeMismatch'));
end

% Initialize P to zero.
p = zeros(size(x),'like',internal.stats.dominantType(x,b)); % single if x or b is

k0 = b > 0 & x <= 0;
if uflag == true && any(k0(:))
    p(k0)=1;
end

k = b > 0 & x > 0;
if any(k(:))
    xk = x(k);
    bk = b(k);
    if uflag == true
        p(k) = exp(-xk .^ 2 ./ (2*bk .^ 2));
    else
        p(k) = - expm1(-xk .^ 2 ./ (2*bk .^ 2));
    end
end

p(~(k0 | k)) = NaN;
