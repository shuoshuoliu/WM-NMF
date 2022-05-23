function [a,di] = iwishrnd(tau,df,di)
%IWISHRND Generate inverse Wishart random matrix
%   W=IWISHRND(TAU,DF) generates a random matrix W from the inverse
%   Wishart distribution with parameters TAU and DF.  The inverse of W
%   has the Wishart distribution with covariance matrix inv(TAU) and DF
%   degrees of freedom.
%
%   W=IWISHRND(TAU,DF,DI) expects DI to be lower triangular so that
%   DI'*DI = INV(TAU), i.e., the transpose of inverse of the Cholesky
%   factor of TAU. If you call IWISHRND multiple times using the same
%   value of TAU, it's more efficient to supply DI instead of computing
%   it each time.
%
%   [W,DI]=IWISHRND(TAU,DF) returns DI so it can be used again in
%   future calls to IWISHRND.
%
%   Note that different sources use different parameterizations for the
%   inverse Wishart distribution.  This function defines the parameter
%   TAU so that the mean of the output matrix is TAU/(DF-K-1), where
%   K is the number of rows and columns in TAU.
%
%   See also WISHRND.

%   Copyright 1993-2007 The MathWorks, Inc.


% Error checking
if nargin<2
   error(message('stats:iwishrnd:TooFewInputs'));
end

[n,m] = size(tau);
if n~=m
   error(message('stats:iwishrnd:BadShapedCovariance'));
end

% Factor tau unless that has already been done
if nargin<3
   [d,p] = cholcov(tau,0);
   if p~=0
      error(message('stats:iwishrnd:BadCovariance'));
   end
   if nargout > 1
      di = d' \ eye(size(d));
   end
elseif ~isempty(tau)
   if ~isequal(size(di),size(tau))
      error(message('stats:iwishrnd:BadCovFactor'))
   end
else
   n = size(di,2);
end

if (~isscalar(df)) || (df<=0)
   error(message('stats:iwishrnd:BadDf'));
elseif (df<n) % require this to ensure invertibility
   error(message('stats:iwishrnd:DfTooSmall'));
end

% For small degrees of freedom, generate the matrix using the definition
% of the Wishart distribution; see Krzanowski for example
if (df <= 81+n) && (df==round(df))
   x = randn(df,n);

% Otherwise use the Smith & Hocking procedure
else
   % Load diagonal elements with square root of chi-square variates
   x = diag(sqrt(chi2rnd(df-(0:n-1))));

   % Load upper triangle with independent normal (0, 1) variates
   x(itriu(n)) = randn(n*(n-1)/2,1);
end

% Desired random matrix is INV(DI'*(X'*X)*DI) = D'*INV(X'*X)*D

% Use Cholesky factor for TAU, D ...
if nargin<3
   [~,R] = qr(x,0);
   T = d' / R;
   
% ... or use the Cholesky factor for TAU, DI
else
   [~,R] = qr(x*di,0);
   T = R \ eye(size(R,2));
end

a = T*T';


% --------- get indices of upper triangle of p-by-p matrix
function d = itriu(p)

d = ones(p*(p-1)/2,1);
d(1+cumsum(0:p-2)) = p+1:-1:3;
d = cumsum(d);
