function x = finv(p,v1,v2)
%FINV   Inverse of the F cumulative distribution function.
%   X=FINV(P,V1,V2) returns the inverse of the F distribution 
%   function with V1 and V2 degrees of freedom, at the values in P.
%
%   The size of X is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also FCDF, FPDF, FRND, FSTAT, ICDF.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.6.2

%   Copyright 1993-2015 The MathWorks, Inc.


if nargin <  3, 
    error(message('stats:finv:TooFewInputs')); 
end

[errorcode, p, v1, v2] = distchck(3,p,v1,v2);

if errorcode > 0
    error(message('stats:finv:InputSizeMismatch'));
end

% Weed out any out of range parameters or probabilities.
okV = (0 < v1 & v1 < Inf) & (0 < v2 & v2 < Inf);
k = (okV & (0 <= p & p <= 1));
allOK = all(k(:));

% Fill in NaNs for out of range cases.
if ~allOK
    x = nan(size(k),'like',internal.stats.dominantType(p,v1,v2)); % Single if p, v1, or v2 is.

    % Special cases - v1 and/or v2 = Inf
    
    % If v1 = Inf but v2 is not Inf AND both v1, v2 > 0 AND 0 <= p <= 1
    % then finv( theta, v1, v2 ) = v2./chi2inv( 1-theta, v2 )       
    kinf1 = ( v2 > 0 & v2 < Inf ) & (v1 == Inf) & ( 0 <= p & p <= 1 );
    if any(kinf1(:))
        x(kinf1) = v2(kinf1) ./ chi2inv( 1-p(kinf1), v2(kinf1) ); 
    end
    clear kinf1;
    
    % If v1 is not Inf but v2 = Inf AND both v1, v2 > 0 AND 0 <= p <= 1
    % then finv( theta, v1, v2 ) = (1/v1)*chi2inv( theta, v1 )   
    kinf2 = ( v1 > 0 & v1 < Inf ) & (v2 == Inf) & ( 0 <= p & p <= 1 );
    if any(kinf2(:))
        x(kinf2) = chi2inv( p(kinf2), v1(kinf2) ) ./ v1(kinf2);
    end
    clear kinf2;
    
    % If v1 = Inf AND v2 = Inf AND 0 <= p <= 1
    % Define: finv( t, Inf, Inf ) = 0 if t = 0 and 1 if t > 0 and t <= 1 
    kinf3 = (v1 == Inf) & (v2 == Inf) & ( 0 <= p & p <= 1 );    
    if any(kinf3(:))
        % set all kfin3 values to 1
        x(kinf3) = 1;    
        % set values for p = 0 to 0
        x( kinf3 & (p==0) ) = 0;
    end
    clear kinf3;    
    
    % Remove the out of range cases.  If there's nothing remaining, return.
    if any(k(:))
        if numel(p) > 1, p = p(k); end
        if numel(v1) > 1, v1 = v1(k); end
        if numel(v2) > 1, v2 = v2(k); end
    else
        return;
    end
end

% Get the smaller of z or 1-z to give the best precision
up = (p > betainc(.5,v1/2,v2/2,'lower'));
t = zeros(size(p),'like',internal.stats.dominantType(p,v1,v2));
zUp = betaincinv(p(up),v2(up)/2,v1(up)/2,'upper');
t(up) = (1-zUp) ./ zUp;
lo = ~up;
zLo = betaincinv(p(lo),v1(lo)/2,v2(lo)/2,'lower');
t(lo) = zLo ./ (1-zLo);
xk = t .* v2 ./ v1;

% Broadcast the values to the correct place if need be.
if allOK
    x = xk;
else
    x(k) = xk;
end
