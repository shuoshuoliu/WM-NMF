function p = nctcdf(x,nu,delta,uflag)
%NCTCDF Noncentral T cumulative distribution function (cdf).
%   P = NCTCDF(X,NU,DELTA) Returns the noncentral T cdf with NU
%   degrees of freedom and noncentrality parameter, DELTA, at the values
%   in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   P = NCTCDF(X,NU,DELTA,'upper') Returns the upper tail probability of
%   the noncentral T distribution with NU degrees of freedom and noncentrality
%   parameter, DELTA, at the values in X.
%
%   See also NCTINV, NCTPDF, NCTRND, NCTSTAT, TCDF, CDF.

%   References:
%      [1]  Johnson, Norman, and Kotz, Samuel, "Distributions in
%      Statistics: Continuous Univariate Distributions-2", Wiley
%      1970 p. 205.
%      [2]  Evans, Merran, Hastings, Nicholas and Peacock, Brian,
%      "Statistical Distributions, Second Edition", Wiley
%      1993 pp. 147-148.

%   Copyright 1993-2016 The MathWorks, Inc.


if nargin > 3
    uflag = convertStringsToChars(uflag);
end

if nargin <  3
    error(message('stats:nctcdf:TooFewInputs'));
end

[errorcode, x, nu, delta] = distchck(3,x,nu,delta);

if errorcode > 0
    error(message('stats:nctcdf:InputSizeMismatch'));
end

if nargin > 3
    %Compute upper tail
    if ~strcmpi(uflag,'upper')
        error(message('stats:cdf:UpperTailProblem'));
    else
        uppertail = true;
    end
else
    uppertail = false;
end

% Initialize p to zero.
if isa(x,'single') || isa(nu,'single') || isa(delta,'single')
    p = zeros(size(x),'single');
    seps = eps('single');
else
    p = zeros(size(x));
    seps = eps;
end

% Find NaN's
hasnan = isnan(x) | isnan(nu) | isnan(delta);
p(hasnan) = NaN;

% Special cases for delta==0 and x<0; and x = Inf.
f1 = (nu <= 0 | isinf(delta)) & ~hasnan;
f0 = delta == 0 & ~f1 & ~hasnan;
fn = x < 0 & ~f0 & ~f1 & ~hasnan;
fInf = x == Inf & ~f0 & ~f1 & ~hasnan;
bignu = nu>2e6 & ~f0 & ~f1 & ~fInf & ~hasnan; % use normal approximation, Johnson & Kotz eq 26.7.10

flag1 = any(f1(:));
flag0 = any(f0(:));
flagn = any(fn(:));
flagInf = any(fInf(:));
flagbignu = any(bignu(:));

if (flag1 || flag0 || flagn || flagInf || flagbignu)
    fp = ~(f1 | f0 | fn | fInf | bignu);
    if flag1,        p(f1) = NaN; end
    if flag0
        if uppertail==true
            p(f0) = tcdf(x(f0),nu(f0),'upper');
        else
            p(f0) = tcdf(x(f0),nu(f0));
        end
    end
    if flagInf
        if uppertail==true
            p(fInf)=0;
        else
            p(fInf) = 1;
        end
    end
    if flagbignu
        s = 1 - 1./(4*nu);
        d = sqrt(1+x.^2./(2*nu));
        if uppertail==true
            p(bignu) = normcdf(x(bignu).*s(bignu),delta(bignu),d(bignu),'upper');
        else
            p(bignu) = normcdf(x(bignu).*s(bignu),delta(bignu),d(bignu));
        end
    end
    if any(fp(:))
        if uppertail==true
            p(fp) = nctcdf(x(fp), nu(fp), delta(fp),'upper');
        else
            p(fp) = nctcdf(x(fp), nu(fp), delta(fp));
        end
    end
    if flagn
        if uppertail==true
            p(fn) = nctcdf(-x(fn), nu(fn), -delta(fn));
        else
            p(fn) = nctcdf(-x(fn), nu(fn), -delta(fn),'upper');
        end
    end
    return
end

%Value passed to Incomplete Beta function.
xsq = x.^2;
denom = nu + xsq;
P = xsq ./ denom;
Q = nu ./ denom;   % Q = 1-P but avoid roundoff if P is close to 1

% Set up for infinite sum.
dsq = delta.^2;

% Compute probability P[t<0] + P[0<t<x], starting with 1st term
if uppertail==true
    fx0 = x == 0 & ~hasnan;
    if any(fx0(:))
        fx = normcdf(-delta,0,1,'upper');
        p(fx0)= fx(fx0);
    end
else
    p(~hasnan) = normcdf(-delta(~hasnan),0,1);
end

% Now sum a series to compute the second term
k0 = find(x~=0 & ~hasnan);
if any(k0(:))
    P = P(k0);
    Q = Q(k0);
    nu = nu(k0);
    dsq = dsq(k0);
    signd = sign(delta(k0));
    subtotal = zeros(size(k0));
    
    % Start looping over term jj and higher, this should be near the
    % peak of the E part of the term (see below)
    jj = 2 * floor(dsq/2);
    
    % Compute an infinite sum using Johnson & Kotz eq 9, or new
    % edition eq 31.16, each term having this form:
    %      B  = betainc(P,(j+1)/2,nu/2);
    %      E  = (exp(0.5*j*log(0.5*delta^2) - gammaln(j/2+1)));
    %      term = E .* B;
    %
    % We'll compute betainc at the beginning, and then update using
    % recurrence formulas (Abramowitz & Stegun 26.5.16).  We'll sum the
    % series two terms at a time to make the recurrence work out.
    
    E1 =          exp(0.5* jj   .*log(0.5*dsq) - dsq/2 - gammaln( jj   /2+1));
    E2 = signd .* exp(0.5*(jj+1).*log(0.5*dsq) - dsq/2 - gammaln((jj+1)/2+1));
    
    % Use either P or Q, whichever is more accurately computed
    t = (P < 0.5);   % or maybe < dsq./(dsq+nu)
    B1 = zeros(size(P));
    B2 = zeros(size(P));
    if uppertail==true
        if any(t)
            B1(t) = betainc(P(t),(jj(t)+1)/2,nu(t)/2,'upper');
            B2(t) = betainc(P(t),(jj(t)+2)/2,nu(t)/2,'upper');
        end
        t = ~t;
        if any(t)
            B1(t) = betainc(Q(t),nu(t)/2,(jj(t)+1)/2,'lower');
            B2(t) = betainc(Q(t),nu(t)/2,(jj(t)+2)/2,'lower');
        end
    else
        if any(t)
            B1(t) = betainc(P(t),(jj(t)+1)/2,nu(t)/2,'lower');
            B2(t) = betainc(P(t),(jj(t)+2)/2,nu(t)/2,'lower');
        end
        t = ~t;
        if any(t)
            B1(t) = betainc(Q(t),nu(t)/2,(jj(t)+1)/2,'upper');
            B2(t) = betainc(Q(t),nu(t)/2,(jj(t)+2)/2,'upper');
        end
    end
    R1 = exp(gammaln((jj+1)/2+nu/2) - gammaln((jj+3)/2) - gammaln(nu/2) + ...
        ((jj+1)/2) .* log(P) + (nu/2) .* log(Q));
    R2 = exp(gammaln((jj+2)/2+nu/2) - gammaln((jj+4)/2) - gammaln(nu/2) + ...
        ((jj+2)/2) .* log(P) + (nu/2) .* log(Q));
    E10 = E1; E20 = E2; B10 = B1; B20 = B2; R10 = R1; R20 = R2; j0 = jj;
    todo = true(size(dsq));
    while(true)
        %Probability that t lies between 0 and x (x>0)
        twoterms = E1(todo).*B1(todo) + E2(todo).*B2(todo);
        subtotal(todo) = subtotal(todo) + twoterms;
        % Convergence test.
        todo(todo) = (abs(twoterms) > (abs(subtotal(todo))+seps)*seps);
        if (~any(todo))
            break;
        end
        
        % Update for next iteration
        jj = jj+2;
        
        E1(todo) = E1(todo) .* dsq(todo) ./ (jj(todo));
        E2(todo) = E2(todo) .* dsq(todo) ./ (jj(todo)+1);
        
        if uppertail==true
            B1(todo) = betainc(P(todo),(jj(todo)+1)/2,nu(todo)/2,'upper');
            B2(todo) = betainc(P(todo),(jj(todo)+2)/2,nu(todo)/2,'upper');
        else
            B1(todo) = B1(todo) - R1(todo);
            B2(todo) = B2(todo) - R2(todo);
            
            R1(todo) = R1(todo) .* P(todo) .* (jj(todo)+nu(todo)-1) ./ (jj(todo)+1);
            R2(todo) = R2(todo) .* P(todo) .* (jj(todo)+nu(todo)  ) ./ (jj(todo)+2);
        end
    end
    
    % Go back to the peak and start looping downward as far as necessary.
    E1 = E10; E2 = E20; B1 = B10; B2 = B20; R1 = R10; R2 = R20;
    jj = j0;
    todo = (jj>0);
    while any(todo)
        JJ = jj(todo);
        E1(todo) = E1(todo) .* (JJ  ) ./ dsq(todo);
        E2(todo) = E2(todo) .* (JJ+1) ./ dsq(todo);
        
        R1(todo) = R1(todo) .* (JJ+1) ./ ((JJ+nu(todo)-1) .* P(todo));
        R2(todo) = R2(todo) .* (JJ+2) ./ ((JJ+nu(todo))   .* P(todo));
        
        if uppertail==true
            B1(todo) = betainc(P(todo),(JJ-1)/2,nu(todo)/2,'upper');
            B2(todo) = betainc(P(todo),JJ/2,nu(todo)/2,'upper');
        else
            B1(todo) = B1(todo) + R1(todo);
            B2(todo) = B2(todo) + R2(todo);
        end
        
        twoterms = E1(todo).*B1(todo) + E2(todo).*B2(todo);
        subtotal(todo) = subtotal(todo) + twoterms;
        
        jj = jj - 2;
        todo(todo) = (abs(twoterms) > (abs(subtotal(todo))+seps)*seps) & ...
            (jj(todo) > 0);
    end
    p(k0) = min(1, max(0, p(k0) + subtotal/2));
end
