function logy = getGLMLogPDF(distr,x,p1,p2)
%getGLMLogPDF Compute log pdf for generalized linear models.

%   Copyright 2016 The MathWorks, Inc.

switch(distr)
    case 'binomial'
        n = p1;
        p = p2;
        logy = gammaln(n+1) - gammaln(x+1) - gammaln(n-x+1) + x.*log(p) + (n-x).*log1p(-p);
    case 'gamma'
        a = p1;
        b = p2;
        z = x./b;
        logy = (a-1).*log(z) - z - gammaln(a) - log(b);
    case 'inverse gaussian'
        mu = p1;
        lambda = p2;
        logy = .5*log(lambda./(2.*pi.*x.^3)) + (-0.5.*lambda.*(x./mu - 2 + mu./x)./mu);
    case 'normal'
        mu = p1;
        sigma = p2;
        logy = (-0.5 * ((x - mu)./sigma).^2) - log(sqrt(2*pi) .* sigma);
    case 'poisson'
        lambda = p1;
        
        logy = (-lambda + x.*log(lambda) - gammaln(x+1));
end
end