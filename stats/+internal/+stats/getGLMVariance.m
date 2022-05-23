function [link,estdisp,sqrtvarFun,devFun,muLims] = ...
    getGLMVariance(distr,estdisp,link,N,y,dataClass)
%getGLMVariance Get variance function and related things for fitglm.

% sqrtvarFun is a function that produces the square root of the variance
%            function
% devFun     is the deviance function
% link       is the link function; if 'canonical' on input it is set here
% estdisp    is set to true for distributions that require it
% muLims     is used to enforce limits on mu to guard against an inverse
%            link that doesn't map into the support of the distribution.
%            binomial: mu is a probability, so order one is the natural
%                      scale, and eps is a reasonable lower limit
%            other non-normal: we don't know the natural scale for mu, so
%                      make the lower limit small to keep mu^4 from
%                      underflowing. No upper limit.

%   Copyright 2016-2020 The MathWorks, Inc.

if isempty(estdisp)
    % Prepare default for binomial and Poisson, others will override
    estdisp = false;
end
try
    switch lower(distr)
        case 'normal'
            link = chooseLink(link,'identity');
            estdisp = true;
            if nargout>=3
                sqrtvarFun = @(mu) ones(size(mu),"like",mu);
                devFun = @(mu,y) (y - mu).^2;
                muLims = [];
            end
        case 'binomial'
            link = chooseLink(link,'logit');
            if nargout>=3
                if isempty(N), N = 1; end
                sqrtN = sqrt(N);
                sqrtvarFun = @(mu) sqrt(mu).*sqrt(1-mu) ./ sqrtN;
                devFun = @(mu,y) 2*N.*(y.*log((y+(y==0))./mu) + (1-y).*log((1-y+(y==1))./(1-mu)));
                muLims = [eps(dataClass) 1-eps(dataClass)];
            end
        case 'poisson'
            link = chooseLink(link,'log');
            if nargout>=3
                if any(y < 0)
                    error(message('stats:glmfit:BadDataPoisson'));
                end
                if nargout>=3
                    sqrtvarFun = @(mu) sqrt(mu);
                    devFun = @(mu,y) 2*(y .* (log((y+(y==0)) ./ mu)) - (y - mu));
                    muLims = realmin(dataClass).^.25;
                end
            end
        case 'gamma'
            link = chooseLink(link,'reciprocal');
            estdisp = true;
            if nargout>=3
                if any(y <= 0)
                    error(message('stats:glmfit:BadDataGamma'));
                end
                sqrtvarFun = @(mu) mu;
                devFun = @(mu,y) 2*(-log(y ./ mu) + (y - mu) ./ mu);
                muLims = realmin(dataClass).^.25;
            end
        case 'inverse gaussian'
            link = chooseLink(link,-2);
            estdisp = true;
            if nargout>=3
                if any(y <= 0)
                    error(message('stats:glmfit:BadDataInvGamma'));
                end
                sqrtvarFun = @(mu) mu.^(3/2);
                devFun = @(mu,y) ((y - mu)./mu).^2 ./  y;
                muLims = realmin(dataClass).^.25;
            end
        otherwise
            error(message('stats:glmfit:BadDistribution'));
    end
catch ME
    throwAsCaller(ME)
end

function link = chooseLink(link,clink)
if isequal(link, 'canonical') || isempty(link)
    link = clink;
end
