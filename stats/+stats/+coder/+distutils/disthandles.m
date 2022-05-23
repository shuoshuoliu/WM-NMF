function str = disthandles(metric,varargin) %#codegen
%MATLAB Code Generation Private Function
% Helper for pdist/pdist2
%DISTHANDLES returns function handles required for distance computation
% distfun is elementwise function handle for a distance metric, these
% function handles can be used for all metrics except for the mahalanobis
% metric, the calculation for which is more complex, inside of the for
% loops

% postfun is elementwise function handle that works on the entire 
%   Copyright 2018 The MathWorks, Inc.
coder.inline('always');
postprocess = false;
twosums = false;

if strncmpi(metric,'min',3)
    param = varargin{1};
elseif strncmpi(metric,'ham',3)
    param = varargin{1};
end
    
switch metric
    
    case {'euc','seu'}
        distfun = @(x,y) (x-y)*(x-y);
        postprocess = true;
        postfun = @(temp) sqrt(temp);
        
    case 'squ'
        distfun = @(x,y) (x-y)*(x-y);
        
    case 'cit'
        distfun = @(x,y) abs(x-y);
        
    case 'min'
        distfun = @(x,y) abs(x-y)^param;
        postprocess = true;
        postfun = @(temp) temp.^(1/param);
        
    case {'cos','cor','spe'}
        distfun = @(x,y) x*y;
        postprocess = true;
        postfun = @(temp) (1-temp).*(temp<1) ;
        
    case 'ham'
        distfun = @(x,y) (x~=y);
        postprocess = true;
        postfun = @(temp) temp./param;
        
    case 'jac'
        twosums = true;
        distfun = @(x,y) ((x~=y) & (x~=0 | y~=0));
        distfunden = @(x,y) (x~=0 | y~=0);
end

str.distfun = distfun;

if twosums
    str.distfunden = distfunden;
else
    str.distfunden = [];
end

if postprocess
    str.postfun = postfun;
else
    str.postfun = [];
end
end
