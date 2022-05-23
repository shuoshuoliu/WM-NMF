function [fout0,xout,u,ksinfo] = ksdensity(yData,varargin)
%KSDENSITY Compute one or two dimensional kernel density or distribution
%estimate
%   [F,XI]=KSDENSITY(X) computes a probability density estimate of the
%   sample in the vector or two-column matrix X. KSDENSITY evaluates the
%   density estimate at 100 points (900 points for bivariate data) covering
%   the range of the data. XI is the set of 100 (or 900) points. For
%   bivariate data, XI is created using MESHGRID from 30 equally spaced
%   points in each dimension. F is the vector of density values. The
%   estimate is based on the normal kernel function, using the window
%   parameter (bandwidth) that is a function of the number of points and
%   dimension in X.
%
%   F=KSDENSITY(X,XI) specifies the vector (or two-column matrix XI for
%   bivariate data) of values where the density estimate is to be
%   evaluated.
%
%   [F,XI,U]=KSDENSITY(...) also returns the bandwidth(s) of the kernel
%   smoothing window(s).
%
%   KSDENSITY(...) without output arguments produces a plot of the results.
%
%   KSDENSITY(AX,...) plots into axes AX instead of GCA.
%
%   [...]=KSDENSITY(...,'PARAM1',val1,'PARAM2',val2,...) specifies
%   parameter name/value pairs to control the density estimation.  Valid
%   parameters are the following:
%
%      Parameter    Value
%      'Kernel'     The type of kernel smoother to use, chosen from among
%                   'normal' (default), 'box', 'triangle', and
%                   'epanechnikov'. For Bivariate data, the same kernel is
%                   applied to each dimension.
%      'Support'    Either 'unbounded' (default) if the density can extend
%                   over the whole real line, or 'positive' to restrict it
%                   to positive values, or a two-element vector for
%                   univariate data giving lower and upper limits for the
%                   support of the density, or a two-by-two matrix for
%                   bivariate data, where the first row contains the lower
%                   limits and the second row contains the upper limits.
%      'Weights'    Vector of the same length as X, giving the weight to
%                   assign to each X value (default is equal weights).
%      'Bandwidth'  A scalar value of the bandwidth of the kernel smoothing
%                   window for each dimension. For bivariate data, it can
%                   also be a 2-element vector.  The default is optimal for
%                   estimating normal densities, but you may want to choose
%                   a smaller value to reveal features such as multiple
%                   modes.
%      'Function'   The function type to estimate, chosen from among 'pdf',
%                   'cdf', or 'survivor'for the density, cumulative
%                   probability, or survivor, respectively. For univariate
%                   data, the function type can also be 'icdf' for inverse
%                   cumulative probability or 'cumhazard' for cumulative
%                   hazard functions. For 'icdf',
%                   F=KSDENSITY(X,YI,...,'function','icdf') computes the
%                   estimated inverse CDF of the values in X, and evaluates
%                   it at the probability values specified in YI.
%      'BoundaryCorrection'  The method for boundary correction for 'pdf'
%                   kernel density estimation with bounded support. Choices
%                   are 'log' and 'reflection'. The 'log' method first
%                   applies a log transformation to convert bounded data to
%                   be unbounded, then transforms back to the original
%                   scale after performing the density estimation. The
%                   'reflection' method adds the reflection of all points
%                   at the boundaries, then computes the density estimate
%                   inside the bounded range including contributions from
%                   the augmented data. Default is 'log'.
%
%   The parameters below are only valid for one dimensional data:
% 
%      'Censoring'  A logical vector of the same length of X, indicating
%                   which entries are censoring times (default is no
%                   censoring).
%      'NumPoints'  The number of equally-spaced points in XI. 
%
%   The parameter below is only valid for two dimensional data:
% 
%      'PlotFcn'    The function used to create the kernel density plot,
%                   chosen from 'surf'(default), 'contour', 'plot3', and
%                   'surfc'.
%
%   In place of the kernel functions listed above, you can specify another
%   kernel function by using @ (such as @normpdf) or quotes (such as
%   'normpdf'). KSDENSITY calls the function with a single argument that is
%   an array containing distances between data values in X and locations in
%   XI where the density is evaluated normalized by the bandwidth.  The
%   function must return an array of the same size containing corresponding
%   values of the kernel function. When the 'function' parameter value is
%   'pdf', this kernel function should return density values, otherwise it
%   should return cumulative probability values. Specifying a custom kernel
%   when the 'function' parameter value is 'icdf' is an error.
%
%   For univariate data, if the 'support' parameter is 'positive',
%   KSDENSITY transforms X using a log function, estimates the density of
%   the transformed values, and transforms back to the original scale.  If
%   'support' is a vector [L U], KSDENSITY uses the transformation
%   log((X-L)/(U-X)). For bivariate data, 'support' can also be a
%   combination of positive, unbounded, or bounded variables specified as
%   [0 -Inf;Inf Inf] or [0 L;Inf U], and KSDENSITY transforms each
%   dimension of X in the same way as the univariate data.   The
%   'bandwidth' parameter and U outputs are on the scale of the transformed
%   values.
%
%   Example: generate a mixture of two normal distributions, and
%   plot the estimated density.
%      x = [randn(30,1); 5+randn(30,1)];
%      [f,xi] = ksdensity(x);
%      plot(xi,f);
%
%   Example: generate a mixture of two normal distributions, and plot the
%   estimated cumulative distribution at a specified set of values.
%      x = [randn(30,1); 5+randn(30,1)];
%      xi = linspace(-10,15,201);
%      f = ksdensity(x,xi,'function','cdf');
%      plot(xi,f);
%
%   Example: generate a mixture of two normal distributions, and plot the
%   estimated inverse cumulative distribution function at a specified set
%   of values.
%      x = [randn(30,1); 5+randn(30,1)];
%      yi = linspace(.01,.99,99);
%      g = ksdensity(x,yi,'function','icdf');
%      plot(yi,g);
%
%   Example: generate a mixture of two bivariate normal distributions,
%   and plot the estimated density.
%      gridx1 = -0.25:.05:1.25;
%      gridx2 = 0:.1:15;
%      [x1,x2] = meshgrid(gridx1, gridx2);
%      x1 = x1(:);
%      x2 = x2(:);
%      xi = [x1 x2];
%      X = [0+.5*rand(20,1) 5+2.5*rand(20,1);
%           .75+.25*rand(10,1) 8.75+1.25*rand(10,1)];
%      ksdensity(X,xi);
%
%   Example: generate a uniform distribution within range [0, 1], and plot 
%   the estimated density with reflection boundary correction
%       x = rand(1000, 1);
%       ksdensity(x, linspace(0,1,1000), 'support', [0,1], ...
%           'BoundaryCorrection', 'Reflection');
%
%   Example: generate samples by taking positive values from a mixture of 
%   two bivariate normal distributions, and plot the estimated density 
%   with reflection boundary correction
%       rng(1)
%       mu = [0, 0];
%       sigma = [1, 0.5; 0.5, 1];
%       t = mvnrnd(mu, sigma, 10000);
%       x = t(t(:,1)>0 & t(:,2)>0,:);
%       gridx1 = 0:0.03:max(x(:,1));
%       gridx2 = 0:0.03:max(x(:,2));
%       [x1, x2] = meshgrid(gridx1, gridx2);
%       x1 = x1(:);
%       x2 = x2(:);
%       xi = [x1, x2];
%       ksdensity(x, xi, 'support', 'positive', 'BoundaryCorrection',...
%           'Reflection');
%   See also MVKSDENSITY, MESHGRID, HIST, @.

%   Undocumented KSINFO output is used for the ProbDistUnivKernel class and
%   is subject to change in a future release.

% If there is any censoring, we would like to estimate the density up
% to the last non-censored observation.  Say this is XMAX.  Without
% censoring, the density estimate near XMAX would consist of contributions
% from kernels centered above and below XMAX.  We can't compute the
% contributions above XMAX, though, because we have no data.  Using only
% the kernels centered below XMAX makes the density estimate biased.
%
% In an attempt to reduce bias, we will compute the contributions
% from kernels centered below XMAX, and fold their values around XMAX.
% The result should be good if the density is nearly flat in this area.
% If the density is increasing then the estimate will still be biased
% downward, and if the density is decreasing it will still be biased
% upward, but the bias will be reduced.

% Reference:
% [1]  A.W. Bowman and A. Azzalini (1997), "Applied Smoothing
%      Techniques for Data Analysis," Oxford University Press.
% [2]  B.W. Silverman (1998), "Density Estimation for Statistics and Data 
%      Analysis", Chapman & Hall/CRC.

%   Copyright 1993-2017 The MathWorks, Inc.


if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

[axarg,varargin] = axescheck(yData,varargin{:});
if isempty(axarg)
    axarg = {};
else
    axarg = {axarg};
end
if isempty(varargin)
    yData = [];
else
    yData = varargin{1};
    varargin(1) = [];
end

if isvector(yData) 
    d = 1;
    yData = yData(:);
else
    d = size(yData,2);
end

if isempty(yData) || ~isvector(yData) && ~(ismatrix(yData) && d==2)
    error(message('stats:ksdensity:BadX'));
end

xi = [];
if ~isempty(varargin)
    if ~ischar(varargin{1})
        xi = varargin{1};
        varargin(1) = [];
    end
end

if nargout>=4
    [fout,xout,u,plottype,ksinfo] = mvksdensity(yData,xi,varargin{:});
else
    [fout,xout,u,plottype] = mvksdensity(yData,xi,varargin{:});
end

% Plot the results if they are not requested as return values
if nargout==0
    if d==1
        plot(axarg{:},xout,fout);
    else %d==2
        x1 = xout(:,1);
        x2 = xout(:,2);
        switch plottype           
            case 'contour'
                try
                    [xq,yq,z] = computeGrid(x1,x2,fout);
                    contour(axarg{:},xq,yq,z);
                catch
                    plot3(axarg{:},x1,x2,fout,'x');
                end
                
            case 'plot3'
                plot3(axarg{:},x1,x2,fout,'x');
            case 'surf'
                try
                    [xq,yq,z] = computeGrid(x1,x2,fout);
                    surf(axarg{:},xq,yq,z);
                catch
                    plot3(axarg{:},x1,x2,fout,'x'); % collinear points
                end
            case 'surfc'
                try
                    [xq,yq,z] = computeGrid(x1,x2,fout);
                    h = surfc(axarg{:},xq,yq,z);
                    zt = get(gca,'ZTick');
                    h(2).ContourZLevel = 2*zt(1) - zt(2);
                catch
                    plot3(axarg{:},x1,x2,fout,'x'); % collinear points
                end               
        end
    end
else
    fout0 = fout;
end

function [xq,yq,z] = computeGrid(x1,x2,fout)
x = linspace(min(x1),max(x1));
y = linspace(min(x2),max(x2));
[xq,yq] = meshgrid(x,y);
orig_state = warning;
warning('off','all');
z = griddata(x1,x2,fout,xq,yq);
warning(orig_state);
