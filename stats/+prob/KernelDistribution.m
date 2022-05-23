classdef KernelDistribution < prob.FittableDistribution & prob.TruncatableDistribution & prob.ToolboxDistribution
%KernelDistribution Kernel-smoothed nonparametric probability distribution.
%    An object of the KernelDistribution class represents a probability
%    distribution estimated from a set of data by kernel smoothing.
%    You create it by fitting to data using the FITDIST function.
%
%    KernelDistribution methods:
%       cdf                   - Cumulative distribution function
%       icdf                  - Inverse cumulative distribution function
%       iqr                   - Interquartile range
%       mean                  - Mean
%       median                - Median
%       negloglik             - Negative log likelihood function
%       pdf                   - Probability density function
%       std                   - Standard deviation
%       random                - Random number generation
%       truncate              - Truncation distribution to an interval
%       var                   - Variance
%
%    KernelDistribution properties:    
%       DistributionName      - Name of the distribution
%       Kernel                - Name of kernel function
%       BandWidth             - Width of kernel
%       Truncation            - Two-element vector indicating truncation limits
%       IsTruncated           - Boolean flag indicating if distribution is truncated
%       InputData             - Structure containing data used to fit the distribution
%
%    See also fitdist, makedist, ksdensity.

%    Copyright 2012-2017 The MathWorks, Inc.

    properties(GetAccess=public, Constant)
%DistributionName Distribution name.
%    The DistributionName property indicates the name of the probability
%    distribution.
%
%    See also Kernel, BandWidth.
        DistributionName = getString(message('stats:dfittool:NameKernel'));
    end
    properties
%Kernel Kernel function name.
%    The Kernel property indicates the name of the kernel function used to
%    smooth the data. You can specify the kernel when you fit the distribution
%    using FITDIST. Valid names are 'normal' (default), 'box', 'triangle',
%    and 'epanechnikov'.
%
%    See also BandWidth, FITDIST.
        Kernel = 'normal'

%BandWidth Width of kernel function.
%    The BandWidth property is the width of the kernel smoothing window.
%    The default value is optimal for estimating normal densities, but a
%    smaller value may reveal features such as multiple modes. You can
%    specify the BandWidth when you fit the distribution using FITDIST.
%
%    See also Kernel, FITDIST.
        BandWidth = [];
    end
    properties(Hidden,Access='protected')
        ksinfo = []; % structure with kernel smoothing info
    end
    properties(Hidden)
        Support = [];
    end
   
    methods(Static,Hidden)
        function pd = fit(varargin)
%FIT Fit distribution to data.
%    FIT is a static method that fits the kernel distribution to data.
%    You should call the FITDIST function instead of calling this method
%    directly.
%
%    See also FITDIST.
            if nargin < 1
                error(message('stats:ProbDistUnivKernel:fit:TooFewInputs'));
            end
            
            pd = prob.KernelDistribution(varargin{:});
        end
    end
    methods(Hidden)
        function pd = KernelDistribution(x,varargin)
            narginchk(1,Inf);
            if ~isvector(x) || ~isnumeric(x) || size(x,2)~=1 || isempty(x)
                error(message('stats:ProbDistUnivKernel:BadX'));
            end
            
            % Process other arguments.
            okargs =   {'censoring' 'frequency' 'kernel' 'support'   'bandwidth' 'options'};
            defaults = {[]          []          'normal' 'unbounded' []          ''};
            
            [cens,freq,kernel,spt,width,~,~,extra] = ...
                internal.stats.parseArgs(okargs,defaults,varargin{:});
            if ~isempty(extra)
                width = internal.stats.parseArgs({'width'},{width},extra{:});
            end
            
            kernelnames = {'normal' 'epanechnikov'  'box'    'triangle'};
            kernel = internal.stats.getParamVal(kernel,kernelnames,'kernel');
            
            if ~isempty(cens) && ...
                    ~(   isequal(size(x),size(cens)) ...
                    && (islogical(cens) || all(ismember(cens,0:1))))
                error(message('stats:ProbDistUnivKernel:BadCens'));
            end
            if ~isempty(freq)
                if ~isvector(freq) || ~isnumeric(freq) || any(freq<0)
                    error(message('stats:ProbDistUnivKernel:BadFreq'))
                end
                if isscalar(freq)
                    freq = repmat(freq,size(x));
                elseif ~isequal(size(freq),size(x))
                    error(message('stats:ProbDistUnivKernel:BadFreq'))
                end
            end
            
            if ischar(spt) && size(spt,1)==1
                spt = internal.stats.getParamVal(spt,{'unbounded' 'positive'},'SUPPORT');
            elseif isequal(spt(:),[-Inf;Inf])
                spt = 'unbounded';
            elseif isequal(spt(:),[0;Inf])
                spt = 'positive';
            elseif ~isnumeric(spt) || numel(spt)~=2 ...
                    || ~all(isfinite(spt)) ...
                    || spt(1)>=spt(2)
                error(message('stats:ProbDistUnivKernel:BadSupport'));
            end
            if ~isempty(width)
                if ~isnumeric(width) || ~isscalar(width) ...
                        || ~isfinite(width) || width<=0
                    error(message('stats:ProbDistUnivKernel:BadWidth'));
                end
            end
            
            % Remove entries with NaN or with 0 frequency
            freq(freq==0) = NaN;
            [~,~,x,cens,freq]=internal.stats.removenan(x,cens,freq);
            if isempty(x)
                error(message('stats:ProbDistUnivKernel:BadReducedX'));
            end
            
            % Get ks info, including default width if needed
            if ~isempty(x)
                xi = x(1);
            elseif isnumeric(spt)
                xi = sum(spt)/2;
            else
                xi = 1;
            end
            [~,~,defaultwidth,ksinf] = ...
                ksdensity(x,xi,'cens',cens,'weight',freq,....
                'support',spt,'width',width);
            if isempty(width)
                width = defaultwidth;
            end
            
            % Fill in object properties
            pd.Kernel = kernel;
            pd.BandWidth = width;
            pd.Support = struct('range',spt, 'closedbound',[false false], 'iscontinuous',true);
            pd.InputData.data = x;
            pd.InputData.cens = cens;
            pd.InputData.freq = freq;
            pd.ksinfo = ksinf;
        end
    end
    methods(Access=protected)
        function y=icdffun(obj,p)
            %ICDF Inverse cumulative distribution function.
            
            % Check for valid input
            if nargin ~= 2
                error(message('stats:ProbDistUnivKernel:icdf:TooFewInputs'));
            end
            
            info = obj.ksinfo;
            y = dfswitchyard('statkscompute','icdf',p,true,length(p),obj.BandWidth, ...
                info.L,info.U,info.weight,[],obj.Kernel,...
                info.ty,obj.InputData.data,info.foldpoint,info.maxp);
        end
        function y=cdffun(obj,x,uflag)
            %CDF Cumulative distribution function.
            
            % Check for valid input
            if (nargin ~= 2) && (nargin ~= 3)
                error(message('stats:ProbDistUnivKernel:cdf:TooFewInputs'));
            end
            
            info = obj.ksinfo;
            y = dfswitchyard('statkscompute','cdf',x,true,length(x),obj.BandWidth, ...
                info.L,info.U,info.weight,[],obj.Kernel,...
                info.ty,[],info.foldpoint,info.maxp);
            if nargin==3
                if ~strcmpi(uflag,'upper')
                    error(message('stats:cdf:UpperTailProblem'));
                else
                    y = 1 - y;
                end           
            end
        end
        function y=pdffun(obj,x)
            %PDF Probability density function.
            
            % Check for valid input
            if nargin ~= 2
                error(message('stats:ProbDistUnivKernel:pdf:TooFewInputs'));
            end
            
            info = obj.ksinfo;
            y = dfswitchyard('statkscompute','pdf',x,true,length(x),obj.BandWidth, ...
                info.L,info.U,info.weight,[],obj.Kernel,...
                info.ty,[],info.foldpoint,info.maxp);
        end
        function y = randomfun(obj,varargin)
            %RANDOM Random number generation.
            
            info = obj.ksinfo;
            if info.maxp<1
                error(message('stats:ProbDistUnivKernel:LowProbability', num2str( info.maxp )))
            end
            
            % Randomly select a point where the ecdf jumps
            u = rand(varargin{:});
            sz = size(u);
            edges = min([0; cumsum(info.weight(:))],1); % guard histc against accumulated round-off;
            edges(end) = 1; % get the upper edge exact
            [~,bin] = histc(u(:),edges);
            y = info.ty(bin);
            y = reshape(y,sz);
            
            % Add random noise
            y = y + obj.BandWidth*noise(obj.Kernel,size(y));
            
            % Transform if necessary
            L = info.L;
            U = info.U;
            if L==0 && U==Inf
                y = exp(y);
            elseif isfinite(info.U)
                ey = exp(y);
                y = (U*ey + L) ./ (1 + ey);
            end
        end
        function displayCallback(this)
            if ~isempty(this)
                if ischar(this.Kernel)
                    fprintf('    %s = %s\n','Kernel',this.Kernel);
                else
                    fprintf('    %s = %s\n','Kernel',func2str(this.Kernel));
                end
                fprintf('    %s = %g\n','Bandwidth',this.BandWidth);
                if ischar(this.Support.range)
                    fprintf('    %s = %s\n','Support',this.Support.range);
                else
                    fprintf('    %s = (%g, %g)\n','Support',this.Support.range);
                end
            end
        end
    end
    
    methods(Static,Hidden)
        function info = getInfo
            info = getInfo@prob.ToolboxDistribution('prob.KernelDistribution');
            info.name = prob.KernelDistribution.DistributionName;
            info.code = 'kernel';
            info.censoring = true;
         end
    end
    
    properties(Hidden,Dependent=true,GetAccess='public',SetAccess='protected')
        % For compatibility with earlier versions of probability objects
        DistName
        NLogL
    end
    methods % get/set methods
        function a = get.DistName(this)
            info = this.getInfo;
            a = info.code;  % code, not name, for backward compatibility
        end
        function a = get.NLogL(this)
            a = negloglik(this);
        end
    end
    methods
        function m = mean(this)
            requireScalar(this)
            requireFull(this) % require no censoring at largest observation
            m = meanvar(this);
        end
        function v = var(this)
            requireScalar(this)
            requireFull(this)
            [~,v] = meanvar(this);
        end
    end
    methods(Access=protected)
        function [m,v] = meanvar(this)
            % Get saved info from ksinfo property
            L = this.ksinfo.L;           % lower limit
            U = this.ksinfo.U;           % upper limit
            ty = this.ksinfo.ty(:);      % y transformed, L<=ty<=U
            wt = this.ksinfo.weight(:)'; % weights
            bw = this.BandWidth;         % kernel width
            kvar = bw^2;     % variance is the square of the width because
                             % all kernels are normalized to variance 1
            if this.IsTruncated || U<Inf || (L~=0 && L~=-Inf)
                % distribution defined between fixed finite limits - handle
                % by numerical integration as for a truncated distribution
                m = truncatedMoment(this,1);
                if nargout>=2
                    v = truncatedMoment(this,2);
                end
            elseif L==-Inf && U==Inf
                % unbounded distribution
                m = wt*ty;
                if nargout>=2
                    bias = ty-m;
                    v = wt*(kvar + bias.^2);
                end
            else   % L==0 && U==Inf
                % positive distribution
                switch(this.Kernel)
                    case 'normal'
                        obsmean = exp(ty + .5*kvar);
                        m = wt*obsmean;
                        if nargout>=2
                            obsvar = exp(2*ty + kvar) .* (exp(kvar)-1);
                            bias = obsmean - m;
                            v = wt*(obsvar + bias.^2);
                        end
                    case 'triangle'
                        h = bw*sqrt(6);
                        obsmean = exp(ty)*(exp(h)+exp(-h)-2)/h^2;
                        m = wt*obsmean;
                        if nargout>=2
                            obsvar = (64*sinh(h/2)^2*exp(2*ty) + 64*sinh(h/2)^4*exp(2*ty) ...
                            - 2*h^2*exp(2*ty) - 32*sinh(h/2)^2*exp(2*ty)*exp(h) ...
                            - 32*sinh(h/2)^2*exp(-h)*exp(2*ty) + h^2*exp(-2*h)*exp(2*ty) ...
                            + h^2*exp(2*h)*exp(2*ty) + 2*h^3*exp(-2*h)*exp(2*ty) ...
                            - 2*h^3*exp(2*h)*exp(2*ty) + 4*h^2*ty.*exp(2*ty) ...
                            - 2*h^2*ty*exp(-2*h).*exp(2*ty) - 2*h^2*ty*exp(2*h).*exp(2*ty) ...
                            + 4*h^3*exp(2*ty)*exp(h)*sinh(h) + 4*h^3*exp(-h)*exp(2*ty)*sinh(h) ...
                            + 4*h^2*ty.*exp(2*ty)*exp(h)*sinh(h) - 4*h^2*ty*exp(-h).*exp(2*ty)*sinh(h))/(4*h^4);
                            bias = obsmean - m;
                            v = wt*(obsvar + bias.^2);
                        end

                    case 'box'
                        h = bw*sqrt(3);
                        obsmean = (1/(2*h)) * exp(ty) * (exp(h)-exp(-h));
                        m = wt*obsmean;
                        if nargout>=2
                            obsvar = (1/(4*h^2))*(exp(2*ty-2*h)*(exp(2*h)-1)...
                                *(h - exp(2*h) + h*exp(2*h) + 1));
                            bias = obsmean - m;
                            v = wt*(obsvar + bias.^2);
                        end
                    case 'epanechnikov'
                        h = bw*sqrt(5);
                        obsmean = (3*exp(-h)*exp(ty)*(h - exp(2*h) + h*exp(2*h) + 1))/(2*h^3);
                        m = wt*obsmean;
                        if nargout>=2
                            obsvar = -(36*exp(2*ty - 2*h) - 72*exp(2*ty) ...
                                + 36*exp(2*h + 2*ty) + 72*h*exp(2*ty - 2*h) ...
                                - 72*h*exp(2*h + 2*ty) + 72*h^2*exp(2*ty) ...
                                + 36*h^2*exp(2*ty - 2*h) + 36*h^2*exp(2*h + 2*ty) ...
                                - 3*h^3*exp(2*ty - 2*h) + 3*h^3*exp(2*h + 2*ty) ...
                                - 6*h^4*exp(2*ty - 2*h) - 6*h^4*exp(2*h + 2*ty))/(16*h^6);
                            bias = obsmean - m;
                            v = wt*(obsvar + bias.^2);
                        end
                end
            end
        end
    end
end % classdef

% ----------------
function requireFull(this)
if this.ksinfo.maxp<1
    error(message('stats:probdists:NotFullDistribution'))
end
end
function z = noise(kernel,sz)
% Generates random values for kernels defined in the function
% stats/private/statkskernelinfo

switch(kernel)
    case 'normal'
        z = randn(sz);
    case 'box'
        z = -sqrt(3) + 2*sqrt(3)*rand(sz);
    case 'triangle'
        u = -1 + 2*rand(sz);
        z = sqrt(6) * sign(u) .* (1-sqrt(abs(u)));
    case 'epanechnikov'
        % Need to compute the icdf by solving the quadratic equation
        %   x^2 + (2-4p)x + 1 = 0
        b = 2 - 4*rand(sz);
        
        % Get one root
        r = (-b + sqrt(b.^2 - 4)) / 2;
        
        % The desired result is z=(s-t) where s and t are defined as
        % below.  Get t using the right cube root of r.  Compute s from t.
        % Remove any imaginary part of z caused by roundoff error, and
        % scale by the standard kernel width sqrt(5).
        t = exp(2*1i*2*pi/3) * r.^(1/3);
        s = -1./t;
        z = sqrt(5) * real(s - t);
end
end

