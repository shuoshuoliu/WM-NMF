function [IncludeTies,K,NSMethod,Distance,AdditionalArg,BucketSize] = validateNSParams(varargin)
%#codegen
% VALIDATENSPARAMS helper function to parse and validate inputs for
% knnsearch and rangesearch.

%   Copyright 2017-2018 The MathWorks, Inc.

coder.inline('always');
coder.internal.prefer_const(varargin);
% MATLAB path needed for ML object to codegen object redirection, which
% happens in MATLAB.
if coder.target('MATLAB') 
    [includeties,k,nsmethod,distance,additionalarg,bucketsize] = parseOptionalInputsML(varargin{:});
else
    
    [includeties,k,nsmethod,distance,additionalarg,bucketsize] = parseOptionalInputsCG(varargin{:});
end

IncludeTies = includeties;
K = k;
NSMethod = nsmethod;
Distance = distance;
AdditionalArg = additionalarg;
BucketSize = bucketsize;

end

function [includeties,k,nsmethod,distance,additionalarg,bucketsize] = parseOptionalInputsCG(varargin)
% PARSEOPTIONALINPUTS  Parse optional PV pairs

coder.inline('always');
coder.internal.prefer_const(varargin);
params = struct( ...
    'IncludeTies', uint32(0), ...
    'K',     uint32(0),...
    'NSMethod',     uint32(0),...
    'Distance',     uint32(0),...
    'P', uint32(0),...
    'Cov',uint32(0),...
    'Scale',uint32(0),...
    'DistExtra',uint32(0),...
    'DistParamExtra',uint32(0),...
    'XSize',uint32(0),...
    'BucketSize',uint32(0));

popts = struct( ...
    'CaseSensitivity', false, ...
    'StructExpand',    false, ...
    'PartialMatching', 'unique',...
    'IgnoreNulls',true);

optarg = eml_parse_parameter_inputs(params, popts, varargin{:});

includeties = eml_get_parameter_value(optarg.IncludeTies,false, varargin{:});
coder.internal.errorIf(~coder.internal.isConst(includeties),...
    'stats:coder:pdist2:ExpectedConstant','IncludeTies');
validateattributes(includeties,{'logical'},{'nonempty','scalar'},mfilename,'IncludeTies');
% includeTies needs to be a compile-time constant

k = eml_get_parameter_value(optarg.K, coder.internal.indexInt(1), varargin{:});
coder.internal.errorIf(~coder.internal.isConst(isscalar(k)) || ~isscalar(k),...
    'stats:coder:pdist2:ExpectedScalar','K');
validateattributes(k,{'numeric'},{'nonempty','integer','positive','nonnan','finite','real'},mfilename,'K');

NSMethod = eml_get_parameter_value(optarg.NSMethod, '', varargin{:});
% coder.internal.errorIf(~coder.internal.isConst(NSMethod) || ~strcmpi(NSMethod,...
%     'exhaustive'),'stats:coder:knnsearch:KDTreeNotAllowed');

NSMethodIsOK = true;
nsmethodSpecified = true;
if strcmpi(NSMethod, 'exhaustive')
    nsmethod = EXHAUSTIVE;
elseif strcmpi(NSMethod, 'kdtree')
    nsmethod = KDTREE;
elseif isempty(NSMethod)
    % Decide which searcher to use
    nsmethodSpecified = false;
    nsmethod = EXHAUSTIVE;
else
    nsmethod = EXHAUSTIVE;
    NSMethodIsOK = false;
end

coder.internal.errorIf(~coder.internal.isConst(NSMethod) || ~NSMethodIsOK,...
     'stats:coder:knnsearch:UnsupportedNSMethod'); 

bucketsize = eml_get_parameter_value(optarg.BucketSize,[], varargin{:});

validateattributes(bucketsize,{'numeric'},{'integer','positive','nonnan','finite','real'},mfilename,'BucketSize');

disttemp = eml_get_parameter_value(optarg.Distance, '', varargin{:});
p = eml_get_parameter_value(optarg.P,[], varargin{:});
cov = eml_get_parameter_value(optarg.Cov,[], varargin{:});
scale = eml_get_parameter_value(optarg.Scale,[], varargin{:});
xsize = eml_get_parameter_value(optarg.XSize,1, varargin{:});
distextra = eml_get_parameter_value(optarg.DistExtra,'', varargin{:});
distparamextra = eml_get_parameter_value(optarg.DistParamExtra,[], varargin{:});

if isempty(disttemp)
    if isempty(distextra)
        dist = 'euclidean';
    else
        dist = distextra;
    end
else
    dist  = disttemp;
end
distance = stats.coder.distutils.validateDistance(dist);

% If nsmethod is not provided, decide whcih one to use
if ~nsmethodSpecified
    isKdtreeDistance = strncmpi(distance,'euc',3) || strncmpi(distance,'cit',3) || ...
        strncmpi(distance,'che',3) || strncmpi(distance,'min',3);
    if isKdtreeDistance && xsize <= 7
        nsmethod = KDTREE;
    else
        nsmethod = EXHAUSTIVE;
        % Error out if BucketSize is specified
        coder.internal.errorIf(~isempty(bucketsize),'stats:coder:knnsearch:bSizeNotNeededEXH');
    end
else
    % Error out if the specified NSMethod is exhaustive and BucketSize is also specified
    bsizeNotNeeded = ~isempty(bucketsize) && nsmethod == EXHAUSTIVE;
    coder.internal.errorIf(bsizeNotNeeded,'stats:coder:knnsearch:bSizeNotNeeded');
end


% If any pair of p, cov and scale inputs are known to be not compile-time empty
% variables, error out.
if ~ (coder.internal.isConst(isempty(p)) && isempty(p)) 
    coder.internal.errorIf(~coder.internal.isConst(isempty(cov)) || ~isempty(cov) ||...
        ~coder.internal.isConst(isempty(scale)) || ~isempty(scale),'stats:coder:knnsearch:AmbiguousAdditionalArg');
elseif ~( coder.internal.isConst(isempty(cov))  && isempty(cov))
    coder.internal.errorIf(~coder.internal.isConst(isempty(p)) || ~isempty(p) ||...
        ~coder.internal.isConst(isempty(scale)) || ~isempty(scale),'stats:coder:knnsearch:AmbiguousAdditionalArg');
elseif ~( coder.internal.isConst(isempty(scale))  && isempty(scale))
    coder.internal.errorIf(~coder.internal.isConst(isempty(cov)) || ~isempty(cov) ||...
        ~coder.internal.isConst(isempty(p)) || ~isempty(p),'stats:coder:knnsearch:AmbiguousAdditionalArg');
end

% At most one of p, cov or scale arguments is nonempty. If all are empty,
% use extra argument for knnsearch/rangesearch methods. 

if ~isempty(p)
    coder.internal.errorIf(~strncmpi(distance,'min',3),'stats:ExhaustiveSearcher:InvalidMinExp'); 
    additionalarg = p;  
elseif ~isempty(cov)
    coder.internal.errorIf(~strncmpi(distance,'mah',3),'stats:ExhaustiveSearcher:InvalidCov'); 
    additionalarg = cov;
elseif ~isempty(scale)
    coder.internal.errorIf(~strncmpi(distance,'seu',3),'stats:ExhaustiveSearcher:InvalidScale');  
    additionalarg = scale;
elseif ~isempty(distparamextra)
    switch distance
        case 'minkowski'
            if strncmpi(distextra,'min',3)
                additionalarg = distparamextra;
            else
                additionalarg = [];
            end
        case 'seuclidean'
            if  strncmpi(distextra,'seu',3)
                additionalarg = distparamextra;
            else
                additionalarg = [];
            end
        case 'mahalanobis'
            if  strncmpi(distextra,'mah',3)
                additionalarg = distparamextra;
            else
                additionalarg = [];
            end
        otherwise
            additionalarg = [];
    end
else
    additionalarg = [];
end

if ~isempty(additionalarg)
    stats.coder.distutils.validateDistParameter(distance,additionalarg,xsize);
end

end


function [includeties,k,nsmethod,distance,additionalarg,bucketsize] = parseOptionalInputsML(varargin)
% PARSEOPTIONALINPUTS  Parse optional PV pairs

parser = inputParser;
addParameter(parser,'Distance','');
addParameter(parser,'K',1);
addParameter(parser,'P',[]);
addParameter(parser,'Cov',[]);
addParameter(parser,'Scale',[]);
addParameter(parser,'IncludeTies',false);
addParameter(parser,'NSMethod','exhaustive');
addParameter(parser,'XSize',1);
addParameter(parser,'DistExtra','');
addParameter(parser,'DistParamExtra',[]);
addParameter(parser,'BucketSize',[]);

parse(parser,varargin{:});
validateattributes(parser.Results.K,{'numeric'},{'nonempty','scalar','integer','positive','nonnan','finite','real'},mfilename,'K');
k = parser.Results.K;
validateattributes(parser.Results.IncludeTies,{'logical'},{'nonempty','scalar'},mfilename,'IncludeTies');
includeties = parser.Results.IncludeTies;
distextra = parser.Results.DistExtra;
if isempty(parser.Results.Distance)
    if isempty(distextra)
        dist = 'euclidean';
    else
        dist = distextra;
    end
else
    dist  = parser.Results.Distance;
end
distance = stats.coder.distutils.validateDistance(dist);

validateattributes(parser.Results.BucketSize,{'numeric'},{'integer','positive','nonnan','finite','real'},mfilename,'BucketSize');
bucketsize = parser.Results.BucketSize;

NSMethodParsed =  parser.Results.NSMethod;

NSMethodIsOK = true;
if strcmpi(NSMethodParsed, 'exhaustive')
    nsmethod = EXHAUSTIVE;
elseif strcmpi(NSMethodParsed, 'kdtree')
    nsmethod = KDTREE;
else
    nsmethod = EXHAUSTIVE;
    NSMethodIsOK = false;
end

if ~NSMethodIsOK
    error(message('stats:coder:knnsearch:UnsupportedNSMethod'));
end

bsizeNotNeeded = ~isempty(bucketsize) && nsmethod == EXHAUSTIVE;
if bsizeNotNeeded
	error(message('stats:coder:knnsearch:bSizeNotNeeded'));
end

xsize = parser.Results.XSize;
p = parser.Results.P;
cov = parser.Results.Cov;
scale = parser.Results.Scale;
distparamextra = parser.Results.DistParamExtra;

% Error out if NSMethod is kdtree and cov or scale are specified
covScaleNotNeeded = nsmethod == KDTREE && (~isempty(scale) || ~isempty(cov));
if covScaleNotNeeded
	error(message('stats:coder:knnsearch:covScaleNotNeeded'));
end

pnotempty = false;
covnotempty = false;
scalenotempty = false;
distparamextranotempty = false;
NOPTS = uint8(~isempty(p)) + uint8(~isempty(scale)) + uint8(~isempty(cov));
if NOPTS > 1
    error(message('stats:coder:knnsearch:AmbiguousAdditionalArg'));
end
if ~isempty(p)
    pnotempty = true;
end
if ~isempty(scale)
    scalenotempty = true;
end
if ~isempty(cov)
    covnotempty = true;
end
if ~isempty(distparamextra)
    distparamextranotempty = true;
end
if pnotempty
    if ~strncmpi(distance,'min',3)
        error(message('stats:ExhaustiveSearcher:InvalidMinExp'));
    end
    additionalarg = p;
elseif covnotempty
    if ~strncmpi(distance,'mah',3)
        error(message('stats:ExhaustiveSearcher:InvalidCov'));
    end    
    additionalarg = cov;
elseif scalenotempty
    if ~strncmpi(distance,'seu',3)
        error(message('stats:ExhaustiveSearcher:InvalidScale'));
    end    
    additionalarg = scale;
elseif distparamextranotempty
    switch distance
        case 'minkowski'
            if strncmpi(distextra,'min',3)
                additionalarg = distparamextra;
            else
                additionalarg = [];
            end
        case 'seuclidean'
            if  strncmpi(distextra,'seu',3)
                additionalarg = distparamextra;
            else
                additionalarg = [];
            end
        case 'mahalanobis'
            if  strncmpi(distextra,'mah',3)
                additionalarg = distparamextra;
            else
                additionalarg = [];
            end
        otherwise
            additionalarg = [];
    end    
else
    additionalarg = [];
end

if ~isempty(additionalarg)
    stats.coder.distutils.validateDistParameter(distance,additionalarg,xsize);
end

end


function k = EXHAUSTIVE
coder.inline('always');
k = int8(1);
end

 function k = KDTREE
 coder.inline('always');
 k = int8(2);
 end