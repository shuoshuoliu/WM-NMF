function [includeties,k,distance,p,cov,scale,nsmethod,bucketsize] = extractSearchParams(isrange,varargin)
%#codegen
% extractSearchParams helper function to parse inputs for
% knnsearch and rangesearch.

% Copyright 2017-2018 The MathWorks, Inc.

coder.inline('always');
coder.internal.prefer_const(isrange,varargin);
% MATLAB path needed for ML object to codegen object redirection, which
% happens in MATLAB.
[includeties,k,distance,p,cov,scale] = stats.coder.distutils.extractNSParams(varargin{:});
if coder.target('MATLAB')
    [nsmethod,bucketsize] = parseOptionalInputsML(isrange,varargin{:});
else
    [nsmethod,bucketsize] = parseOptionalInputsCG(isrange,varargin{:});
end

end

function [nsmethod,bucketsize] = parseOptionalInputsCG(isrange,varargin)
% PARSEOPTIONALINPUTS  Parse optional PV pairs

coder.inline('always');
coder.internal.prefer_const(varargin);
params = struct( ...
    'Distance', uint32(0),...
    'P', uint32(0),...
    'Cov',uint32(0),...
    'Scale',uint32(0),...
    'NSMethod', uint32(0),...
    'BucketSize', uint32(0));

popts = struct( ...
    'CaseSensitivity', false, ...
    'StructExpand',    false, ...
    'PartialMatching', 'unique',...
    'IgnoreNulls',true);

if ~isrange
    params.IncludeTies = uint32(0);
    params.K = uint32(0);
    optarg = eml_parse_parameter_inputs(params, popts, varargin{:});
else
    optarg = eml_parse_parameter_inputs(params, popts, varargin{:});
end

nsmethod   = eml_get_parameter_value(...
    optarg.NSMethod,'', varargin{:});
bucketsize = eml_get_parameter_value(...
    optarg.BucketSize,[], varargin{:});

end


function [nsmethod,bucketsize] = parseOptionalInputsML(isrange,varargin)
% PARSEOPTIONALINPUTS  Parse optional PV pairs
parser = inputParser;
addParameter(parser,'NSMethod','');
addParameter(parser,'Distance','');
addParameter(parser,'P',[]);
addParameter(parser,'Cov',[]);
addParameter(parser,'Scale',[]);
addParameter(parser,'XSize',1);
addParameter(parser,'DistExtra','');
addParameter(parser,'DistParamExtra',[]);
addParameter(parser,'BucketSize',[]);

if ~isrange
    addParameter(parser,'K',1);
    addParameter(parser,'IncludeTies',false);
end
parse(parser,varargin{:});
nsmethod = parser.Results.NSMethod;
bucketsize = parser.Results.BucketSize;
end
