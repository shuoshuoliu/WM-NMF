function [includeties,k,distance,p,cov,scale,bucketsize] = extractNSParams(varargin)
%#codegen
% VALIDATENSPARAMS helper function to parse and validate inputs for
% knnsearch and rangesearch.
coder.inline('always');
coder.internal.prefer_const(varargin);
% MATLAB path needed for ML object to codegen object redirection, which
% happens in MATLAB.
if coder.target('MATLAB') 
     [includeties,k,distance,p,cov,scale,bucketsize] = parseOptionalInputsML(varargin{:});
else
     [includeties,k,distance,p,cov,scale,bucketsize] = parseOptionalInputsCG(varargin{:});
end

end

function [includeties,k,distance,p,cov,scale,bucketsize] = parseOptionalInputsCG(varargin)
% PARSEOPTIONALINPUTS  Parse optional PV pairs

coder.inline('always');
coder.internal.prefer_const(varargin);
params = struct( ...
    'IncludeTies', uint32(0), ...
    'K',     uint32(0),...
    'Distance',     uint32(0),...
    'P', uint32(0),...
    'Cov',uint32(0),...
    'Scale',uint32(0),...
    'NSMethod', uint32(0),...
    'BucketSize',uint32(0));

% don't extract BucketSize until KDTreeSearcher is supported
%'BucketSize',uint32(0)

popts = struct( ...
    'CaseSensitivity', false, ...
    'StructExpand',    false, ...
    'PartialMatching', 'unique',...
    'IgnoreNulls',true);

optarg           = eml_parse_parameter_inputs(params, popts, ...
    varargin{:});
includeties   = eml_get_parameter_value(...
    optarg.IncludeTies,false, varargin{:});
k         = eml_get_parameter_value(...
    optarg.K, coder.internal.indexInt(1), varargin{:});
distance = eml_get_parameter_value(...
    optarg.Distance, '', varargin{:});
p = eml_get_parameter_value(...
    optarg.P,[], varargin{:});
cov = eml_get_parameter_value(...
    optarg.Cov,[], varargin{:});
scale = eml_get_parameter_value(...
    optarg.Scale,[], varargin{:});

bucketsize = eml_get_parameter_value(...
    optarg.BucketSize,[], varargin{:});

end


function [includeties,k,distance,p,cov,scale,bucketsize] = parseOptionalInputsML(varargin)
% PARSEOPTIONALINPUTS  Parse optional PV pairs

parser = inputParser;
addParameter(parser,'Distance','');
addParameter(parser,'K',1);
addParameter(parser,'P',[]);
addParameter(parser,'Cov',[]);
addParameter(parser,'Scale',[]);
addParameter(parser,'IncludeTies',false);
addParameter(parser,'NSMethod','');
addParameter(parser,'BucketSize',[]);
parse(parser,varargin{:});
k = parser.Results.K;
includeties = parser.Results.IncludeTies;
distance = parser.Results.Distance;
p = parser.Results.P;
cov = parser.Results.Cov;
scale = parser.Results.Scale;
bucketsize = parser.Results.BucketSize;
end
