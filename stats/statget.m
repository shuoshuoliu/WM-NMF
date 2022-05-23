function val = statget(options,name,default,flag)
%STATGET Get STATS options parameter value.
%   VAL = STATGET(OPTIONS,'NAME') extracts the value of the named parameter
%   from statistics options structure OPTIONS, returning an empty matrix if
%   the parameter value is not specified in OPTIONS.  Case is ignored for
%   parameter names, and unique partial matches are allowed.
%   
%   VAL = STATGET(OPTIONS,'NAME',DEFAULT) extracts the named parameter, but
%   returns DEFAULT if the named parameter is not specified (is []) in
%   OPTIONS.
%   
%   See also STATSET.

%   Copyright 1993-2012 The MathWorks, Inc.


if nargin > 1
    name = convertStringsToChars(name);
end

if nargin > 2
    default = convertStringsToChars(default);
end

if nargin > 3
    flag = convertStringsToChars(flag);
end

if nargin < 2
    error(message('stats:statget:TooFewInputs'));
elseif nargin < 3
    default = [];
end

% Undocumented usage for fast access with no error checking.
if nargin == 4 && isequal('fast',flag)
    val = statgetfast(options,name,default);
    return
end

if ~isempty(options) && ~isa(options,'struct')
    error(message('stats:statget:InvalidOptions'));
end

if isempty(options)
    val = default;
    return;
end

names = {'Display'     'MaxFunEvals'  'MaxIter' ...
         'TolBnd'      'TolFun'       'TolTypeFun' ...
         'TolX'        'TolTypeX'     'GradObj' ...
         'Jacobian'    'DerivStep'    'FunValCheck' ...
         'Robust'      'RobustWgtFun' 'WgtFun' ...
         'Tune'        'UseParallel' 'UseSubstreams' ...
         'Streams'     'OutputFcn'};
lowNames = lower(names);
[~,i] = internal.stats.getParamVal(name,lowNames,'parameter name');
name = names{i};
val = options.(name);

% If empty, use default value, but keep the type of the original
% if the default is also empty.
if isempty(val) && nargin>=3  % only if a default was given explicitly
    val = default;
end


%------------------------------------------------------------------
function value = statgetfast(options,name,defaultopt)
%STATGETFAST Get STATS OPTIONS parameter with no error checking.
%   VAL = STATGETFAST(OPTIONS,FIELDNAME,DEFAULTOPTIONS) will get the value
%   of the FIELDNAME from OPTIONS with no error checking or fieldname
%   completion.  If the value is [], it gets the value of the FIELDNAME from
%   DEFAULTOPTIONS, another OPTIONS structure which is  probably a subset
%   of the options in OPTIONS.

if isempty(options)
    value = defaultopt.(name);
else
    value = options.(name);
    if isempty(value)
        value = defaultopt.(name);
    end
end
