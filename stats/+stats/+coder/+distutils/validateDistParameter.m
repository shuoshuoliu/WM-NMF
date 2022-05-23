function validateDistParameter(distance,dp,xsize)
%#codegen
%VALIDATENSPARAMS helper function to validate Distance Parameter for
%   knnsearch and rangesearch.

%   Copyright 2017 The MathWorks, Inc.

coder.inline('always');

% MATLAB path needed for ML object to codegen object redirection, which
% happens in MATLAB.
if coder.target('MATLAB')
    if strncmpi(distance,'min',3)
        validateattributes(dp,{'double','single'},{'scalar','nonempty','positive','nonnan','real'},mfilename,'DistanceParameter');
    elseif strncmpi(distance,'mah',3)
        validateattributes(dp,{'double','single'},{'nonempty','square','nonnan','finite','real'},mfilename,'DistanceParameter');
    else %seuclidean
        validateattributes(dp,{'double','single'},{'vector','nonempty','nonnegative','nonnan','finite','real'},mfilename,'DistanceParameter');
    end
else% not using validateattributes to be able to prevent compile-time checking for run-time variable distance & distance parameter arguments
    % NOOP since distance is enforced to be compile-time constant.
    %if coder.internal.isConst(dp)
    %    coder.internal.assert(coder.internal.isConst(distance),'stats:coder:pdist2:ExpectedConstantArg');
    %end
    validateattributes(dp,{'double','single'},{'nonempty'},mfilename,'DistanceParameter');
    if coder.internal.isConst(distance)
        if strncmpi(distance,'min',3)
            coder.internal.assert((isreal(dp) && isscalar(dp) && dp(1) > 0 && ~isnan(dp(1))), ...
                'stats:pdist2:BadMinExp'); % need a message id
        elseif strncmpi(distance,'mah',3)
            coder.internal.assert((ismatrix(dp) && size(dp,1) == size(dp,2) && size(dp,1) == xsize && allfinite(dp(:)) && allnonnan(dp(:))), ...
                'stats:pdist2:InvalidCov');
        else %seuclidean
            coder.internal.assert((isvector(dp) && allnonneg(dp) && numel(dp) == xsize  && allnonnan(dp(:))), ...
                'stats:pdist2:InvalidWeights');
        end
    else
        %validateattributes(dp,{'numeric'},{'nonempty','nonnan'},mfilename,'DistanceParameter');
        ISMIN = strncmpi(distance,'min',3);
        ISMAH = strncmpi(distance,'mah',3);
        coder.internal.assert(~ISMIN ||(isreal(dp) && isscalar(dp) && dp(1) > 0 && ~isnan(dp(1))), ...
            'stats:pdist2:BadMinExp'); % need a message id
        coder.internal.assert(~ISMAH || ...
            (ismatrix(dp) && size(dp,1) == size(dp,2) && size(dp,1) == xsize && allfinite(dp(:)) && allnonnan(dp(:))), ...
            'stats:pdist2:InvalidCov');
        coder.internal.assert(ISMIN || ISMAH || ...
            (isvector(dp) && allnonneg(dp) && numel(dp) == xsize && allnonnan(dp(:))), ...
            'stats:pdist2:InvalidWeights');
    end
end
end

function p = allnonneg(x)
% local function to validate if all elements are nonnegative
p = true;
for k = 1:numel(x)
    p = p && x(k) >= 0;
end
end

function p = allfinite(x)
% local function to validate if all elements are finite
p = true;
for k = 1:numel(x)
    p = p && isfinite(x(k));
end
end

function p = allnonnan(x)
% local function to validate if all elements are nonnan
p = true;
for k = 1:numel(x)
    p = p && ~isnan(x(k));
end
end