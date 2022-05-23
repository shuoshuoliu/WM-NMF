function varargout = linsolve(A, B, varargin)
%LINSOLVE Solve linear system A*X=B with gpuArray support
%   This is a helper function to work around linsolve not supporting
%   gpuArray inputs.
%
%   For gpuArray inputs, mldivide and rank will be used.  Otherwise the
%   built-in linsolve will be used.

%   Copyright 2020 The MathWorks, Inc.

if isa(A, 'gpuArray') || isa(B, 'gpuArray')
    varargout{1} = A\B;
    if nargout == 2
        varargout{2} = iCalculateSecondOutput(A, varargin{:});
    end
else
    [varargout{1:nargout}] = linsolve(A,B,varargin{:});
end
end

function r = iCalculateSecondOutput(A, varargin)
if nargin > 1
    % If opts is specified, then r is the reciprocal of the condition
    % number of A unless RECT is true and both LT and UT are false, in
    % which case, r gives the rank of A.
    opts = varargin{1};
    LT = isfield(opts, "LT") && opts.LT;
    UT = isfield(opts, "UT") && opts.UT;
    RECT = isfield(opts, "RECT") && opts.RECT;
    if RECT && ~LT && ~UT
        r = rank(A);
    else
        r = cond(A,1);
    end
    return
end

if ~diff(size(X))
    % Square matrix so r is the reciprocal condition number of A.
    r = cond(A,1);
else
    % Rectangular matrix so r is the rank of A.
    r = rank(A);
end
end