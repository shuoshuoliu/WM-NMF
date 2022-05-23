classdef(Sealed = true) qrandstream < handle
%QRANDSTREAM Create a quasi-random stream.
%   Q = QRANDSTREAM(TYPE,D) creates a stream reference object that
%   encapsulates the specified type of point set. TYPE is a string
%   containing the name of a point set and must be one of 'sobol' or
%   'halton'.  D specifies the number of dimensions.
%
%   Q = QRANDSTREAM(TYPE,D,PROP,VAL,...) specifies a set of property-value
%   pairs that are used in creating the point set on which the stream is
%   based.
%
%   Q = QRANDSTREAM(PS) constructs a stream based on a copy of the point
%   set PS.  PS must be a quasi-random point set object, for example
%   sobolset or haltonset.
%
%   A quasi-random stream has the following properties and methods:
%
%   qrandstream properties:
%     Read-only:
%      PointSet     - The point set object that the stream draws from.
%
%     Settable:
%      State        - The index in the point stream of the last point.
%
%   qrandstream methods:
%      qrand,rand   - Draw points from the stream and increment its state.
%      reset        - Set the stream back to its initial state.
%
%   Examples:
%
%      Create a 5-dimensional stream based on a Sobol sequence and draw
%      values from it:
%         Q = qrandstream('sobol',5);
%         qrand(Q,10000)
%
%      Create a stream based on a leaping Halton sequence that skips the
%      initial point:
%         Q = qrandstream('halton',5,'Leap',12,'Skip',1);
%         qrand(Q,100)
%
%      Create a stream based on a scrambled sobol set:
%         S = sobolset(5);
%         S = scramble(S,'MatousekAffineOwen');
%         Q = qrandstream(S);
%         qrand(Q,2^16)
%
%      Create a stream and draw values from it using rand:
%         Q = qrandstream('sobol',5);
%         rand(Q, 100, 5)
%
%   See also QRAND, HALTONSET, SOBOLSET, RAND, RANDSTREAM.

%   Copyright 2007-2017 The MathWorks, Inc.



properties(Dependent = true)
    %STATE Current state of the stream.
    %   The State property of a quasi-random stream contains the index into
    %   the associated point set of the next point that will be drawn in
    %   the stream.  Getting and resetting the State property allows you
    %   to return a stream to a previous state.  The initial value of State
    %   is 1.
    %
    %   Example:
    %      Q = qrandstream('sobol', 5);
    %      s = Q.State;
    %      u1 = qrand(Q, 10)
    %      Q.State = s;
    %      u2 = qrand(Q, 10)   % contains exactly the same values as u1
    %
    %   See also QRAND.
    State;
end

properties(SetAccess='private')
    %POINTSET (Read-only) Point set from which the stream is drawn.
    %   The PointSet property contains a copy of the point set from which
    %   the stream is providing points.  The point set is specified during
    %   construction of a quasi-random stream and cannot subsequently be
    %   altered.
    %
    %   Example: 
    %      Q = qrandstream('sobol', 5, 'Skip', 8);
    %      % Create a new stream based on the same sequence as that in Q
    %      Q2 = qrandstream(Q.PointSet);
    %      u1 = qrand(Q, 10)
    %      u2 = qrand(Q2, 10)   % contains exactly the same values as u1 
    %
    %   See also QRANDSTREAM.
    PointSet;
end

properties(SetAccess='private', GetAccess='private')
    StateObject;
end


% State handling methods
methods
    function set.State(obj, val)
    if (val>length(obj.PointSet))
        error(message('stats:qrandstream:InvalidState'));
    end
    try
        obj.StateObject.resetState(val);
    catch E
        error(message('stats:qrandstream:ResetStateFailure', E.message));
    end
    end
    function val = get.State(obj)
    val = obj.StateObject.Index;
    end
    function reset(obj)
    %RESET Reset the stream to the default state.
    %   RESET(Q) sets the stream Q back to its initial state. Subsequent
    %   points will be the same as those produced by a new stream.
    %
    %   The state of Q can also be accessed and set by using its State
    %   property.
    %
    %   Example:
    %
    %      Q = qrandstream('sobol',5);
    %      X = qrand(Q,5)
    %      X = qrand(Q,5)
    %      reset(Q);
    %      X = qrand(Q,5)
    %
    %   See also QRANDSTREAM, STATE.

    obj.State = 1;
    end
end


% Constructor
methods
    function obj = qrandstream(Type, varargin)
    %QRANDSTREAM Construct a new quasi-random stream.
    %   Q = QRANDSTREAM(TYPE,D) creates a stream reference object that
    %   encapsulates the specified type of point set. TYPE is a string
    %   containing the name of a point set and must be one of 'sobol' or
    %   'halton'.  D specifies the number of dimensions.
    %
    %   Q = QRANDSTREAM(TYPE,D,PROP,VAL,...) specifies a set of
    %   property-value pairs that are applied to the point set before
    %   creating the stream.
    %
    %   Q = QRANDSTREAM(PS) constructs a stream based on the point set PS.
    %   PS must be an object of class qrandset, for example sobolset or
    %   haltonset.
    %
    %   See also QRAND, QRANDSET, RAND, RANDSTREAM.

    if nargin==0
        error(message('stats:qrandstream:MissingPointSetType'));
    end
    Type = convertStringsToChars(Type);
    if iscell(varargin)
        [varargin{:}] = convertStringsToChars(varargin{:});
    end
    
    if ischar(Type)
        switch lower(Type)
            case 'sobol'
                qps = sobolset(varargin{:});
            case 'halton'
                qps = haltonset(varargin{:});
            otherwise
                error(message('stats:qrandstream:InvalidCharPointSetType'));
        end
    elseif isa(Type, 'qrandset')
        if nargin==1
            qps = Type;
        else
            error(message('stats:qrandstream:InvalidNumberOfArguments'));
        end
    else
        error(message('stats:qrandstream:InvalidPointSetType'));
    end
    obj.PointSet = qps;
    obj.StateObject = createStreamState(qps);
    end
end

% Point generation methods
methods
    function vals = qrand(obj, N)
    %QRAND Generate quasi-random points from stream.
    %   QRAND(Q,N) returns an N-by-D matrix containing successive
    %   quasi-random points from the point set that the quasi-random stream
    %   Q contains.  Each row of the matrix contains a multi-dimensional
    %   quasi-random point in D dimensions.
    %
    %   QRAND(Q) returns a 1-by-D vector containing a single point.
    %
    %   Examples:
    %
    %      Generate successive points from a Sobol sequence:
    %         Q = qrandstream('sobol',5);
    %         X = qrand(Q)
    %         X = qrand(Q)
    %
    %      Generate points in blocks of 256:
    %         Q = qrandstream('sobol',5);
    %         X = qrand(Q,256)
    %         X = qrand(Q,256)
    %
    %   See also QRANDSTREAM, RAND, RANDSTREAM/RAND.

    if nargin<2
        N = 1;
    elseif ~checkPosScalar(N)
        error(message('stats:qrandstream:InvalidNumberOfPoints'));
    end

    vals = getStreamPoints(obj.PointSet, obj.StateObject, double(N));
    end

    
    function vals = rand(obj, varargin)
    %RAND Generate quasi-random points from stream.
    %   RAND returns a matrix of quasi-random values and is intended to
    %   allow qrandstream objects to be used in code that contains calls to
    %   the RAND method of the pseudo-random randstream class. Due to the
    %   multi-dimensional nature of quasi-random numbers, only some
    %   syntaxes of RAND are supported by qrandstream, as detailed below.
    %
    %   RAND(Q,N) returns a N-by-N matrix only when N is equal to the
    %   number of dimensions.  Any other value of N produces an error.
    %
    %   RAND(Q) returns a scalar only when the stream is in one dimension.
    %   Having more than one dimension in Q produces an error.
    %
    %   RAND(Q,M,N) or RAND(Q,[M,N]) returns an M-by-N matrix only when N
    %   is equal to the number of dimensions in the stream. Any other value
    %   of N produces an error.
    %
    %   RAND(Q,M,N,P,...) or RAND(Q,[M,N,P,...]) produces an error unless P
    %   and all following dimensions sizes are equal to one.
    %
    %   Example:
    %
    %      Generate the first 256 points from the Sobol sequence:
    %         Q = qrandstream('sobol',5);
    %         X = rand(Q,256,5)
    %
    %   See also QRANDSTREAM, QRAND, RAND, RANDSTREAM/RAND.

    if nargin==1
        Dims = [1 1];
    else
        % Check for a trailing char argument.  This would be a class
        % specifier in rand but is not supported here.
        if ischar(varargin{end})
            error(message('stats:qrandstream:OutputTypeNotSupported'));
        end
        
        % Decided whether the user has specified a single vector of
        % dimensions or separate inputs.
        if length(varargin)==1
            Dims = varargin{1};
            if ~checkPositiveVector(Dims)
                error(message('stats:qrandstream:InvalidDimensions'));
            end
        else
            % Check that all dim arguments are valid scalars
            if ~all(cellfun(@checkPosScalar, varargin))
                error(message('stats:qrandstream:DimensionsNotPositiveScalarInts'));
            end
            
            Dims = double([varargin{:}]);
        end

        if isscalar(Dims)
            Dims = [Dims Dims];
        end
    end
    if length(Dims)>2 && prod(Dims(3:end))>1
        error(message('stats:qrandstream:TooManyDimensions'));
    elseif Dims(2)~=obj.PointSet.Dimensions
        error(message('stats:qrandstream:DimensionsMismatch'));
    else
        vals = qrand(obj, Dims(1));
    end
    end
end


% Display  and other information methods
methods
    function disp(obj)
    %DISP Display a qrandstream object.
    %   DISP(Q) displays the quasi-random stream Q, without printing the
    %   variable name.  DISP prints the type and number of dimensions in
    %   the stream, and follows it with the list of point set properties.
    %
    %   See also QRANDSTREAM.

    if strcmp(get(0,'FormatSpacing'),'loose')
        LooseLine = '\n';
    else
        LooseLine = '';
    end

    fprintf('%s\n',getString(message('stats:qrandstream:DispQuasirandomStream', ...
        obj.PointSet.Type, obj.PointSet.Dimensions)));
    fprintf(LooseLine);
    fprintf('%s:\n',getString(message('stats:qrandstream:DispPointSetProperties')));
    dispProperties(obj.PointSet);
    fprintf(LooseLine);
    end
end


% Override methods to prevent () subscripting and object arrays.
methods(Hidden)
    function varargout = subsref(obj, S)
    %SUBSREF Subscripted reference for qrandstream.
    %   X = SUBSREF(Q,S) is called for the syntax Q(I), Q{I}, or Q.I.  S is
    %   a structure array with the fields:
    %      type -- string containing '()', '{}', or '.' specifying the
    %              subscript type.
    %      subs -- Cell array or string containing the actual subscripts.
    %
    %   See also QRANDSTREAM, SUBSREF.

    if strcmp(S(1).type, '.')
        Identifier = S(1).subs;
        if isPublic(obj, Identifier, 'get') || isPublic(obj, Identifier, 'call')
            % Property access or method call
            [varargout{1:nargout}] = builtin('subsref', obj, S);
        else
            error(message('stats:qrandstream:NoSuchMethodOrField', Identifier));
        end
    elseif strcmp(S(1).type, '()')
        error(message('stats:qrandstream:IndexingNotAllowed'));
    else  % {} index
        error(message('stats:qrandstream:CellRefFromNonCell'));
    end
    end

    function varargout = subsasgn(obj, S, Value)
    %SUBSASGN Subscripted assignment for qrandstream.
    %   Q = SUBSASGN(Q,S,X) is called for the syntax Q(I)=X, Q{I}=X, or
    %   Q.I=X.  S is a structure array with the fields:
    %      type -- string containing '()', '{}', or '.' specifying the
    %              subscript type.
    %      subs -- Cell array or string containing the actual subscripts.
    %
    %   See also QRANDSTREAM, SUBSASGN.
    
    if strcmp(S(1).type, '.')
        Identifier = S(1).subs;
        if isPublic(obj, Identifier, 'set')
            % Set a property
            [varargout{1:nargout}] = builtin('subsasgn', obj, S, Value);

        elseif isPublic(obj, Identifier, 'get')
            % Property is visible so we don't want to pretend
            % it doesn't exist
            error(message('stats:qrandstream:SetProhibited', Identifier));
        else
            error(message('stats:qrandstream:NoSuchField', Identifier));
        end

    elseif strcmp(S(1).type, '()')
        error(message('stats:qrandstream:AssignmentNotAllowed'));
    else % {} assignment
        error(message('stats:qrandstream:CellAssToNonCell'));
    end
    end

    % Override and hide methods for concatenation and general array
    % manipulation.
    function varargout = horzcat(varargin), throwCatError; end
    function varargout = vertcat(varargin), throwCatError; end
    function varargout = cat(varargin),     throwCatError; end
    function varargout = repmat(varargin),  throwCatError; end
    function varargout = transpose(varargin)
        error(message('stats:qrandstream:TransposeNotAllowed'));
    end
    function varargout = ctranspose(varargin)
        error(message('stats:qrandstream:TransposeNotAllowed'));
    end
    function varargout = reshape(varargin)
        error(message('stats:qrandstream:ReshapeNotAllowed'));
    end
    function varargout = permute(varargin)
        error(message('stats:qrandstream:PermuteNotAllowed'));
    end
end
methods(Hidden,Static)
    function obj = empty(varargin)
        error(message('stats:qrandstream:EmptyNotAllowed'));
    end
end
end

function OK = isPublic(obj, Ident, Type)
% Checks whether an identifier is a public property or method.

persistent PublicIdents
if isempty(PublicIdents)
    % Initialise the lists of public identifiers
    C = metaclass(obj);
    P = [C.Properties{:}];
    PropGet = strcmp({P.GetAccess}, 'public');
    PropSet = strcmp({P.SetAccess}, 'public');
    PropNames = {P.Name};

    M = [C.Methods{:}];
    MethAccess = strcmp({M.Access}, 'public');
    MethNames = {M.Name};

    PublicIdents.get = PropNames(PropGet);
    PublicIdents.set = PropNames(PropSet);
    PublicIdents.call = MethNames(MethAccess);
end

OK = any(strcmp(Ident, PublicIdents.(Type)));
end


function throwCatError
m = message('stats:qrandstream:NoCatAllowed');
E = MException(m.Identifier,'%s',getString(m));
E.throwAsCaller;
end


function ok = checkPosScalar(val)
ok = isnumeric(val) && isscalar(val) ...
    && isreal(val) && isfinite(val) && val>0 && val==fix(val);
end

function ok = checkPositiveVector(val)
ok = isnumeric(val) ...
    && size(val,1)==1 && ndims(val)==2 ...
    && isreal(val) && all(isfinite(val)) && all(val>0) && all(val==fix(val));
end
