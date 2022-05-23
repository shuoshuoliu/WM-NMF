function [seq,  states]= hmmgenerate(L,tr,e,varargin)
%HMMGENERATE generate a sequence for a Hidden Markov Model
%   [SEQ, STATES] = HMMGENERATE(LEN,TRANSITIONS,EMISSIONS) generates
%   sequence of emission symbols, SEQ, and a random sequence of states,
%   STATES, of length LEN from a Markov Model specified by transition
%   probability matrix, TRANSITIONS, and EMISSION probability matrix,
%   EMISSIONS. TRANSITIONS(I,J) is the probability of transition from state
%   I to state J. EMISSIONS(K,L) is the probability that symbol L is
%   emitted from state K.
%
%   HMMGENERATE(...,'SYMBOLS',SYMBOLS) allows you to specify the symbols
%   that are emitted. SYMBOLS can be a numeric array or a string/cell array
%   of the names of the symbols. The default symbols are integers 1 through
%   N, where N is the number of possible emissions.
%
%   HMMGENERATE(...,'STATENAMES',STATENAMES) allows you to specify the
%   names of the states. STATENAMES can be a numeric array or a cell array
%   of the names of the states. The default statenames are 1 through M,
%   where M is the number of states.
%
%   This function always starts the model in state 1 and then makes a
%   transition to the first step using the probabilities in the first row
%   of the transition matrix. So in the example given below, the first
%   element of the output states will be 1 with probability 0.95 and 2 with
%   probability .05.
%
%   Examples:
%
%       tr = [0.95,0.05; ...
%             0.10,0.90];
%           
%       e = [1/6,  1/6,  1/6,  1/6,  1/6,  1/6; ...
%            1/10, 1/10, 1/10, 1/10, 1/10, 1/2;];
%
%       [seq, states] = hmmgenerate(100,tr,e)
%
%       [seq, states] = hmmgenerate(100,tr,e,'Symbols',...
%                 {'one','two','three','four','five','six'},...
%                  'Statenames',{'fair';'loaded'})
%
%   See also  HMMVITERBI, HMMDECODE, HMMESTIMATE, HMMTRAIN.

%   Copyright 1993-2011 The MathWorks, Inc.


if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

seq = zeros(1,L);
states = zeros(1,L);

% tr must be square

numStates = size(tr,1);
checkTr = size(tr,2);
if checkTr ~= numStates
    error(message('stats:hmmgenerate:BadTransitions'));
end

% number of rows of e must be same as number of states

checkE = size(e,1);
if checkE ~= numStates
    error(message('stats:hmmgenerate:InputSizeMismatch'));
end

numEmissions = size(e,2);

customSymbols = false;
customStatenames = false;

if nargin > 3
    okargs = {'symbols','statenames'};
    [symbols,statenames] = ...
        internal.stats.parseArgs(okargs, {[] []}, varargin{:});
    
    if ~isempty(symbols)
        numSymbolNames = numel(symbols);
        if ~isvector(symbols) || numSymbolNames ~= numEmissions
            error(message('stats:hmmgenerate:BadSymbols'));
        end
        customSymbols = true;
    end
    if ~isempty(statenames)
        numStateNames = length(statenames);
        if numStateNames ~= numStates
            error(message('stats:hmmgenerate:BadStateNames'));
        end
        customStatenames = true;
    end 
end

% create two random sequences, one for state changes, one for emission
statechange = rand(1,L);
randvals = rand(1,L);

% calculate cumulative probabilities
trc = cumsum(tr,2);
ec = cumsum(e,2);

% normalize these just in case they don't sum to 1.
trc = trc./repmat(trc(:,end),1,numStates);
ec = ec./repmat(ec(:,end),1,numEmissions);

% Assume that we start in state 1.
currentstate = 1;

% main loop 
for count = 1:L
    % calculate state transition
    stateVal = statechange(count);
    state = 1;
    for innerState = numStates-1:-1:1
        if stateVal > trc(currentstate,innerState)
            state = innerState + 1;
            break;
        end
    end
    % calculate emission
    val = randvals(count);
    emit = 1;
    for inner = numEmissions-1:-1:1
        if val  > ec(state,inner)
            emit = inner + 1;
            break
        end
    end
    % add values and states to output
    seq(count) = emit;
    states(count) = state;
    currentstate = state;
end

% deal with names/symbols
if customSymbols
    seq = reshape(symbols(seq),1,L);
end
if customStatenames
    states = reshape(statenames(states),1,L);
end

