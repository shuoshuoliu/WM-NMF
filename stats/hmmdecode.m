function [pStates,pSeq, fs, bs, s] = hmmdecode(seq,tr,e,varargin)
%HMMDECODE calculates the posterior state probabilities of a sequence.
%   PSTATES = HMMDECODE(SEQ,TRANSITIONS,EMISSIONS) calculates the posterior
%   state probabilities, PSTATES, of sequence SEQ from a Hidden Markov
%   Model specified by transition probability matrix,  TRANSITIONS, and
%   EMISSION probability matrix, EMISSIONS. TRANSITIONS(I,J) is the
%   probability of transition from state I to state J. EMISSIONS(K,SYM) is
%   the probability that symbol SYM is emitted from state K. The posterior
%   probability of sequence SEQ is the probability P(state at step i = k |
%   SEQ).  PSTATES is an array with the same length as SEQ and one row for
%   each state in the model. The (i,j) element of PSTATES gives the
%   probability that the model was in state i at the jth step of SEQ.
%
%   [PSTATES, LOGPSEQ] = HMMDECODE(SEQ,TR,E) returns, LOGPSEQ, the log of
%   the probability of sequence SEQ given transition matrix, TR and
%   emission matrix, E.
%
%   [PSTATES, LOGPSEQ, FORWARD, BACKWARD, S] = HMMDECODE(SEQ,TR,E) returns
%   the forward and backward probabilities of the sequence scaled by S. The
%   actual forward probabilities can be recovered by using:
%        f = FORWARD.*repmat(cumprod(s),size(FORWARD,1),1);
%   The actual backward probabilities can be recovered by using:
%       bscale = cumprod(S, 'reverse'); b =
%       BACKWARD.*repmat([bscale(2:end), 1],size(BACKWARD,1),1);
%
%   HMMDECODE(...,'SYMBOLS',SYMBOLS) allows you to specify the symbols that
%   are emitted. SYMBOLS can be a numeric array or a string array cell
%   array of the names of the symbols.  The default symbols are integers 1
%   through M, where N is the number of possible emissions.
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
%       [seq, states] = hmmgenerate(100,tr,e); 
%       pStates = hmmdecode(seq,tr,e);
%
%       [seq, states] = hmmgenerate(100,tr,e,'Symbols',...
%                 {'one','two','three','four','five','six'});
%       pStates = hmmdecode(seq,tr,e,'Symbols',...
%                 {'one','two','three','four','five','six'});
%
%   See also HMMGENERATE, HMMESTIMATE, HMMVITERBI, HMMTRAIN.

%   Reference: Biological Sequence Analysis, Durbin, Eddy, Krogh, and
%   Mitchison, Cambridge University Press, 1998.

%   Copyright 1993-2020 The MathWorks, Inc.


% tr must be square

if nargin > 0
    if isstring(seq)
        seq = cellstr(seq);
    end
end

if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

numStates = size(tr,1);
checkTr = size(tr,2);
if checkTr ~= numStates
    error(message('stats:hmmdecode:BadTransitions'));
end

% number of rows of e must be same as number of states

checkE  = size(e,1);
if checkE ~= numStates
    error(message('stats:hmmdecode:InputSizeMismatch'));
end

numSymbols = size(e,2);

% deal with options
if nargin > 3
    okargs = {'symbols'};
    symbols = internal.stats.parseArgs(okargs, {''}, varargin{:});
    
    if ~isempty(symbols)
        numSymbolNames = numel(symbols);
        if ~isvector(symbols) || numSymbolNames ~= numSymbols
            error(message('stats:hmmdecode:BadSymbols'));
        end
        [~, seq]  = ismember(seq,symbols);
        if any(seq(:)==0)
            error(message('stats:hmmdecode:MissingSymbol'));
        end
    end
end

if ~isnumeric(seq)
    error(message('stats:hmmdecode:MissingSymbolArg'));
end
numEmissions = size(e,2);
if ~(isempty(seq) || isvector(seq)) ...
        || any(seq(:)<1) ...
        || any(seq(:)~=round(seq(:))) ...
        || any(seq(:)>numEmissions) 
     error(message('stats:hmmdecode:BadSequence', numEmissions));
end

% add extra symbols to start to make algorithm cleaner at f0 and b0
if isrow(seq)
    seq = [numSymbols+1, seq];
else
    seq = [numSymbols+1; seq];
end
L = length(seq);

% This is what we'd like to do but it is numerically unstable
% warnState = warning('off');
% logTR = log(tr);
% logE = log(e);
% warning(warnState);
% f = zeros(numStates,L);
% f(1,1) = 1;
% % for count = 2:L
%     for state = 1:numStates
%         f(state,count) = logE(state,seq(count)) + log(sum( exp(f(:,count-1) + logTR(:,state))));
%     end
% end
% f = exp(f);

% so we introduce a scaling factor
fs = zeros(numStates,L);
fs(1,1) = 1;  % assume that we start in state 1.
s = zeros(1,L);
s(1) = 1;
for count = 2:L
    for state = 1:numStates
        fs(state,count) = e(state,seq(count)) .* (sum(fs(:,count-1) .*tr(:,state)));
    end
    % scale factor normalizes sum(fs,count) to be 1. 
    s(count) =  sum(fs(:,count));
    fs(:,count) =  fs(:,count)./s(count);
end

%  The  actual forward and  probabilities can be recovered by using
%   f = fs.*repmat(cumprod(s),size(fs,1),1);


% This is what we'd like to do but it is numerically unstable
% b = zeros(numStates,L);
% for count = L-1:-1:1
%     for state = 1:numStates
%         b(state,count) = log(sum(exp(logTR(state,:)' + logE(:,seq(count+1)) + b(:,count+1)  )));
%     end
% end

% so once again use the scale factor
bs = ones(numStates,L);
for count = L-1:-1:1
    for state = 1:numStates
      bs(state,count) = (1/s(count+1)) * sum( tr(state,:)'.* bs(:,count+1) .* e(:,seq(count+1))); 
    end
end

%  The  actual backward and  probabilities can be recovered by using
%  scales = cumprod(s, 'reverse'); 
%  b = bs.*repmat([scales(2:end), 1],size(bs,1),1);

pSeq = sum(log(s));
pStates = fs.*bs;

% get rid of the column that we stuck in to deal with the f0 and b0 
pStates(:,1) = [];


