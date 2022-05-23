function [x, conf] = fracfact(gen,varargin)
%FRACFACT Fractional factorial design for two-level factors.
%   X = FRACFACT(GEN) produces the fractional factorial design defined by
%   the generator string GEN.  GEN must be a sequence of "words" separated
%   by spaces.  If the generator string consists of P words using K letters
%   of the alphabet, then X has N=2^K rows and P columns.  Each word
%   defines how the corresponding factor's levels are defined as products
%   of generators from a 2^K full-factorial design.  Alternatively, GEN can
%   be a string array or cell array of strings, with one word per cell.
%
%   [X, CONF] = FRACFACT(GEN) also returns CONF, a cell array of
%   strings containing the confounding pattern for the design.
%
%   [...] = FRACFACT(GEN, 'PARAM1',val1, 'PARAM2',val2,...) specifies one
%   or more of the following name/value pairs:
%
%       'MaxInt'      Maximum level of interaction to include in the
%                     confounding output (default 2)
%       'FactorNames' Cell array specifying the name for each factor
%                     (default names are X1, X2, ...)
%
%   Example:
%      x = fracfact('a b c abc')
%
%   produces an 8-run fractional factorial design for four variables, where
%   the first three columns are an 8-run 2-level full factorial design for
%   the first three variables, and the fourth column is the product of the
%   first three columns.  The fourth column is confounded with the
%   three-way interaction of the first three columns.
%
%   See also FF2N, FULLFACT, FRACFACTGEN.


%   Copyright 1993-2011 The MathWorks, Inc.


% Check the input argument for validity
if nargin > 0
    gen = convertStringsToChars(gen);
end

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if iscellstr(gen)
    gen = sprintf('%s ',gen{:});
end

if (~ischar(gen)) || (size(gen,1) > 1)
   error(message('stats:fracfact:BadGenerators'))
end
if (~all(isletter(gen) | isspace(gen) | (gen == '-')))
   error(message('stats:fracfact:BadGenerators'))
end

okargs =   {'maxint' 'factornames'};
defaults = {2        []};
[maxint, names] = internal.stats.parseArgs(okargs,defaults,varargin{:});

% Find out the dimensions of the design
genCell = textscan(strtrim(gen),'%s');
gen = char(genCell{1});
alpha = [('a':'z'),('A':'Z')];
if ~any(ismember(alpha(1:26),gen))
    gen = lower(gen);
end

alpha = alpha(ismember(alpha, gen));
K = length(alpha);
if (K == 0)
   error(message('stats:fracfact:BadGenerators'))
end
alpha = ['-' alpha];
N = 2^K;
P = size(gen, 1);
bitmat = (zeros(P, K+1) > 0);

% Generate a full factorial as the basis, augment with -1's
Full = [-ones(N,1) -1+2*ff2n(K)];
x = ones(N, P);

% Multiply columns as specified by generators; also create bit pattern
ident = eye(K+1);
for j=1:P
   g = deblank(gen(j,:));
   cols = loc(g, alpha);
   cols = sort(cols(cols>0));
   x(:,j) = prod(Full(:,cols), 2);
   bitmat(j,:) = mod(sum(ident(cols,:),1), 2);
end

if (nargout > 1)
   if ~isnumeric(maxint) || ~isreal(maxint) || ~isscalar(maxint) ...
                         || maxint<=0 || maxint~=floor(maxint)
       error(message('stats:fracfact:BadMaxInt'));
   end
   
   if isempty(names)
       names = strcat({'X'}, strjust(num2str((1:P)'), 'left'));
   elseif ~isvector(names) || ~iscellstr(names) || numel(names)~=P
       error(message('stats:fracfact:BadFactorNames'));
   end

   conf = confounding(bitmat, alpha, maxint, names);
end

function conf = confounding(bitmat, alpha, maxint, names)
% CONFOUNDING returns the confounding pattern for the design whose
% generators are encoded in the rows of the matrix bitmat and whose
% words are written in the alphabet alpha.  The first letter of the
% alphabet is '-' and is treated specially.  The display shows
% interactions with up to maxint factors and gets factor names
% from the array names.
nfact = size(bitmat, 1);
K = length(alpha);

% Write factor main effects into a cell array
conf = cell(nfact+1, 3);
conf(1,:) = {'Term', 'Generator', 'Confounding'};
for j=2:nfact+1
   conf{j,1} = names{j-1};
   mask = bitmat(j-1,:);
   if (sum(mask) > 0)
      gen = alpha(mask);
   else
      gen = '1';
   end
   conf{j,2} = gen;
end

% Write next higher order effects into the cell array
for j=2:min(maxint,K)
   
   % First allocate enough additional rows
   group = nchoosek(1:nfact, j);
   base = size(conf,1);
   ninter = size(group,1);
   conf{base+ninter, 1} = [];
   
   % Loop over all combinations
   for j1=1:ninter
      
      % Construct term name
      g = group(j1,:);
      name = sprintf('*%s',conf{1+g});
      name(1) = [];  % remove leading *
      
      % Construct term generator
      mask = (mod(sum(bitmat(g,:),1), 2) > 0);
      jbase = base + j1;
      conf{jbase, 1} = name;
      if (sum(mask) > 0)
         conf{jbase,2} = alpha(mask);
      else
         conf{jbase,2} = '1';
      end
   end
end

% Fill in confounding entries
for j=2:size(conf, 1)
   if (isempty(conf{j,3}))          % only if this row is not yet done
      gen = conf{j,2};              % look at this generator
      if (isempty(gen))
         other = '- -';             % fake that will not match anything
      elseif (gen(1) == '-')
         other = gen(2:end);
      else
         other = ['-', gen];
      end
      
      % Find other terms with +/- the same generator
      rows1 = find(strcmp(conf(:,2), {gen}));
      rows2 = find(strcmp(conf(:,2), {other}));
      rows = sort([rows1(:); rows2(:)]);
      
      % Construct a confounding string for this term
      t = conf{rows(1),1};
      if (~strcmp(conf{rows(1), 2}, gen))
         t = ['-', t];
      end
      for j1=2:length(rows)
         if (strcmp(conf{rows(j1), 2}, gen))
            t = [t, ' + ', conf{rows(j1),1}];
         else
            t = [t, ' - ', conf{rows(j1),1}];
         end
      end
      for j1=rows1'
         conf{j1,3} = t;
      end
   end
end

% -----------------------
function x=loc(src, targ)
%LOC locates the first occurrence of each element of src in targ, returns index
x = zeros(size(src));
t = targ(:);
for j=length(t):-1:1
   x(src==t(j)) = j;
end
