function P = getRR2Perm(Base)
%GETRR2PERM Get HaltonRR2 scrambling permutation vector
%   GETRR2PERM(BASE) returns a permutation vector for the base specified.
%   This is calculated using cached copies of the permutations for bases
%   2^N.
    
%   Copyright 2007-2018 The MathWorks, Inc.


persistent CachedPerms CachedBases

if isempty(CachedPerms)
    % Always cache for primes up to 2^8.  Beyond this they are added as
    % needed
    CachedBases = 2.^(1:10);
    CachedPerms = internal.stats.haltonset.computeRR2Perm(CachedBases);
end

% Check that the cache goes high enough for this base
if Base>CachedBases(end)
    % Need to add to the cache
    NMax = length(CachedBases)+1;
    while Base>(2^NMax)
        NMax = NMax+1;
    end
    NewBases = 2.^((length(CachedBases)+1):NMax);
    CachedBases =[CachedBases, NewBases];
    CachedPerms = [CachedPerms, internal.stats.haltonset.computeRR2Perm(NewBases)];
end

P = internal.stats.haltonset.computeRR2Perm(Base, CachedPerms);
