function r = gammaincratio(x,s)
% GAMMAINCRATIO Ratio of incomplete gamma function values at S and S-1.
%    R=GAMMAINCRATIO(X,S) computes GAMMAINC(X,S)/GAMMAINC(X,S-1). S
%    must be greater or equal to 1. X and S must be of the same size.

%   Copyright 2012 The MathWorks, Inc.

% Initialize
r = zeros(size(s));

% What is small S?
smalls = s<2 | s<=x;

% For small S, use the ratio computed directly
if any(smalls(:))
    r(smalls) = gammainc(x(smalls),s(smalls))./gammainc(x(smalls),s(smalls)-1);
end

% For large S, estimate numerator and denominator using 'scaledlower' option
if any(~smalls(:))
    idx = find(~smalls);
    x = x(idx);
    s = s(idx);
    r(idx) = gammainc(x,s,'scaledlower')./gammainc(x,s-1,'scaledlower').*x./s;
end
end
