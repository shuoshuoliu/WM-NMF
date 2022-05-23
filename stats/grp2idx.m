function varargout = grp2idx(s)
% GRP2IDX  Create index vector from a grouping variable.
%   [G,GN] = GRP2IDX(S) creates an index vector G from the grouping
%   variable S. S can be a categorical, numeric, logical, datetime or 
%   duration vector; a cell vector of strings; or a character matrix with 
%   each row representing a group label. The result G is a vector taking 
%   integer values from 1 up to the number K of distinct groups. GN is a 
%   cell array of strings representing group labels. GN(G) reproduces S 
%   (aside from any differences in type).
%
%   Type "help groupingvariable" for more information about grouping
%   variables.
%
%   [G,GN,GL] = GRP2IDX(S) returns a column vector GL representing the
%   group levels. The set of groups and their order in GL and GN are the
%   same, except that GL has the same type as S. If S is a character
%   matrix, GL(G,:) reproduces S, otherwise GL(G) reproduces S.
%
%   GRP2IDX treats NaNs (numeric, duration or logical), empty strings (char
%   or cell array of strings), <undefined> values (categorical), or NaTs 
%   (datetime) in S as missing values and returns NaNs in the corresponding
%    rows of G. GN and GL don't include entries for missing values.
%
%   See also GROUPINGVARIABLE, GRPSTATS, GSCATTER.

%   Copyright 1999-2015 The MathWorks, Inc.

[varargout{1:nargout}] = statslib.internal.grp2idx(s);

end
