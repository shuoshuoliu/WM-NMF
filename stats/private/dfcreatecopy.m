function [new, fittype]=dfcreatecopy(original);


%   Copyright 2003-2004 The MathWorks, Inc.

fittype = original.fittype;

new = copyfit(original);
new = java(new);

