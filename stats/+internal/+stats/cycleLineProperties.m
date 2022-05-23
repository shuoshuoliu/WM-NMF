function varargout = cycleLineProperties(nG,varargin)
%Cycle through the line properties by expansion/shrinkage
%   varargout = CYCLELINEPROPERTIES(nG,varargin) expands/shrinks the values
%   in each line property given in varargin to the number of groups nG.
%   Used by plotting functions with a grouping variable such as
%   scatterhist, gscatter. Can be used on the following line properties:
%      Property           value
%      - 'LineWidth'      scalar or a numerical vector
%      - 'LineSytle'      string or cell string
%      - 'MarkerSize'     scalar or a numerical vector
%      - 'MarkerSymbol'   string or cell string
%      - 'Color'          string , cell string or a RGB color matrix. 
%                        
%   Example:
%   >> clr = 'rgb';
%   >> ls = {'-','-.',':'};
%   >> lw  = 2;
%   >> [clr, ls, lw] = internal.stats.cycleLineProperties(5,clr,ls,lw);
%   The results are:   
%   clr = 'rgbrg'
%   ls =  {'-'    '-.'    ':'    '-'    '-.'}
%   lw =  [2     2     2     2     2]

%   Copyright 2012-2015 The MathWorks, Inc.

[varargout{1:nargout}] = statslib.internal.cycleLineProperties(nG,varargin{:});

end
