function [table, stats] = gagerr(y,group,varargin)
% GAGERR Gage repeatability and reproducibility (R&R) study
%     GAGERR(Y,{PART,OPERATOR}) performs a gage R&R study on measurements
%     in Y collected by OPERATOR on PART.  Y is a column vector containing
%     the measurements on different parts.  PART and OPERATOR are
%     categorical variables, numeric vectors, string array, character
%     matrices, or cell arrays of strings.  The number of elements in PART
%     and OPERATOR should be the same as in Y.
%
%     A table is printed in the command window in which the decomposition 
%     of variance, standard deviation, study var (5.15*standard deviation)
%     are listed with respective percentages for different sources.
%     Summary statistics are printed below the table giving the number of 
%     distinct categories (NDC) and the percentage of Gage R&R of total 
%     variations (PRR).
%
%     A bar graph is also plotted showing the percentage of different
%     components of variations. Gage R&R, repeatability, reproducibility,
%     and part to part variations are plotted as four vertical bars.
%     Variance and study var are plotted as two groups.
%
%     The guideline to determine the capability of a measurement system
%     using NDC is the following:
%           (1) If NDC > 5, the measurement system is capable
%           (2) If NDC < 2, the measurement system is not capable
%           (3) Otherwise, the measurement system may be acceptable
%
%     The guideline to determine the capability of a measurement
%     system using PRR is the following:
%           (1) If PRR < 10%, the measurement system is capable
%           (2) If PRR > 30%, the measurement system is not capable
%           (3) Otherwise, the measurement system may be acceptable
% 
%     GAGERR(Y,GROUP) performs a gage R&R study on measurements in Y
%     with PART and OPERATOR represented in GROUP. GROUP is a numeric 
%     matrix whose first and second columns specify different parts and
%     operators respectively. The number of rows in GROUP should be the
%     same as the number of elements in Y.
% 
%     GAGERR(Y,PART) performs a gage R&R study  on measurements in Y
%     without operator information. The assumption is that all variability 
%     is contributed by PART. 
% 
%     GAGERR(...,'PARAM1',val1,'PARAM2',val2,...) performs a gage R&R study
%     using  one or more of the following parameter name/value pairs:
%
%       Parameter       Value
%
%       'spec'          A two element vector which defines the lower and 
%                       upper limit of the process, respectively. In this 
%                       case, summary statistics printed in the command 
%                       window include Precision-to-Tolerance Ratio (PTR). 
%                       Also, the bar graph includes an additional group, 
%                       the percentage of tolerance.
%
%                       The guideline to determine the capability of a
%                       measurement system using PTR is the following: 
%                         (1) If PTR < 0.1, the measurement system is
%                         capable
%                         (2) If PTR > 0.3, the measurement system is not
%                         capable
%                         (3) Otherwise, the measurement system may be
%                         acceptable
%
%      'printtable'     A string with a value 'on' or 'off' which indicates
%                       whether the tabular output should be printed in the
%                       command window or not. The default value is 'on'.
%
%      'printgraph'     A string with a value 'on' or 'off' which indicates
%                       whether the bar graph should be plotted or not. The 
%                       default value is 'on'.
%
%      'randomoperator' A logical value, true or false, which indicates
%                       whether the effect of OPERATOR is random or not. 
%
%      'model'          The model to use, specified by one of:
%                         'linear' -- Main effects only (default)
%                         'interaction' -- Main effects plus two-factor 
%                                          interactions
%                         'nested' -- Nest OPERATOR in PART
%     
%    [TABLE, STATS] = GAGERR(...) returns a 6x5 matrix TABLE and a 
%    structure STATS. The columns of TABLE, from left to right, represent 
%    variance, percentage of variance, standard deviations, study var, and 
%    percentage of study var. The rows of TABLE, from top to bottom, 
%    represent different sources of variations: gage R&R, repeatability, 
%    reproducibility, operator, operator and part interactions, and part. 
%    STATS is a structure containing summary statistics for the performance 
%    of the measurement system. The fields of STATS are:
%          ndc -- Number of distinct categories
%          prr -- Percentage of gage R&R of total variations
%          ptr -- Precision-to-tolerance ratio. The value is NaN if the 
%                 parameter 'spec' is not given. 
%
%  Example 
%    Conduct a gage R&R study for a simulated measurement system using a 
%    mixed ANOVA model without interactions:
%       y = randn(100,1);                                % measurements
%       part = ceil(3*rand(100,1));                      % parts
%       operator = ceil(4*rand(100,1));                  % operators
%       gagerr(y,{part, operator},'randomoperator',true) % analysis
%

% Copyright 2006-2011 The MathWorks, Inc.
if nargin > 1
    if isstring(group)
        group = cellstr(group);
    end
end

if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

if nargin<2
    error(message('stats:gagerr:FewInput'));
end

if isa(group,'categorical')  % may be a single categorical variable
    group = {group};
elseif isnumeric(group)      % may be a matrix with one or two columns
    group = num2cell(group,1);
end

if length(group)>2
    error(message('stats:gagerr:BadGroup'));
end;
if cellfun(@length,group) ~= size(y,1);
    error(message('stats:gagerr:MismatchGroup'));
end
 
% Measurements must be a numeric vector
if ~isvector(y) || ~isnumeric(y)
    error(message('stats:gagerr:BadY'));
end
 
% Determine whether the operator factor is given
nooperators = (length(group)==1);
 
 
args =   {'model', 'randomoperator','spec', 'printtable','printgraph'};
defaults = {'linear',true, [],'on', 'on'};
[model,randomoperator,spec,printtable,printgraph] ...
                   =  internal.stats.parseArgs(args,defaults,varargin{:});
 
if ~ischar(model) || size(model,1)>1 
    error(message('stats:gagerr:BadModelString')) 
end 
model = lower(model);
if ~any(strcmp(model,{'linear','interaction','nested'}))
    error(message('stats:gagerr:BadModel')) 
end;
 
if randomoperator ~= 1 && randomoperator ~= 0
    error(message('stats:gagerr:BadRandomoperator')) 
end 
 
switch model 
    case 'linear'
        if randomoperator
            modelnum = 1;
        else
            modelnum = 4;
        end;
    case 'interaction'
        if randomoperator
            modelnum = 2;
        else
            modelnum = 5;
        end;
    case 'nested'
         modelnum = 3;
end;
 
if  ~isnumeric(spec) || (isnumeric(spec)&& numel(spec)~=2 && numel(spec)~=0)
    error(message('stats:gagerr:BadSpec')) 
end
 
% check argument printtable
if ~ischar(printtable) || size(printtable,1)>1
    error(message('stats:gagerr:BadPrinttableString')) 
end
printtable = lower(printtable);
outputtable = strcmp(printtable, 'on');
if ~outputtable && ~strcmp(printtable, 'off');
    error(message('stats:gagerr:BadPrinttable'))
end;
    
% check argument printgraph
if ~ischar(printgraph)|| size(printgraph,1)>1
    error(message('stats:gagerr:BadPrintgraphString')) 
end
printgraph = lower(printgraph);
outputgraph = strcmp(printgraph, 'on');
if ~outputgraph && ~strcmp(printgraph, 'off');
    error(message('stats:gagerr:BadPrintgraph'))
end;
 
if nooperators 
    [p,atab,anovastats] = anovan(y,group,'random',1,'display','off'); 
else
    group = group(:);
    switch modelnum
        case 1 % cross model without interaction
            [p,atab,anovastats] = anovan(y,group, 'random',[1 2],'display','off');
        case 2 % cross model with interaction
            [p,atab,anovastats] = anovan(y,group, ...
                'model','interaction','random',[1 2],'display','off');
        case 3 % nested model 
            [p,atab,anovastats] = anovan(y,group, ...
                'nested',[0 0; 1 0],'random',[1 2],'display','off');
        case 4 % mixed model without interaction
            [p,atab,anovastats] = anovan(y,group, 'random',1,'display','off');     
        case 5 % mixed model with interaction
            [p,atab,anovastats] = anovan(y,group, ...
                'model','interaction','random',1,'display','off');     
    end;          
end
 
% Variance breakdown
repeatabilityvar = anovastats.varest(end);
parttopartvar = max(0, anovastats.varest(1));
 
if nooperators || modelnum>=4            % operator is not random
    operatorvar = 0;
else 
    operatorvar = max(0, anovastats.varest(2));
end; 
 
if any(modelnum ==[2,5])
    interactionvar = max(0, anovastats.varest(end-1));
else
    interactionvar = 0;
end;
 
% variation components 
allvar = [repeatabilityvar operatorvar interactionvar parttopartvar];
%      Gage R&R         |Repeatability | Reproducibility | Operator |Operator*Part |Part to Part
vars = [sum(allvar(1:3)) allvar(1)    sum(allvar(2:3))   allvar(2:end)];    % variance components
totalvar = sum(allvar);                 % total variance 
varspercentage = 100* vars/totalvar;    % percentage of variance
stds =  sqrt(vars);                     % standard deviation components;
totalstd =sqrt(totalvar);               % total standard deviation  
studyvar = 5.15*stds;                    % study var components
stdspercentage = 100* stds/totalstd;    % percentage of standard deviations
 
% summary statistics 
ndc = round(sqrt(2*vars(6)/vars(1)));           % NDC = 1.41*Part to Part/gage RR  
prr = stdspercentage(1);                     % PRR = gage RR / total
% PTR can only be calculated if we have spec.
if ~isempty(spec) 
    ptr = 5.15*stds(1)/abs(spec(2)-spec(1)); % PTR = 5.15*gage RR /diff(spec)
else 
    ptr = NaN;
end;
 
% print results in command window
if outputtable
    colheaders = {getString(message('stats:gagerr:ColHeadSource')),      ...
                  getString(message('stats:gagerr:ColHeadVariance')),    ...
                  getString(message('stats:gagerr:ColHeadPctVariance')), ...
                  getString(message('stats:gagerr:ColHeadSigma')),       ...
                  getString(message('stats:gagerr:ColHead515Sigma')),    ...
                  getString(message('stats:gagerr:ColHeadPct515Sigma'))};
    rowheaders = {getString(message('stats:gagerr:XLabelGageRR')),                 ...
                  ['  ' getString(message('stats:gagerr:XLabelRepeatability'))],   ...
                  ['  ' getString(message('stats:gagerr:XLabelReproducibility'))], ...
                  ['   ' getString(message('stats:gagerr:RowHeadOperator'))],      ...
                  ['   ' getString(message('stats:gagerr:RowHeadPartOperator'))],  ...
                  getString(message('stats:gagerr:RowHeadPart'))}';
    varscol = num2cell(vars');
    varspercentagecol = num2cell(varspercentage');
    stdscol = num2cell(stds');
    studyvarcol = num2cell(studyvar');
    stdspercentagecol = num2cell(stdspercentage');

    tblcore = horzcat(rowheaders, varscol, varspercentagecol, stdscol, studyvarcol, stdspercentagecol);

    % Remove rows in they do not show up
    partop_row = 5;
    if nooperators || modelnum >= 4% operator item does not show up 
        tblcore(4, :) = [];
        partop_row = partop_row - 1;
    end
    if all(modelnum ~= [2, 5]) % interaction item only shows up if the model includes that item 
        tblcore(partop_row, :) = [];
    end
    
    totalrow = {getString(message('stats:gagerr:RowHeadTotal')), totalvar, 100, totalstd, 5.15*totalstd, '' };

    tbl = vertcat(colheaders, tblcore, totalrow);
    disp(tbl)
    disp([getString(message('stats:gagerr:TableNoteNDC')), num2str(ndc)]);  
    disp([getString(message('stats:gagerr:TableNoteTotalVar')), sprintf('%6.2f', prr)]);  
    if ~isempty(spec)    % PTR is printed if we have spec.
          disp([getString(message('stats:gagerr:TableNotePTR')), sprintf('%6.2f', ptr)]);
    end;
    disp(getString(message('stats:gagerr:TableNoteLastCol')))
end;
 
if outputgraph 
    % bar plot of variation decomposition
    idx = [1,2,3,6];  % indices for Gage R&R, Repeatability, Reproducibility, and Part to Part
    % PTR is plotted only if we have spec.
    if isempty(spec)        
        h = bar([varspercentage(idx);stdspercentage(idx)]','grouped');
        legend(getString(message('stats:gagerr:LegendVariance')),getString(message('stats:gagerr:LegendStudyVar')))
    else
        tolerancepercentage = 100 * 5.15*stds/abs(spec(2)-spec(1));
        h = bar([varspercentage(idx);stdspercentage(idx);tolerancepercentage(idx)]','grouped');
        legend(getString(message('stats:gagerr:LegendVariance')),getString(message('stats:gagerr:LegendStudyVar')),getString(message('stats:gagerr:LegendTolerance')))
        set(h(3),'Tag','tol');
    end;
    set(gca,'XTickLabel',{getString(message('stats:gagerr:XLabelGageRR')), getString(message('stats:gagerr:XLabelRepeatability')),getString(message('stats:gagerr:XLabelReproducibility')),getString(message('stats:gagerr:XLabelPartToPart'))},'XTick',1:4)    
    ylabel(getString(message('stats:gagerr:LabelPercent')))
    set(h(1),'Tag','var');
    set(h(2),'Tag','std');
end;
 
if nargout>0
    table = [vars; varspercentage; stds; studyvar;stdspercentage]';    
end;
 
if nargout>1
    stats.ndc  = ndc; 
    stats.prr  = prr; 
    stats.ptr  = ptr;
end;
