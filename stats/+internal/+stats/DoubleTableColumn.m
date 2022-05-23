classdef DoubleTableColumn < double
% DoubleTableColumn
%     The DoubleTableColumn class represents a column vector of doubles
%     to be used in a table display in a dataset. It provides a method
%     to indicate that some values should be considered "absent" and not
%     displayed in the table. An example might be an anova table, where
%     there is a column of F statistics but the "Error" row should not
%     have an F statistic.
%
%     See also DoubleTableColumn.DoubleTableColumn.
    

%   Copyright 2011-2020 The MathWorks, Inc.

    properties
        absent
    end
    
    methods(Access='public')
        function d = DoubleTableColumn(data,absent)
% DoubleTableColumn.DoubleTableColumn Constructor
%     C = DoubleTableColumn(VALS,ABSENT) accepts a column vector VAL
%     of double values and a logical vector ABSENT of the same size
%     as V. The result C displays blanks in places where ABSENT has
%     the value TRUE.
%
%     The output C is suitable for use as a variable in a dataset
%     when some variable values should not be present.

            if nargin==0
                data = [];
            end
            
            if ~isfloat(data) || ~iscolumn(data)
                error(message('stats:internal:DoubleTableColumn:BadData'));
            end
            if nargin<2 || isempty(absent)
                absent = false(size(data));
            elseif ~isequal(size(data),size(absent)) || ~islogical(absent)
                error(message('stats:internal:DoubleTableColumn:BadAbsent'));
            end
            
            % Protect against passing gpuArrays into the superclass
            % constructor. See g2163082.
            if isa(data, 'gpuArray')
                data = gather(data);
            end
            
            d = d@double(data);
            d.absent = absent;
        end
        
        function disp(d)
            % Display numbers in the usual format, but blank out entries
            % marked as absent
            dblFmt = getFloatFormats();
            varChars = num2str(double(d),dblFmt);
            varChars(d.absent,:) = ' ';
            disp(varChars)
        end
        
        function a = num2str(d,varargin)
            a = num2str(double(d),varargin{:});
            a(d.absent,:) = ' ';
        end
        
        function sref = subsref(obj,s)
            switch s(1).type
                case '.'
                    switch s(1).subs
                        case 'data'
                            sref = double(obj);
                        case 'absent'
                            sref = obj.absent;
                            if length(s)>1 && strcmp(s(2).type, '()')
                                sref = subsref(sref,s(2:end));
                            end
                        otherwise
                            error(message('stats:internal:DoubleTableColumn:BadField'));
                    end
                case '()'
                    d = double(obj);
                    if ~isscalar(s)
                        error(message('stats:internal:DoubleTableColumn:BadParenIndex'));
                    end
                    sref = subsref(d,s);
% *** it may be better to return a double
%                     sref = DoubleTableColumn(subsref(d,s),...
%                                     subsref(obj.absent,s));
                otherwise
                    error(message('stats:internal:DoubleTableColumn:BadCellIndex'));
            end
        end
        function idx = subsidx(s)
            idx = double(s);
        end
        
        function newobj = horzcat(varargin)
            c1 = cellfun(@(x)double(x),varargin,'UniformOutput',false);
            newobj = horzcat(c1{:});  % result is a double
        end
        function newobj = vertcat(varargin)
            c1 = cellfun(@(x)double(x),varargin,'UniformOutput',false);
            d = vertcat(c1{:});
            c2 = cellfun(@(x)x.absent,varargin,'UniformOutput',false);
            a = vertcat(c2{:});
            newobj = internal.stats.DoubleTableColumn(d,a);
        end
        function newobj = cat(dim,varargin)
            if isequal(dim,1)
                newobj = vertcat(varargin{:});
            else
                error(message('stats:internal:DoubleTableColumn:BadDim'));
            end
        end
    end
end

function [dblFmt,snglFmt] = getFloatFormats()
% Display for double/single will follow 'format long/short g/e' or 'format bank'
% from the command window. 'format long/short' (no 'g/e') is not supported
% because it often needs to print a leading scale factor.
switch lower(matlab.internal.display.format)
    case {'short' 'shortg' 'shorteng'}
        dblFmt  = '%.5g    ';
        snglFmt = '%.5g    ';
    case {'long' 'longg' 'longeng'}
        dblFmt  = '%.15g    ';
        snglFmt = '%.7g    ';
    case 'shorte'
        dblFmt  = '%.4e    ';
        snglFmt = '%.4e    ';
    case 'longe'
        dblFmt  = '%.14e    ';
        snglFmt = '%.6e    ';
    case 'bank'
        dblFmt  = '%.2f    ';
        snglFmt = '%.2f    ';
    otherwise % rat, hex, + fall back to shortg
        dblFmt  = '%.5g    ';
        snglFmt = '%.5g    ';
end
end
