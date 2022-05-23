function [err,d,c,f,wmsg]=checkselections(data,censoring,frequency,dval,cval,fval)

% For use by DFITTOOL


%   Copyright 2003-2011 The MathWorks, Inc.

err = '';
wmsg = '';
d_l = 0;
c_l = 0;
f_l = 0;
NONE = '(none)'; % saved in files, not to be translated
d = [];
c = [];
f = [];

dblwarn = false;  % warn about conversion to double

% If there's no data vector input, and not even a name or expression, then
% there's really no data, and it's an error.
if (nargin<4) && (isempty(data) || isequal(data, NONE) || all(isspace(data(:))))
    if isequal(data, NONE)
        err = getString(message('stats:dfstrings:sprintf_InvalidDataChoice', NONE));
    else
        err = getString(message('stats:dfstrings:sprintf_DataMustBeSpecified'));
    end
    
% Otherwise, there may be a data vector and no name/expression, and we give it
% a default name; or there may be a name/expression and no vector, and we will
% eval to get the data.
else
    if isempty(data)
        dataname = getString(message('stats:dfstrings:assignment_DataVariable'));
    else
        dataname = getString(message('stats:dfstrings:sprintf_DataVariable',data));
    end
    try
        if nargin<4
            dval = evalin('base',data);
        end
        if ~isa(dval,'double')
            d = double(dval);
            dblwarn = true;
        else
            d = dval;
        end
        if isvector(d) && (length(d) > 1)
            if any(isinf(d))
               err = getString(message('stats:dfstrings:sprintf_CannotContainInfOrInf', dataname));
            elseif ~isreal(d)
               err = getString(message('stats:dfstrings:sprintf_CannotBeComplex', dataname));
            elseif sum(~isnan(d))==0
               err = getString(message('stats:dfstrings:sprintf_ContainsAllNaNValues', dataname));
            else
               d_l = length(d);
            end
        else
            if ~isvector(d)
               err = getString(message('stats:dfstrings:sprintf_IsNotAVector', dataname));
            else
               err = getString(message('stats:dfstrings:sprintf_DoesNotContainAtLeast2Observations', dataname));
            end
        end
    catch ME
        err = [err getString(message('stats:dfstrings:sprintf_InvalidExpression', dataname, ME.message))];
    end
end

% If there's no censoring vector input, and not even a name or expression,
% then there's really no censoring.
if (nargin<5) && (isempty(censoring) || isequal(censoring, NONE) || all(isspace(censoring(:))))
    c_l = -1;
    
% Otherwise, there may be a censoring vector and no name/expression, and we
% give it a default name; or there may be a name/expression and no vector, and
% we will eval to get the vector.
else
    if isempty(censoring)
        censname = getString(message('stats:dfstrings:assignment_CensoringVariable'));
    else
        censname = getString(message('stats:dfstrings:sprintf_CensoringVariable',censoring));
    end
    try
        if nargin<5
            cval = evalin('base',censoring);
        end
        if ~isa(cval,'double')
            if ~isa(cval,'logical')
               dblwarn = true;
            end
            c = double(cval);
        else
            c = cval;
        end
        if isempty(c)
           c_l = -1;
        elseif isvector(c) && (length(c) > 1)
            if ~all(ismember(c, 0:1))
                err = [err getString(message('stats:dfstrings:sprintf_MustBeALogicalVector',censname))];
            elseif any(isinf(c))
                err = [err getString(message('stats:dfstrings:sprintf_CannotContainInfOrInf', censname))];
            elseif ~isreal(c)
                err = [err getString(message('stats:dfstrings:sprintf_CannotBeComplex', censname))];
            else
                c_l = length(c);
            end
        else
            err = [err getString(message('stats:dfstrings:sprintf_IsNotAVector', censname))];
        end
    catch ME
        err = [err getString(message('stats:dfstrings:sprintf_InvalidExpression',censname, ME.message))];
    end
end

% If there's no frequency vector input, and not even a name or expression,
% then there's really no frequencies.
if (nargin<6) && (isempty(frequency) || isequal(frequency, NONE) || all(isspace(frequency(:))))
    f_l = -1;

% Otherwise, there may be a frequency vector and no name/expression, and we
% give it a default name; or there may be a name/expression and no vector, and
% we will eval to get the vector.
else
    if isempty(frequency)
        freqname = getString(message('stats:dfstrings:assignment_FrequencyVariable'));
    else
        freqname = getString(message('stats:dfstrings:sprintf_FrequencyVariable',frequency));
    end
    try
        if nargin<6
            fval = evalin('base',frequency);
        end
        if ~isa(fval,'double')
            if ~isa(fval,'logical')
               dblwarn = true;
            end
            f = double(fval);
        else
            f = fval;
        end
        if isempty(f)
           f_l = -1;
        elseif isvector(f) && (length(f) > 1)
            if any(f<0) || any(f~=round(f) & ~isnan(f))
               err = [err getString(message('stats:dfstrings:sprintf_MustBeNonnegativeIntegers',freqname))];
            elseif any(isinf(f))
               err = [err getString(message('stats:dfstrings:sprintf_CannotContainInfOrInf', freqname))];
            elseif ~isreal(f)
               err = [err getString(message('stats:dfstrings:sprintf_CannotBeComplex', freqname))];
            else
               f_l = length(f);
            end
        else
            err = [err getString(message('stats:dfstrings:sprintf_IsNotAVector', freqname))];
        end
    catch ME
        err = [err getString(message('stats:dfstrings:sprintf_InvalidExpression', freqname, ME.message))];
    end
end

% Check lengths if no other errors
if isequal(err, '')
    if ((c_l ~= -1) && (c_l ~= d_l)) || ((f_l ~= -1) && (f_l ~= d_l))
        err = getString(message('stats:dfstrings:sprintf_LengthsMustBeEqual'));
        err = [err getString(message('stats:dfstrings:sprintf_DataLength', d_l))];
        if (c_l ~= -1) && (c_l ~= d_l)
            err = [err getString(message('stats:dfstrings:sprintf_CensoringLength', c_l))];
        end
        if (f_l ~= -1) && (f_l ~= d_l)
            err = [err getString(message('stats:dfstrings:sprintf_FrequencyLength', f_l))];
        end
    end
end

% Must have some non-censored data
if isempty(err) && c_l~=-1
   if (f_l==-1 && all(c==1)) || (f_l~=-1 && all(c(f>0)==1))
      err = getString(message('stats:dfstrings:assignment_AllCensored'));
   end
end

% Warn if we had to convert to double
if isempty(err) && dblwarn
   wmsg = getString(message('stats:dfstrings:sprintf_RequiresDouble'));
end
      
