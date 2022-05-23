function [y,N,classname] = getGLMBinomialData(distr,y,N)
%getGLMBinomialData Convert data to standard form, separate binomial y/N.

%   Copyright 2016-2020 The MathWorks, Inc.

classname = '';
if nargin<3 || isempty(N)
    N = ones(size(y,1),1);
end

if strcmpi(distr, 'binomial')
    % Categorical responses allowed for 'binomial'
    if isa(y,'categorical')
        [y, textname, classname] = grp2idx(y);
        nc = length(textname);
        if nc > 2
            error(message('stats:glmfit:TwoLevelCategory'));
        end
        y(y==1) = 0;
        y(y==2) = 1;
    end
    
    if ~islogical(y)
        internal.stats.checkSupportedNumeric('Y',y,false,false,false,true); %gpu okay
    end
    try
        % Separate y and N
        wasvector = true;
        if ~(isnumeric(y) || islogical(y))
            error(message('stats:glmfit:BadDataBinomialFormat'))
        elseif size(y,2) == 2
            wasvector = false;
            N = y(:,2);
            y = y(:,1);
        elseif size(y,2) ~= 1
            error(message('stats:glmfit:MatrixOrBernoulliRequired'));
        end
        if any(N<0) || any(~isnan(N) & N~=round(N))
            if wasvector
                error(message('stats:GeneralizedLinearModel:BadNValues'))
            else
                error(message('stats:GeneralizedLinearModel:BadBinomialResponse'))
            end
        elseif any(y < 0 | y > N)
            error(message('stats:glmfit:BadDataBinomialRange'));
        end
        y = y./N;
    catch ME
        throwAsCaller(ME)
    end
else
    if ~islogical(y)
        internal.stats.checkSupportedNumeric('Y',y,false,false,false,true); %gpu okay
    end

end
end
