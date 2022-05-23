classdef cvpartitionInMemoryImpl < internal.stats.cvpartitionImpl
   
methods
        function cv = cvpartitionInMemoryImpl(N,method,T,varargin)
            
            if nargin>0 && isstruct(N) && isfield(N,'BckCmpBackdoorConstructor')
                % only loadobj uses this backdoor constructor for backward
                % compatibility
                N = rmfield(N,'BckCmpBackdoorConstructor');
                fn = fieldnames(N);
                for ii=1:numel(fn)
                    cv.(fn{ii}) = N.(fn{ii});
                end
                return
            end
            
            hasstratify = false;
                
            if isempty(varargin)
                s = RandStream.getGlobalStream;
                stdargin = nargin;
            else
                if strncmpi(method,{'resubstitution'},length(method)) % resubstitution method should not combined with 'Stratify'
                    error(message('stats:cvpartition:UnsupportedCV')); % verification here to keep error message consistent
                end
                
                if strncmpi(method,{'leaveout'},length(method)) && T(1)~=1 %leaveout method should not combined with 'Stratify'
                    error(message('stats:cvpartition:UnsupportedLeaveout'));
                end
                
                if length(varargin)==1 % must be a RandStream object
                    stdargin = nargin-1;
                    s = varargin{1};
                    if ~isa(s,'RandStream')
                        error(message('stats:cvpartition:cvpartBadArg'));
                    end
                elseif length(varargin)==2 % should be 'Stratify' NVP
                    stra = varargin{1};
                    s = RandStream.getGlobalStream;
                    stdargin = nargin;
                    if ~strncmpi(stra,{'stratify'},length(stra))
                        error(message('stats:cvpartition:UnknownThirdInput'));
                    end
                    stratified = varargin{2};
                    hasstratify = true;
                    if ~islogical(stratified)
                        error(message('stats:cvpartition:InvalidStratifyType'));
                    end
                elseif length(varargin)==3 % should be 'Stratify' NVP plus RandStream object
                    stdargin = nargin-1;
                    tmp1 = varargin{1};
                    tmp2 = varargin{2};
                    tmp3 = varargin{3};
                    if strncmpi(tmp1,{'stratify'},length(tmp1)) && isa(tmp3,'RandStream')
                        stratified = varargin{2};
                        s = tmp3;
                    elseif strncmpi(tmp2,{'stratify'},length(tmp2)) && isa(tmp1,'RandStream')
                        stratified = varargin{3};
                        s = tmp1;
                    else
                        error(message('stats:cvpartition:UnknownThirdInput'));
                    end
                    hasstratify = true;
                    if ~islogical(stratified)
                        error(message('stats:cvpartition:InvalidStratifyType'));
                    end
                else % >3, too many inputs
                    error(message('stats:cvpartition:cvpartTooManyArg'));
                end
            end
            
            if stdargin < 2
                error(message('stats:cvpartition:TooFewInputs'));
            end

            if ischar(method) && size(method,1) == 1
                methodNames = {'kfold','holdout','leaveout','resubstitution'};
                j = find(strncmpi(method,methodNames,length(method)));
                if length(j) > 1
                    error(message('stats:cvpartition:AmbiguousMethod', method));
                elseif isempty(j)
                    error(message('stats:cvpartition:UnknownMethod'));
                end
            else
                error(message('stats:cvpartition:InvalidType'));
            end

            cv.Type = methodNames{j};
            
            switch cv.Type
                case 'kfold'
                    if hasstratify && isscalar(N) && stratified % bad combination
                        error(message('stats:cvpartition:ScalarStratify'));
                    elseif hasstratify && isscalar(N) && ~stratified % do whatever currently doing
                    elseif hasstratify && ~isscalar(N) && stratified % do whatever currently doing
                    elseif hasstratify && ~isscalar(N) && ~stratified % do non-stratified cvpartition
                        if ischar(N)
                            N = cellstr(N);
                        end
                        tmpN = ismissing(N).*0;
                        tmpN(ismissing(N)) = nan;
                        N = tmpN;
                    end
                case 'leaveout'
                    if hasstratify
                        error(message('stats:cvpartition:LeaveOutStratify'));
                    end
                case 'resubstitution'
                    if hasstratify
                        error(message('stats:cvpartition:ResubstitutionStratify'));
                    end
                case 'holdout'
                    if hasstratify && isscalar(N) && stratified
                        error(message('stats:cvpartition:ScalarStratify'));
                    elseif hasstratify && isscalar(N) && ~stratified % do whatever currently doing
                    elseif hasstratify && ~isscalar(N) && stratified % do whatever currently doing
                    elseif hasstratify && ~isscalar(N) && ~stratified % do non-stratified cvpartition
                        if ischar(N)
                            N = cellstr(N);
                        end
                        tmpN = ismissing(N).*0;
                        tmpN(ismissing(N)) = nan;
                        N = tmpN;
                    end
            end
                        
            if isscalar(N)
                if ~isnumeric(N) || N <= 1 || N ~= round(N) || ~isfinite(N)
                    error(message('stats:cvpartition:BadNX'));      
                end
                cv.N = N;
            else
                cv.Group = grp2idx(N);
                cv.N = length(cv.Group); % the number of observations including NaNs
                [~,wasnan,cv.Group] = internal.stats.removenan(cv.Group);
                hadNaNs = any(wasnan);
                if hadNaNs
                    warning(message('stats:cvpartition:MissingGroupsRemoved'));
                    if length (cv.Group) <= 1
                        error(message('stats:cvpartition:BadNGrp'));
                    end
                end
            end

            dftK = 10; % the default number of subsamples(folds) for Kfold
            P  = 1/10; % the default holdout ratio

            switch cv.Type
                case 'kfold'                    
                    if stdargin == 2 || isempty(T)
                        T = dftK;
                    elseif ~isscalar(T) || ~isnumeric(T) || T <= 1 || ...
                            T ~= round(T) || ~isfinite(T)
                        error(message('stats:cvpartition:BadK'));
                    end

                    if  isempty(cv.Group) && T > cv.N
                        warning(message('stats:cvpartition:KfoldGTN'));
                        T = cv.N;
                    elseif ~isempty(cv.Group) && T > length(cv.Group)
                        warning(message('stats:cvpartition:KfoldGTGN'));
                        T = length(cv.Group);
                    end

                    cv.NumTestSets = T; %set the number of fold
                    cv = cv.rerandom(s);

                case 'leaveout'
                    if stdargin > 2 && T(1) ~= 1 %fix bug when compare char array with 1
                        error(message('stats:cvpartition:UnsupportedLeaveout'));
                    end

                    if isempty(cv.Group)
                        cv.NumTestSets = cv.N;
                    else
                        cv.NumTestSets = length(cv.Group);
                    end

                    [~,cv.indices] = sort(rand(s,cv.NumTestSets,1));

                    cv.TrainSize = (cv.NumTestSets-1) * ones(1,cv.NumTestSets);
                    cv.TestSize = ones(1,cv.NumTestSets);

                case 'resubstitution'
                    if stdargin > 2 
                        error(message('stats:cvpartition:UnsupportedCV'));
                    end

                    if isempty(cv.Group)
                        numObs = N;
                    else
                        numObs = length(cv.Group);
                    end

                    cv.indices = (1: numObs)';
                    cv.NumTestSets = 1;
                    cv.TrainSize =  numObs;
                    cv.TestSize =  numObs;

                case 'holdout'
                    if stdargin == 2 || isempty(T)
                        T = P;
                    elseif ~isscalar(T) || ~ isnumeric(T) || T <= 0 || ~isfinite(T)
                        error(message('stats:cvpartition:BadP'));
                    end

                    if T >= 1 %hold-T observations out
                        if T ~=round(T)
                            error(message('stats:cvpartition:BadP'));
                        end
                        if isempty(cv.Group)
                            if T >= cv.N
                                error(message('stats:cvpartition:PNotLTN'));
                            end
                        else
                            if T>= length(cv.Group)
                                error(message('stats:cvpartition:PNotLTGN'));
                            end
                        end
                    else
                        if (isempty(cv.Group) && floor(cv.N *T) == 0) ||...
                                (~isempty(cv.Group) && floor(length(cv.Group) * T) == 0)
                            error(message('stats:cvpartition:PTooSmall'));

                        end
                    end

                    cv.holdoutT = T;
                    cv = cv.rerandom(s);
                    cv.NumTestSets = 1;
            end

            %add NaNs back
            if ~isempty(cv.Group) && hadNaNs
                [cv.indices, cv.Group] =...
                    internal.stats.insertnan(wasnan, cv.indices, cv.Group);
            end
        end % cvpartition constructor    
    
end
    
end