classdef ExhaustiveSearcher < NeighborSearcher
%ExhaustiveSearcher Neighbor search object using exhaustive search.
%   An ExhaustiveSearcher object performs KNN (K-nearest-neighbor) search
%   or radius search using exhaustive search. You can create an
%   ExhaustiveSearcher object based on X using either of the following
%   syntaxes:
%
%   CREATENS function:
%        NS = CREATENS(X,'NSMethod','exhaustive')
%   ExhaustiveSearcher constructor:
%        NS = ExhaustiveSearcher(X)
%
%   Rows of X correspond to observations and columns correspond to
%   variables. When one of the above methods creates an ExhaustiveSearcher
%   object, it saves X. The KNNSEARCH method and RANGESEARCH method compute
%   the distance values from all the points in X to the query points to
%   find requested points.
%
%   ExhaustiveSearcher properties:
%       X               - Data used to create the object.
%       Distance        - The distance metric.
%       DistParameter   - The additional parameter for the distance metric.
%
%   ExhaustiveSearcher methods:
%       ExhaustiveSearcher - Construct an ExhaustiveSearcher object
%       knnsearch          - Find nearest neighbors for query points.
%       rangesearch        - Find points within a given radius for query
%                            points.
%
%   Example:
%      % Create an ExhaustiveSearcher object for data X with the cosine 
%      % distance. The kd-tree method does not support nearest neighbors 
%      % search for the cosine distance, therefore CREATENS creates an
%      % ExhaustiveSearcher object.
%      X = randn(100,5);
%      Y = randn(25,5);
%      ns = createns(X, 'distance', 'cosine');
%
%      % Find 5 nearest neighbors in X for each point in Y.
%      idx = knnsearch(ns,Y,'k',5);
%
%      % Find the points in X whose distance are not greater than 0.3 to
%      % the points in Y.
%      idx = rangesearch(ns,Y,0.3);
%
%      % Choose a different distance metric from the one used when creating
%      % this ExhaustiveSearcher object.
%      idx2 = knnsearch(ns,Y,'k',5, 'distance','correlation');
%
%  See also  CREATENS, KDTreeSearcher, STATS/KNNSEARCH, STATS/RANGESEARCH.

%   Copyright 2009-2017 The MathWorks, Inc.


  properties(GetAccess=protected, SetAccess=private) 
      checkNegativeDistance = false;
   end
   

  methods (Access = protected)
        function this = setDistance(this, distMetric)
            distMetric = convertStringsToChars(distMetric);
            if ischar(distMetric) 
                methods = { 'minkowski';'seuclidean'; 'mahalanobis'; ...
                     'euclidean';'cityblock'; 'chebychev'; ...
                     'cosine'; 'correlation'; ...
                    'spearman'; 'hamming'; 'jaccard'};
                i =  find(strncmpi(distMetric, methods, length(distMetric)));
                if length(i) > 1
                    error(message('stats:ExhaustiveSearcher:AmbiguousDistance', distMetric));
                elseif isempty(i)
                    error(message('stats:ExhaustiveSearcher:UnrecognizedDistance',distMetric));
                end
                this.PrivDistance = methods{i};
                %Set default distance parameter if needed
                if i == 1
                    this.PrivDistParameter = 2;
                elseif i == 2
                    this.PrivDistParameter = nanstd(this.X,[],1);
                elseif i == 3
                    this.PrivDistParameter = nancov(this.X);
                else 
                    this.PrivDistParameter = [];
                end
            elseif isa(distMetric,'function_handle')
                this.PrivDistance = distMetric;
                this.PrivDistParameter = [];
            else
                error(message('stats:ExhaustiveSearcher:BadDistance'));
            end
 
        end
        
        function this = setDistParameter(this, para)
            i = find(strncmpi(this.Distance, {'minkowski','mahalanobis','seuclidean'},3));
            if isempty(i)
                 error(message('stats:ExhaustiveSearcher:InvalidDistanceParam'));
            elseif i == 1 %'min'
                if ~(isscalar(para) && isnumeric(para) && para > 0 )
                    error(message('stats:ExhaustiveSearcher:BadMinExp'))
                end
                              
            elseif i==2  %mah
                 nDims = size(this.X,2);
                if ~isequal(size(para),[nDims,nDims])
                    error(message('stats:ExhaustiveSearcher:BadCov'));
                end
                %use cholcov because we also need to check whether the matrix is symmetric
                [~,flag] = cholcov(para,0);
                if flag ~= 0
                    error(message('stats:ExhaustiveSearcher:SingularCov'));
                end
              
            else %secu
                nDims = size(this.X,2);
                if ~(isvector(para) && length(para) == nDims...
                        && all(para >= 0))
                    error(message('stats:ExhaustiveSearcher:BadScale'));
                end
           
            end
            this.PrivDistParameter= para;
        end
    end

methods
    function obj = ExhaustiveSearcher(X,varargin)
%ExhaustiveSearcher Construct an ExhaustiveSearcher object.
%   NS = ExhaustiveSearcher(X,'NAME1',VALUE1,...,'NAMEN',VALUEN) creates
%   an ExhaustiveSearcher object, specifying the following optional
%   argument name/value pairs:
%
%     Name         Value
%     'Distance'   A string or a function handle specifying the default
%                  distance metric when you call the KNNSEARCH method.
%                  The value can be one of the following:
%                  'euclidean'   - Euclidean distance (default).
%                  'seuclidean'  - Standardized Euclidean distance. Each
%                                  coordinate difference between X and a
%                                  query point is scaled by dividing by a
%                                  scale value S. The default value
%                                  of S is the standard deviation computed
%                                  from X, S=NANSTD(X). To specify another
%                                  value for S, use the 'Scale' argument.
%                 'cityblock'   -  City Block distance.
%                 'chebychev'   -  Chebychev distance (maximum coordinate
%                                  difference).
%                 'minkowski'   -  Minkowski distance. The default exponent
%                                  is 2. To specify a different exponent,
%                                  use the 'P' argument.
%                 'mahalanobis' -  Mahalanobis distance, computed using a
%                                  positive definite covariance matrix C.
%                                  The default value of C is the sample
%                                  covariance matrix, as computed by
%                                  NANCOV(X). To specify another value for
%                                  C, use the 'COV' argument.
%                 'cosine'      -  One minus the cosine of the included
%                                  angle between observations (treated as
%                                  vectors).
%                 'correlation' -  One minus the sample linear
%                                  correlation between observations
%                                  (treated as sequences of values).
%                 'spearman'    -  One minus the sample Spearman's rank
%                                  correlation between observations
%                                  (treated as sequences of values).
%                 'hamming'     -  Hamming distance, percentage of
%                                  coordinates that differ.
%                 'jaccard'     -  One minus the Jaccard coefficient, the
%                                  percentage of nonzero coordinates that
%                                  differ.
%                  distance     -  A distance function specified using @
%                                  (for example @DISTFUN). A distance
%                                  function must be of the form:
%
%                                  function D2 = DISTFUN(ZI, ZJ),
%
%                                  taking as arguments a 1-by-N vector ZI
%                                  containing a single row from X or from
%                                  the query points Y, and an M2-by-N
%                                  matrix ZJ containing multiple rows of X
%                                  or Y, and returning an M2-by-1 vector of
%                                  distances D2, whose Jth element is the
%                                  distance between the observations ZI and
%                                  ZJ(J,:).
%
%    'P'          A positive scalar indicating the exponent for Minkowski
%                 distance. This argument is only valid when 'Distance' is
%                 'minkowski'. Default is 2.
%
%    'Cov'        A positive definite matrix indicating the covariance
%                 matrix when computing the Mahalanobis distance. This
%                 argument is only valid when 'Distance' is 'mahalanobis'.
%                 Default is NANCOV(X).
%
%    'Scale'      A vector S containing non-negative values, with length
%                 equal to the number of columns in X. Each coordinate
%                 difference between X and a query point is scaled by the
%                 corresponding element of S when computing the
%                 standardized Euclidean distance. This argument is only
%                 valid when 'Distance' is 'seuclidean'. Default is
%                 NANSTD(X).
%
%   See also ExhaustiveSearcher, CREATENS.

        if nargin == 0
            error(message('stats:ExhaustiveSearcher:NoDataInput'));
        end
        pnames = {'distance'    'p'  'cov' 'scale' 'checknegativedistance'};
        dflts =  { 'euclidean'  []    []    []      false};
        
        [varargin{:}] = convertStringsToChars(varargin{:});
        [dist,minExp, mahaCov, seucInvWgt,checkNegDist] = ...
             internal.stats.parseArgs(pnames, dflts, varargin{:});
         
           
        if isempty(dist)
            dist = 'euclidean' ;
            obj.PrivDistance = dist;
        else
            obj.Distance = dist;
        end

        if ischar(dist) %  %built-in distance
            %Integer/logical/char/anything data will be converted to double.
            %Complex floating point data can't be handled by a built-in
            %distance function.
            
            if ~isfloat(X)
                warning(message('stats:ExhaustiveSearcher:DataConversion', class( X )));
                X = double(X);
            end
            if ~isreal(X)
                error(message('stats:ExhaustiveSearcher:ComplexData'));
            end
        end
        obj.X = X;   
        
        if strncmp(obj.Distance,'minkowski',3)
            if  ~isempty(minExp)
                obj.DistParameter = minExp;
            else
                obj.PrivDistParameter = 2;
            end
        elseif ~isempty(minExp)
            error(message('stats:ExhaustiveSearcher:InvalidMinExp'));
        end
        
        [nx, nDims]= size(X);
        if strncmp(obj.Distance,'seuclidean',3)
            if ~isempty(seucInvWgt)
                obj.DistParameter = seucInvWgt;
            else
                obj.PrivDistParameter = nanstd(X,[],1);
            end
            
        elseif ~isempty(seucInvWgt)
            error(message('stats:ExhaustiveSearcher:InvalidScale'));
        end
        
        if strncmp(obj.Distance, 'mahalanobis',3)
            if ~isempty(mahaCov)
                if ~isequal(size(mahaCov),[nDims,nDims])
                    error(message('stats:ExhaustiveSearcher:BadCov'));
                end
                %use cholcov because we also need to check whether the matrix is symmetric
                [~,flag] = cholcov(mahaCov,0);
                if flag ~= 0
                    error(message('stats:ExhaustiveSearcher:SingularCov'));
                end
                obj.DistParameter = mahaCov;
            else
                if nx == 1
                    error(message('stats:ExhaustiveSearcher:TooFewXRowsForMah'));
                    
                end
                obj.DistParameter = nancov(X);
            end
        elseif ~isempty(mahaCov)
            error(message('stats:ExhaustiveSearcher:InvalidCov'));
        end
        
        if ~islogical(checkNegDist) || ~isscalar(checkNegDist)
            error(message('stats:ExhaustiveSearcher:BadCheckDist'));
        else
            obj.checkNegativeDistance = checkNegDist;
        end
               
    end %  ExhaustiveSearcher constructor
    
    
end %for methods

    methods(Hidden)
       
        function s = toStruct(this)
            warnState  = warning('query','all');
            warning('off','MATLAB:structOnObject');
            cleanupObj = onCleanup(@() warning(warnState));
            s = struct;
            s.checkNegativeDistance = this.checkNegativeDistance;
            if isa(this.Distance,'function_handle')
                error(message('stats:NeighborSearcher:searcherToStruct:CustomDistanceMetricNotSupported'));
            end
            s.Distance = this.Distance;
            s.X = this.X;
            s.DistParameterValue = this.DistParameter;
            if strcmpi(s.Distance,'minkowski')
                s.DistParameterName = 'P';
            elseif strcmpi(s.Distance,'mahalanobis')
                s.DistParameterName = 'Cov';
            elseif strcmpi(s.Distance,'seuclidean')
                s.DistParameterName = 'Scale';
            end
            s.FromStructFcn = 'ExhaustiveSearcher.fromStruct';
        end  
    end
    methods(Hidden, Static)
        function name = matlabCodegenRedirect(~)
            name = 'stats.coder.searcher.ExhaustiveSearcher';
        end
        
        function obj = fromStruct(s)
            if ~isempty(s.DistParameterValue)
               obj = ExhaustiveSearcher(s.X,'Distance',s.Distance,s.DistParameterName,...
                            s.DistParameterValue,'checkNegativeDistance',s.checkNegativeDistance); 
            else
               obj = ExhaustiveSearcher(s.X,'Distance',s.Distance,'checkNegativeDistance',s.checkNegativeDistance); 
            end
            
        end
    end
end % classdef
