classdef CompactTreeBagger  < classreg.learning.internal.DisallowVectorOps
%COMPACTTREEBAGGER Compact ensemble of decision trees grown by bagging
%    The CompactTreeBagger object is a lightweight object that contains
%    the trees grown by TreeBagger.  CompactTreeBagger does not
%    preserve any information about how the decision trees were grown.
%    It does not contain the input data used for growing trees, nor
%    does it contain training parameters such as minimal leaf size or
%    number of variables sampled for each decision split at random.
%    The CompactTreeBagger object can only be used for predicting the
%    response of the trained ensemble given new data X, and other
%    related functions.
%
%    CompactTreeBagger lets you save the trained ensemble to disk, or use
%    it in any other way, while discarding  training data and various
%    parameters of the training configuration irrelevant for predicting
%    response of the fully grown ensemble.  This reduces storage and memory
%    requirements, especially for ensembles trained on large datasets.
%
%    When you use the TreeBagger constructor to grow trees, it creates a
%    CompactTreeBagger object. You can obtain the compact object from the
%    full TreeBagger object using TREEBAGGER/COMPACT method. You do not
%    create an instance of CompactTreeBagger directly.
%
%    CompactTreeBagger properties:
%      Trees                       - Decision trees in the ensemble.
%      NumTrees                    - Number of decision trees in the ensemble.
%      ClassNames                  - Names of classes.
%      PredictorNames              - Predictor names.
%      Method                      - Method used by trees (classification or regression).
%      DeltaCriterionDecisionSplit - Split criterion contributions for each predictor.
%      NumPredictorSplit           - Number of decision splits on each predictor.
%      SurrogateAssociation        - Predictor associations.
%      DefaultYfit                 - Default value returned by PREDICT.
%
%    CompactTreeBagger methods:
%      combine        - Combine two ensembles.
%      error          - Error (misclassification probability or MSE).
%      margin         - Classification margin.
%      mdsprox        - Multidimensional scaling of proximity matrix.
%      meanMargin     - Mean classification margin per tree.
%      outlierMeasure - Outlier measure for data.
%      predict        - Predict response.
%      proximity      - Proximity matrix for data.
%      setDefaultYfit - Set the default value for PREDICT.
%
%   See also TREEBAGGER, TREEBAGGER/TREEBAGGER, TREEBAGGER/COMPACT,
%   ClassificationTree, RegressionTree,
%   classreg.learning.classif.CompactClassificationTree,
%   classreg.learning.regr.CompactRegressionTree.
    
%   Copyright 2008-2020 The MathWorks, Inc.

    properties(SetAccess=protected,GetAccess=public)
        %TREES Decision trees in the ensemble.
        %   The Trees property is a cell array of size NumTrees-by-1
        %   containing the trees in the ensemble.
        %
        %   See also COMPACTTREEBAGGER, NumTrees.
        Trees;
        
        %CLASSNAMES Names of classes.
        %   The ClassNames property is a cell array containing the class names for
        %   the response variable Y that was supplied to TreeBagger.  This property
        %   is empty for regression trees.
        %
        %   See also COMPACTTREEBAGGER.
        ClassNames;

        %DEFAULTYFIT Default value returned by PREDICT.
        %   The DefaultYfit property controls what predicted value is
        %   returned when no prediction is possible, for example when the
        %   PREDICT method needs to predict for an observation which has
        %   only false values in the matrix supplied through
        %   'UseInstanceForTree' argument.  For classification, you can set
        %   this property to either '' or 'MostPopular'.  If you choose
        %   'MostPopular' (default), the property value becomes the name of
        %   the most probable class in the training data. For regression,
        %   you can set this property to any numeric scalar.  The default
        %   is the mean of the response for the training data.
        %
        %   See also COMPACTTREEBAGGER, PREDICT, SETDEFAULTYFIT.
        DefaultYfit;
    end
    
    properties(SetAccess=protected,GetAccess=public,Dependent=true)
        %NUMTREES Number of decision trees in the ensemble.
        %   The NumTrees property is a scalar equal to the number of
        %   decision trees in the ensemble.
        %
        %   See also COMPACTTREEBAGGER, TREES.
        NumTrees;
        
        %PREDICTORNAMES Predictor names.
        %   The PredictorNames property is a cell array containing the names of the
        %   predictor variables (features). For an ensemble trained on a table,
        %   these are taken from the variable names in the table. For an ensemble
        %   trained on a matrix, these names are taken from the optional 'names'
        %   parameter, and the default names are 'x1', 'x2', etc. 
        %
        %   See also COMPACTTREEBAGGER.
        PredictorNames;
                
        %METHOD Method used by trees (classification or regression).
        %   The Method property is 'classification' for classification ensembles
        %   and 'regression' for regression ensembles.
        %
        %   See also COMPACTTREEBAGGER.
        Method;

        %DELTACRITERIONDECISIONSPLIT Split criterion contributions for each predictor.
        %   The DeltaCriterionDecisionSplit property is a numeric array of size
        %   1-by-Nvars of changes in the split criterion summed over splits on 
        %   each variable, averaged across the entire ensemble of grown trees.
        %
        % See also COMPACTTREEBAGGER, ClassificationTree/predictorImportance,
        % RegressionTree/predictorImportance.
        DeltaCriterionDecisionSplit;
        
        %NUMPREDICTORSPLIT Number of decision splits for each predictor.
        %   The NumPredictorSplit property is a numeric array of size 1-by-Nvars, where
        %   every element gives a number of splits on this predictor summed over
        %   all trees.
        %
        %   See also COMPACTTREEBAGGER.
        NumPredictorSplit;
        
        %SURROGATEASSOCIATION Predictor associations.
        %   The SurrogateAssociation property is a matrix of size Nvars-by-Nvars with
        %   predictive measures of predictor association, averaged across the entire
        %   ensemble of grown trees. If you grew the ensemble setting 'surrogate'
        %   to 'on', this matrix for each tree is filled with predictive measures
        %   of association averaged over the surrogate splits. If you grew the
        %   ensemble setting 'surrogate' to 'off' (default), SurrogateAssociation is diagonal.
        %
        % See also COMPACTTREEBAGGER, ClassificationTree/surrogateAssociation,
        % RegressionTree/surrogateAssociation.
        SurrogateAssociation ;
    end
    
    properties(SetAccess=public,GetAccess=public,Hidden=true)
        ClassProb;
        DefaultScore;
    end
    
    properties(SetAccess=protected,GetAccess=public,Hidden=true)
        NTrees
        VarNames;
    end
    
    properties(SetAccess=protected,GetAccess=public,Dependent=true,Hidden=true)
        DeltaCritDecisionSplit;
        NVarSplit
        VarAssoc
    end
    
    
    methods
        function n = get.NumTrees(bagger)
            n = bagger.NTrees;
        end
        
        function vnames = get.PredictorNames(bagger)
            vnames = bagger.VarNames;
        end
        
        function meth = get.Method(bagger)
            if isempty(bagger.ClassNames)
                meth = 'regression';
            else
                meth = 'classification';
            end
        end

        function deltacrit = get.DeltaCritDecisionSplit(bagger)
            nvar = length(bagger.VarNames);
            deltacrit = zeros(1,nvar);
            nTrees = bagger.NTrees;
            if nTrees==0
                return;
            end
            if isa(bagger.Trees{1},'classregtree')
                for it=1:nTrees
                    deltacrit = deltacrit + varimportance(bagger.Trees{it});
                end
            else
                for it=1:nTrees
                    deltacrit = deltacrit + predictorImportance(bagger.Trees{it});
                end
            end
            deltacrit = deltacrit/nTrees;
        end
             
        function deltacrit = get.DeltaCriterionDecisionSplit(bagger)
            deltacrit = bagger.DeltaCritDecisionSplit;
        end
              
        function nsplit = get.NVarSplit(bagger)
            nvar = length(bagger.VarNames);
            nsplit = zeros(1,nvar);
            nTrees = bagger.NTrees;
            if nTrees==0
                return;
            end
            if isa(bagger.Trees{1},'classregtree')
                if isempty(bagger.Trees{1}.nvarsplit)
                    nsplit = [];
                else
                    for it=1:nTrees
                        nsplit = nsplit + bagger.Trees{it}.nvarsplit;
                    end
                end
            else
                for it=1:nTrees
                    [~,nsplitOneTree] = predictorImportance(bagger.Trees{it}.Impl);
                    nsplit = nsplit + nsplitOneTree;
                end
            end
        end
        
        function nsplit = get.NumPredictorSplit(bagger)
            nsplit = bagger.NVarSplit;
        end
        
        function assoc = get.VarAssoc(bagger)
            nvar = length(bagger.VarNames);
            assoc = zeros(nvar);
            nTrees = bagger.NTrees;
            if nTrees==0
                return;
            end
            if isa(bagger.Trees{1},'classregtree')
                for it=1:nTrees
                    assoc = assoc + meansurrvarassoc(bagger.Trees{it});
                end
            else
                for it=1:nTrees
                    assoc = assoc + meanSurrVarAssoc(bagger.Trees{it});
                end
            end
            assoc = assoc/nTrees;
        end
        
        function assoc = get.SurrogateAssociation(bagger)
            assoc = bagger.VarAssoc;
        end
        
        function bagger = set.DefaultScore(bagger,score)
            if ~isnumeric(score)
                error(message('stats:CompactTreeBagger:DefaultScore:NotNumeric'));
            end
            if     bagger.Method(1)=='c' && ...
                    length(score)~=length(bagger.ClassNames) %#ok<MCSUP>
                error(message('stats:CompactTreeBagger:DefaultScore:SizeMismatch'));
            elseif bagger.Method(1)=='r' && ~isscalar(score) %#ok<MCSUP>
                error(message('stats:CompactTreeBagger:DefaultScore:MustBeScalarForRegr'));
            end
            if size(score,1)~=1
                score = score';
            end
            bagger.DefaultScore = score;
        end
    end
       
    methods
        function [yfit,scores,stdevs] = predict(bagger,X,varargin)
            %PREDICT Predict response.
            %   YFIT = PREDICT(B,X) computes the predicted response of the trained
            %   ensemble B for predictors X.  X should be a table if B was originally
            %   trained on a table, or a matrix if B was originally trained on a
            %   matrix. The output YFIT has one prediction for each row of X. YFIT is a
            %   cell array of strings for classification and a numeric array for
            %   regression. By default, PREDICT takes a democratic (non-weighted)
            %   average vote from all trees in the ensemble.
            %
            %   For regression, [YFIT,STDEVS] = PREDICT(B,X) also returns standard
            %   deviations of the computed responses over the ensemble of the grown
            %   trees.
            %
            %   For classification, [YFIT,SCORES] = PREDICT(B,X) returns scores for all
            %   classes. SCORES is a matrix with one row per observation and one column
            %   per class. The order of the columns is given by ClassNames property.
            %   For each observation and each class, the score generated by each tree
            %   is the probability of this observation originating from this class
            %   computed as the fraction of observations of this class in a tree leaf.
            %   These scores are averaged over all trees in the ensemble.
            %
            %   [YFIT,SCORES,STDEVS] = PREDICT(B,X) also returns standard deviations of
            %   the computed scores for classification. STDEVS is a matrix with one row
            %   per observation and one column per class, with standard deviations
            %   taken over the ensemble of the grown trees.
            %
            %   Y = PREDICT(B,X,'PARAM1',val1,'PARAM2',val2,...) specifies optional
            %   parameter name/value pairs:
            %
            %      'Trees'     Array of tree indices to use for computation of
            %                  responses.  Default is 'all'.
            %      'TreeWeights'  Array of NumTrees weights for weighting votes from the
            %                  specified trees.
            %      'UseInstanceForTree'  Logical matrix of size Nobs-by-NumTrees 
            %                  indicating which trees should be used to
            %                  make predictions for each observation.  By
            %                  default all trees are used for all
            %                  observations. If no trees are used for an
            %                  observation (the row of the matrix has no
            %                  true values), PREDICT returns
            %                  the default prediction, most popular class
            %                  for classification or sample mean for
            %                  regression. The default prediction can be
            %                  changed using the DefaultYfit property.
            %
            %   See also COMPACTTREEBAGGER, CLASSNAMES, ClassificationTree/predict,
            %   RegressionTree/predict, DEFAULTYFIT, SETDEFAULTYFIT, Method.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            if istall(X)
               % Process inputs
               
               %Disallow useinstancefortree, trees  or treeweights for tall
               args = {{'useifort' 'useinstancefortree'} 'trees' 'treeweights'};
               defs = {                               []   []            []};
               [useIforT,useTrees,treeWeights,~, extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
               if ~isempty(useIforT)||~isempty(useTrees)||~isempty(treeWeights)
                     error(message('stats:tall:TreeBagger:NonSupportOpt'))
               end
               
               [yfit,scores,stdevs] = hSlicefun(@bagger.predict,X,extraArgs{:});
               return;
            end
               
            % Process inputs
            args = {{'useifort' 'useinstancefortree'} 'trees' 'treeweights'};
            defs = {                            'all'   'all'            []};
            [useIforT,useTrees,treeWeights] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Compute response
            if nargout>2 || (bagger.Method(1)=='r' && nargout>1)
                [yfit,scores,stdevs,scoreweights] = ...
                    predictAccum(bagger,X,'useifort',useIforT,...
                    'trees',useTrees,'treeweights',treeWeights);
                stdevs = bsxfun(@rdivide,stdevs,scoreweights);
                stdevs = sqrt(max(stdevs,0));
                if bagger.Method(1)=='r'
                    scores = stdevs;
                end
            else
                [yfit,scores] = ...
                    predictAccum(bagger,X,'useifort',useIforT,...
                    'trees',useTrees,'treeweights',treeWeights);
            end
        end
        
        function [err,sumW] = error(bagger,X,varargin)
            %ERROR Error (misclassification probability or MSE).
            %   ERR = ERROR(B,X,Y) computes the misclassification probability (for
            %   classification trees) or mean squared error (MSE, for regression trees)
            %   for predictors X given true response Y. X must be a table if B was
            %   originally trained on a table, or a matrix if B was originally trained
            %   on a matrix. For classification, Y can be either a numeric vector, cell
            %   array of strings, categorical vector, or logical vector. For
            %   regression, Y must be a numeric vector. Y can be omitted if X is a
            %   table that includes the response variable. ERR is a vector with one
            %   error measure for each of the NumTrees trees in the ensemble B.
            %
            %   ERR = ERROR(B,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
            %   parameter name/value pairs.
            %
            %     'Mode'     String indicating how errors are computed. If set to
            %                'cumulative' (default), cumulative errors are computed and
            %                ERR is a vector of length NumTrees, where the 1st element
            %                gives error from tree 1, 2nd element gives error from
            %                trees 1:2 etc, up to 1:NumTrees. If set to 'individual',
            %                ERR is a vector of length NumTrees, where each element is
            %                an error from each tree in the ensemble. If set to
            %                'ensemble', ERR is a scalar showing the cumulative error
            %                for the entire ensemble.
            %     'Weights'  Vector of observation weights to use for error
            %                averaging. By default the weight of every observation
            %                is set to 1. The length of this vector must be equal
            %                to the number of rows in X.
            %     'Trees'    Vector of indices indicating what trees to include
            %                in this calculation. By default, this argument is set to
            %                'all' and ERROR uses all trees. If 'Trees' is a numeric
            %                vector, the method returns a vector of length NumTrees for
            %                'cumulative' and 'individual' modes, where NumTrees is the
            %                number of elements in the input vector, and a scalar for
            %                'ensemble' mode. For example, in the 'cumulative' mode,
            %                the 1st element gives an error from trees(1), the 2nd
            %                element gives an error from trees(1:2) etc.
            %     'TreeWeights' Vector of tree weights. This vector must have the same
            %                length as the 'Trees' vector. ERROR uses these weights to
            %                combine output from the specified trees by taking a
            %                weighted average instead of the simple non-weighted
            %                majority vote. You cannot use this argument in the
            %                'individual' mode.
            %     'UseInstanceForTree' Logical matrix of size Nobs-by-NumTrees
            %                indicating which trees ERROR should use to
            %                make predictions for each observation.  By
            %                default the method uses all trees for all
            %                observations.
            %
            % See also COMPACTTREEBAGGER, PREDICT, NumTrees.
           
            [varargin{:}] = convertStringsToChars(varargin{:});
            if  istall(X) || any(cellfun(@(x) istall(x),varargin)) 
               err = errorTall(bagger,X,varargin{:});
               return;
            end
            
            % Get Y from input or table as required 
            [Y,varargin] = inferResponse(bagger,X,varargin{:});

            % Process inputs
            args = {      'mode' 'weights' {'useifort' 'useinstancefortree'} ...
                'trees' 'treeweights'};
            defs = {'cumulative'        []                             'all' ...
                  'all'            []};
            [mode,W,useIforT,useTrees,treeWeights] = ...
                internal.stats.parseArgs(args,defs,varargin{:});

            % Process mode
            if ~ischar(mode) || ...
                    ~any(strncmpi(mode,{'cumulative' 'individual' 'ensemble'},length(mode)))
                error(message('stats:CompactTreeBagger:error:BadMode'));
            end
            
            % Check observation weights
            N = size(X,1);
            if isempty(W)
                W = ones(N,1);
            end
            if ~isnumeric(W) || length(W)~=N || any(W<0) || all(W==0)
                error(message('stats:CompactTreeBagger:predictAccum:BadW', N));
            end

            % Convert logical to double
            if islogical(Y)
                Y = double(Y);
            end
            
            % For classification, convert labels to strings
            if bagger.Method(1)=='c' && ~iscellstr(Y)
                Y = cellstr(classreg.learning.internal.ClassLabel(Y));
            end
            
            % Check trees
            if ~strcmpi(useTrees,'all') && ~isnumeric(useTrees)
                error(message('stats:CompactTreeBagger:error:BadUseTrees'));
            end
            
            % How many trees?
            if ischar(useTrees)
                useTrees = 1:bagger.NTrees;
            end
            nTrees = length(useTrees);
            if strncmpi(mode,'ensemble',length(mode))
                nTrees = 1;
            end

            % Process weights
            if isempty(treeWeights)
                treeWeights = ones(length(useTrees),1);
            else
                if strncmpi(mode,'individual',length(mode))
                    warning(message('stats:CompactTreeBagger:error:TreeWeightsNotAllowedForIndivMode'));
                end
                if ~isnumeric(treeWeights)
                    error(message('stats:CompactTreeBagger:error:TreeWeightsNotNumeric'));
                end
            end

            % Init
            nCols = max(length(bagger.ClassNames),1);
            sfit = NaN(N,nCols);
            stdfit = NaN(N,nCols);
            wfit = zeros(N,1);
            sumW = zeros(nTrees,1);
            % Init errors
            err = NaN(nTrees,1);
            
            % Compute
            if     strncmpi(mode,'individual',length(mode))
                for it=1:nTrees
                    yfit = predict(bagger,X,'useifort',useIforT,'trees',useTrees(it));
                    if bagger.Method(1)=='c' % classification
                        filled = ~cellfun(@isempty,yfit);
                        Wtot = sum(W(filled));
                        if Wtot>0
                            err(it) = dot(W(filled),~strcmp(Y(filled),yfit(filled)))/Wtot;
                            sumW(it) = Wtot;
                        end
                    else
                        filled = ~arrayfun(@isnan,yfit);
                        Wtot = sum(W(filled));
                        if Wtot>0
                            err(it) = dot(W(filled),(Y(filled)-yfit(filled)).^2)/Wtot;
                            sumW(it) = Wtot;
                        end
                    end
                end
            elseif strncmpi(mode,'cumulative',length(mode))
                for it=1:nTrees
                    [yfit,sfit,stdfit,wfit] = ...
                        predictAccum(bagger,X,'useifort',useIforT,...
                        'trees',useTrees(it),'treeweights',treeWeights(it),...
                        'scores',sfit,'stdevs',stdfit,'scoreweights',wfit);
                    if bagger.Method(1)=='c' % classification
                        filled = ~cellfun(@isempty,yfit);
                        Wtot = sum(W(filled));
                        if Wtot>0
                            err(it) = dot(W(filled),~strcmp(Y(filled),yfit(filled)))/Wtot;
                            sumW(it) = Wtot;
                        end
                    else
                        filled = ~arrayfun(@isnan,yfit);
                        Wtot = sum(W(filled));
                        if Wtot>0
                            err(it) = dot(W(filled),(Y(filled)-yfit(filled)).^2)/Wtot;
                            sumW(it) = Wtot;
                        end
                    end
                end
            elseif strncmpi(mode,'ensemble',length(mode))
                yfit = predict(bagger,X,'useifort',useIforT,...
                    'trees',useTrees,'treeweights',treeWeights);
                if bagger.Method(1)=='c' % classification
                    filled = ~cellfun(@isempty,yfit);
                    Wtot = sum(W(filled));
                    if Wtot>0
                        err = dot(W(filled),~strcmp(Y(filled),yfit(filled)))/Wtot;
                        sumW = Wtot;
                    end
                else
                    filled = ~arrayfun(@isnan,yfit);
                    Wtot = sum(W(filled));
                    if Wtot>0
                        err = dot(W(filled),(Y(filled)-yfit(filled)).^2)/Wtot;
                        sumW = Wtot;
                    end
                end
            end
        end
        
        function mar = margin(bagger,X,varargin)
            %MARGIN Classification margin.
            %   MAR = MARGIN(B,X,Y) computes the classification margins for predictors
            %   X given true response Y. X must be a table if B was originally trained
            %   on a table, or a matrix if B was originally trained on a matrix. Y can
            %   be either a numeric vector, cell array of strings, categorical vector,
            %   or logical vector.  Y can be omitted if X is a table that includes the
            %   response variable. MAR is a numeric array of size Nobs-by-NumTrees,
            %   where Nobs is the number of rows of X and Y, and NumTrees is the number
            %   of trees in the ensemble B.  For observation I and tree J, MAR(I,J) is
            %   the difference between the score for the true class and the largest
            %   score for other classes.  This method is available for classification
            %   ensembles only.
            %
            %   MAR = MARGIN(B,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
            %   parameter name/value pairs.
            %
            %     'Mode'     String indicating how MARGIN computes errors. If set to
            %                'cumulative' (default), the method computes cumulative
            %                margins and MAR is an Nobs-by-NumTrees matrix, where the 1st
            %                column gives margins from tree 1, 2nd column gives margins
            %                from trees 1:2 etc, up to 1:NumTrees. If set to
            %                'individual', MAR is an Nobs-by-NumTrees matrix, where each
            %                column gives margins from each tree in the ensemble. If
            %                set to 'ensemble', MAR is a single column of length Nobs
            %                showing the cumulative margins for the entire ensemble.
            %     'Trees'    Vector of indices indicating what trees to include
            %                in this calculation. By default, this argument is set to
            %                'all' and MARGIN uses all trees. If 'Trees' is a numeric
            %                vector, the method returns an Nobs-by-NumTrees matrix for
            %                'cumulative' and 'individual' modes, where NumTrees is the
            %                number of elements in the input vector, and a single
            %                column for 'ensemble' mode. For example, in the
            %                'cumulative' mode, the 1st column gives margins from
            %                trees(1), the 2nd column gives margins from trees(1:2)
            %                etc.
            %     'TreeWeights' Vector of tree weights. This vector must have the same
            %                length as the 'Trees' vector. MARGIN uses these weights to
            %                combine output from the specified trees by taking a
            %                weighted average instead of the simple non-weighted
            %                majority vote. This argument cannot be used in the
            %                'individual' mode.
            %     'UseInstanceForTree' Logical matrix of size Nobs-by-NumTrees 
            %                indicating which trees should be used to make
            %                predictions for each observation.  By default
            %                all trees are used for all observations.
            %
            % See also COMPACTTREEBAGGER, PREDICT, GROUPINGVARIABLE.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            % Get Y from input or table as required 
            if   istall(X) || any(cellfun(@(x) istall(x),varargin)) 
       
                [Y,varargin] = inferResponseTall(bagger,X,varargin{:});
                if ~istall(Y)
                    error(message('stats:tall:classreg:ResponseNotTall'))
                end
                
               %Disallow useinstancefortree, trees  or treeweights for tall
               args = {{'useifort' 'useinstancefortree'} 'trees' 'treeweights'};
               defs = {                               []   []            []};
               [useIforT,useTrees,treeWeights,~, extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
               if ~isempty(useIforT)||~isempty(useTrees)||~isempty(treeWeights)
                     error(message('stats:tall:TreeBagger:NonSupportOpt'))
               end
   
               mar = hSlicefun(@(x,y,varargin)localMargin(bagger,x,y,varargin{:}), X,Y,extraArgs{:});

               return;
            end
            [Y,varargin] = inferResponse(bagger,X,varargin{:});
            
            % Process inputs
            args = {      'mode' {'useifort' 'useinstancefortree'} 'trees' 'treeweights'};
            defs = {'cumulative'                             'all'   'all'            []};
            [mode,useIforT,useTrees,treeWeights] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Process mode
            if ~ischar(mode) || ...
                    ~any(strncmpi(mode,{'cumulative' 'individual' 'ensemble'},length(mode)))
                error(message('stats:CompactTreeBagger:margin:BadMode'));
            end
            
            % This Method is valid for classification only
            if bagger.Method(1)~='c'
                error(message('stats:CompactTreeBagger:margin:InvalidOperation'));
            end
            
            % For classification, convert labels to strings
            if islogical(Y)
                Y = double(Y);
            end
            Y = cellstr(classreg.learning.internal.ClassLabel(Y));
            
            % Check trees
            if ~strcmpi(useTrees,'all') && ~isnumeric(useTrees)
                error(message('stats:CompactTreeBagger:margin:BadUseTrees'));
            end
            
            % How many trees?
            if ischar(useTrees)
                useTrees = 1:bagger.NTrees;
            end
            nTrees = length(useTrees);
            if strncmpi(mode,'ensemble',length(mode))
                nTrees = 1;
            end

            % Process weights
            if isempty(treeWeights)
                treeWeights = ones(length(useTrees),1);
            else
                if strncmpi(mode,'individual',length(mode))
                    warning(message('stats:CompactTreeBagger:margin:TreeWeightsForIndivModeNotAllowed'));
                end
                if ~isnumeric(treeWeights)
                    error(message('stats:CompactTreeBagger:margin:TreeWeightsNotNumeric'));
                end
            end

            % Init
            N = size(X,1);
            nCols = length(bagger.ClassNames);
            sfit = NaN(N,nCols);
            stdfit = NaN(N,nCols);
            wfit = zeros(N,1);
            
            % Init margin
            mar = NaN(N,nTrees);
            
            % Get positions of true classes in the scores matrix
            tf = ismember(Y,bagger.ClassNames);
            if any(~tf)
                error(message('stats:CompactTreeBagger:margin:BadLabels'));
            end
            
            % Get a matrix of class counts
            C = classreg.learning.internal.classCount(bagger.ClassNames,Y);
            
            % Compute
            if     strncmpi(mode,'individual',length(mode))
                for it=1:nTrees
                    [~,sfit] = predict(bagger,X,...
                        'useifort',useIforT,'trees',useTrees(it));
                    mar(:,it) = classreg.learning.loss.classifmargin(C,sfit);
                end
            elseif strncmpi(mode,'cumulative',length(mode))
                for it=1:nTrees
                    [~,sfit,stdfit,wfit] = ...
                        predictAccum(bagger,X,'useifort',useIforT,...
                        'trees',useTrees(it),'treeweights',treeWeights(it),...
                        'scores',sfit,'stdevs',stdfit,'scoreweights',wfit);
                    mar(:,it) = classreg.learning.loss.classifmargin(C,sfit);
                end
            elseif strncmpi(mode,'ensemble',length(mode))
                 [~,sfit] = predict(bagger,X,'useifort',useIforT,...
                    'trees',useTrees,'treeweights',treeWeights);
                mar = classreg.learning.loss.classifmargin(C,sfit);
            end
        end
        
        function mar = meanMargin(bagger,X,varargin)
            %MEANMARGIN Mean classification margin.
            %   MAR = MEANMARGIN(B,X,Y) computes average classification margins for
            %   predictors X given true response Y. X must be a table if B was
            %   originally trained on a table, or a matrix if B was originally trained
            %   on a matrix. Y can be either a numeric vector, cell array of strings,
            %   categorical vector, or logical vector.  Y can be omitted if X is a
            %   table that includes the response variable. MEANMARGIN averages the
            %   margins over all observations (rows) in X for each tree.  MAR is a
            %   matrix of size 1-by-NumTrees, where NumTrees is the number of trees in
            %   the ensemble B. This method is available for classification ensembles
            %   only.
            %
            %   MAR = MEANMARGIN(B,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies
            %   optional parameter name/value pairs.
            %
            %     'Mode'     String indicating how MEANMARGIN computes errors. If set
            %                to 'cumulative' (default), the method computes cumulative
            %                errors and MAR is a vector of length NumTrees, where the 1st
            %                element gives mean margin from tree 1, 2nd element gives
            %                mean margin from trees 1:2 etc, up to 1:NumTrees. If set to
            %                'individual', MAR is a vector of length NumTrees, where each
            %                element is a mean margin from each tree in the ensemble.
            %                If set to 'ensemble', MAR is a scalar showing the
            %                cumulative mean margin for the entire ensemble.
            %     'Weights'  Vector of observation weights to use for margin
            %                averaging. By default the weight of every observation
            %                is set to 1. The length of this vector must be equal
            %                to the number of rows in X.
            %     'Trees'    Vector of indices indicating what trees to include
            %                in this calculation. By default, this argument is set to
            %                'all' and all trees are used. If 'Trees' is a numeric
            %                vector, the method returns a vector of length NumTrees for
            %                'cumulative' and 'individual' modes, where NumTrees is the
            %                number of elements in the input vector, and a scalar for
            %                'ensemble' mode. For example, in the 'cumulative' mode,
            %                the 1st element gives mean margin from tree trees(1), the
            %                2nd element gives mean margin from trees trees(1:2) etc.
            %     'TreeWeights' Vector of tree weights. This vector must have the same
            %                length as the 'Trees' vector. These weights are used to
            %                combine output from the specified trees by taking a
            %                weighted average instead of the simple non-weighted
            %                majority vote. This argument cannot be used in the
            %                'individual' mode.
            %     'UseInstanceForTree' Logical matrix of size Nobs-by-NumTrees 
            %                indicating which trees should be used to make
            %                predictions for each observation.  By default
            %                all trees are used for all observations.
            %
            %   See also COMPACTTREEBAGGER, PREDICT.
           
            [varargin{:}] = convertStringsToChars(varargin{:});
            if   istall(X) || any(cellfun(@(x) istall(x),varargin)) 
               mar = meanMarginTall(bagger,X,varargin{:});
               return;
            end
            % Get Y from input or table as required 
            [Y,varargin] = inferResponse(bagger,X,varargin{:});

            % Process inputs
            args = {'weights'};
            defs = {       []};
            [W,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check observation weights
            N = size(X,1);
            if isempty(W)
                W = ones(N,1);
            end
            if ~isnumeric(W) || length(W)~=N || any(W<0) || all(W==0)
                error(message('stats:CompactTreeBagger:meanMargin:BadW', N));
            end
            
            % Compute mean margins
            m = margin(bagger,X,Y,extraArgs{:});
            Ncol = size(m,2);
            mar = NaN(1,Ncol);
            for ncol=1:Ncol
                tf = ~isnan(m(:,ncol));
                Wtot = sum(W(tf));
                if Wtot>0
                    mar(ncol) = dot(W(tf),m(tf,ncol))/Wtot;
                end
            end
        end

        function prox = proximity(bagger,X,varargin)
            %PROXIMITY Proximity matrix for data.
            %   PROX = PROXIMITY(B,X) computes a numeric matrix of size Nobs-by-Nobs of
            %   proximities for data X, where Nobs is the number of observations (rows)
            %   in X.  Proximity between any two observations in the input data is
            %   defined as a fraction of trees in the ensemble B for which these two
            %   observations land on the same leaf.  This is a symmetric matrix with
            %   1's on the diagonal and off-diagonal elements ranging from 0 to 1.
            %
            %   See also COMPACTTREEBAGGER.
             if   istall(X) || any(cellfun(@(x) istall(x),varargin)) 
                 error('PROXIMITY does not support tall arrays.');
             end
            prox = squareform(flatprox(bagger,X,varargin{:}));
            N = size(X,1);
            prox(1:N+1:end) = 1;
        end
        
        function outlier = outlierMeasure(bagger,data,varargin)
            %OUTLIERMEASURE Outlier measure for data.
            %   OUT = OUTLIERMEASURE(B,X) computes outlier measures for predictors X
            %   using trees in the ensemble B.  X must be a table if B was originally
            %   trained on a table, or a matrix if B was originally trained on a
            %   matrix. outlierMeasure computes the outlier measure for a given
            %   observation by taking an inverse of the average squared proximity
            %   between this observation and other observations. OUTLIERMEASURE then
            %   normalizes these outlier measures by subtracting the median of their
            %   distribution, taking the absolute value of this difference, and
            %   dividing by the median absolute deviation. A high value of the outlier
            %   measure indicates that this observation is an outlier.
            %
            %   You can supply the proximity matrix directly by using the 'Data'
            %   parameter.
            %
            %   OUT = OUTLIERMEASURE(B,X,'PARAM1',val1,'PARAM2',val2,...) specifies
            %   optional parameter name/value pairs:
            %
            %     'Data'     Flag indicating how to treat the X input argument.  If
            %                set to 'predictors' (default), the method assumes X is a
            %                matrix of predictors and uses it for computation of the
            %                proximity matrix. If set to 'proximity', the method treats
            %                X as a proximity matrix returned by the PROXIMITY method.
            %                If you do not supply the proximity matrix, OUTLIERMEASURE
            %                computes it internally.  If you use the PROXIMITY method
            %                to compute a proximity matrix, supplying it as input to
            %                OUTLIERMEASURE reduces computing time.
            %     'Labels'   Vector of true class labels for a classification
            %                ensemble. True class labels can be either a numeric
            %                vector, character matrix, cell array of strings,
            %                categorical vector or logical vector (see help for
            %                groupingvariable). When you supply this parameter, the
            %                method performs the outlier calculation for any
            %                observations using only other observations from the same
            %                class.  This parameter must specify one label for each
            %                observation (row) in X.
            %
            %   See also COMPACTTREEBAGGER, PROXIMITY, GROUPINGVARIABLE.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            % Process inputs
              if   istall(data) 
                 error('OUTLIERMEASURE does not support tall arrays.');
              end
             
            args = {      'data' 'labels'};
            defs = {'predictors'       []};
            [datatype,Y] = internal.stats.parseArgs(args,defs,varargin{:});

            % Check if Y input is used for regression
            if bagger.Method(1)=='r' && ~isempty(Y)
                error(message('stats:CompactTreeBagger:outlierMeasure:LabelsForRegressionNotAllowed'));
            end
            
            % Check the type of input data
            if ~ischar(datatype) || ...
                    ~any(strncmpi(datatype,{'predictors' 'proximity'},length(datatype)))
                error(message('stats:CompactTreeBagger:outlierMeasure:BadDataArg'));
            end

            % Get proximity and predictors
            prox = [];
            X = [];
            if strncmpi(datatype,'predictors',length(datatype))
                X = data;
                N = size(X,1);
            end
            if strncmpi(datatype,'proximity',length(datatype))
                prox = data;
                [n1,n2] = size(prox);
                if n1~=n2
                    error(message('stats:CompactTreeBagger:outlierMeasure:ProxNotSquare'));
                end
                N = n1;
            end
            
            % For classification, convert labels to strings
            if islogical(Y)
                Y = double(Y);
            end
            if bagger.Method(1)=='c' && ~isempty(Y) && ~iscellstr(Y)
                Y = cellstr(classreg.learning.internal.ClassLabel(Y));
            end
            
            % Init
            outlier = zeros(N,1);

            % Check size
            if bagger.Method(1)=='c' && ~isempty(Y) && N~=length(Y)
                error(message('stats:CompactTreeBagger:outlierMeasure:XYSizeMismatch'));
            end
            
            % Split into classes
            if bagger.Method(1)=='c' && ~isempty(Y)
                [tf,truepos] = ismember(Y,bagger.ClassNames);
                if any(~tf)
                    error(message('stats:CompactTreeBagger:outlierMeasure:BadLabels'));
                end
            else
                truepos = ones(N,1);
            end
            nclass = max(truepos);
            
            % Process input proximity matrix if any was supplied
            if ~isempty(prox)
                [n1,n2] = size(prox);
                if n1~=n2 || n1~=N
                    error(message('stats:CompactTreeBagger:outlierMeasure:ProxSizeMismatch', N, N));
                end
            end
            
            % Compute proximities if they were not supplied
            if isempty(prox)
                prox = proximity(bagger,X);
            end

            % Process proximities for each class separately
            for c=1:nclass
                isc = truepos==c;
                nobs = sum(isc);
                if nobs>2 % Can't compute outliers for two or less distances
                    outlier(isc) = nobs./sum(prox(isc,isc).^2,2);
                    m = median(outlier(isc));
                    madm = mad(outlier(isc),1);% median abs dev
                    if madm>0
                        outlier(isc) = abs(outlier(isc)-m)/madm;
                    end
                end
            end
        end
        
        function [sc,eigen] = mdsprox(bagger,data,varargin)
            %MDSPROX Multidimensional scaling of proximity matrix.
            %   [SC,EIGEN] = MDSPROX(B,X) applies classical multidimensional scaling to
            %   the proximity matrix computed for the data in X, and returns scaled
            %   coordinates SC and eigenvalues EIGEN of the scaling transformation. X
            %   must be a table if B was originally trained on a table, or a matrix if
            %   B was originally trained on a matrix. MDSPROX applies multidimensional
            %   scaling to the matrix of distances defined as 1-PROX, where PROX is the
            %   proximity matrix returned by the PROXIMITY method.
            %
            %   You can supply the proximity matrix directly by using the 'Data'
            %   parameter.
            %
            %   [SC,EIGEN] = MDSPROX(B,X,'PARAM1',val1,'PARAM2',val2,...) specifies
            %   optional parameter name/value pairs:
            %
            %     'Data'   Flag indicating how the method treats the X input argument.
            %              If set to 'predictors' (default), MDSPROX assumes X to be a
            %              table or matrix of predictors to be used for computation of
            %              the proximity matrix. If set to 'proximity', the method
            %              treats X as a proximity matrix returned by the PROXIMITY
            %              method.
            %     'Colors' If you supply this argument, MDSPROX makes overlaid scatter
            %              plots of two scaled coordinates using specified colors for
            %              different classes. You must supply the colors as a string
            %              with one character for each color.  If there are more
            %              classes in the data than characters in the supplied string,
            %              MDSPROX only plots the first C classes, where C is the
            %              length of the string. For regression or if you do not
            %              provide the vector of true class labels, MDSPROX uses the
            %              first color for all observations in X.
            %     'Labels' Vector of true class labels for a classification
            %              ensemble. True class labels can be either a numeric vector,
            %              character matrix, cell array of strings, categorical vector
            %              or logical vector (see help for groupingvariable). If
            %              supplied, this vector must have as many elements as there
            %              are observations (rows) in X. This argument has no effect
            %              unless 'Colors' argument is supplied.
            %     'MDSCoordinates' Indices of the two or three scaled coordinates to be
            %              plotted. By default, MDSPROX makes a scatter plot of the 1st
            %              and 2nd scaled coordinates which correspond to the two
            %              largest eigenvalues.  You can specify any other two or three
            %              indices not exceeding the dimensionality of the scaled data.
            %              This argument has no effect unless you also supply the
            %              'Colors' argument.
            %
            %   See also COMPACTTREEBAGGER, PROXIMITY, CMDSCALE, GROUPINGVARIABLE.
            [varargin{:}] = convertStringsToChars(varargin{:});
             if   istall(data) 
                 error('MDSPROX does not support tall arrays.');
             end
            [sc,eigen] = mdsProx(bagger,data,varargin{:});    
        end
                
        function b1 = combine(b1,b2)
            %COMBINE Combine two ensembles.
            %   B1 = COMBINE(B1,B2) appends decision trees from ensemble B2 to those
            %   stored in B1 and returns ensemble B1. This method requires that the
            %   class and variable names be identical in both ensembles.
            %
            %   See also COMPACTTREEBAGGER.
            
            % Check classes
            if length(b1.ClassNames)~=length(b2.ClassNames) ...
                    || any(~strcmp(b1.ClassNames,b2.ClassNames))
                error(message('stats:CompactTreeBagger:combine:IncompatibleClasses'));
            end
            
            % Check vars
            if length(b1.VarNames)~=length(b2.VarNames) ...
                    || any(~strcmp(b1.VarNames,b2.VarNames))
                error(message('stats:CompactTreeBagger:combine:IncompatibleVarNames'));
            end

            % Add trees
            b1 = addTrees(b1,b2.Trees);
        end
        
        function bagger = setDefaultYfit(bagger,yfit)
            %SETDEFAULTYFIT Set the default value for PREDICT.
            %   B = SETDEFAULTYFIT(B,YFIT) sets the default prediction for ensemble B
            %   to YFIT. The default prediction must be a character variable for
            %   classification or a numeric scalar for regression. This setting
            %   controls what predicted value is returned when no prediction is
            %   possible, for example when the PREDICT method needs to predict for an
            %   observation which has only false values in the matrix supplied through
            %   'UseInstanceForTree' argument.
            %
            %   See also COMPACTTREEBAGGER, PREDICT, DEFAULTYFIT.
            
            yfit = convertStringsToChars(yfit);
            if     bagger.Method(1)=='c'
                if ~ischar(yfit) || ...
                        ~ismember(lower(yfit),{'' 'mostpopular'})
                    error(message('stats:CompactTreeBagger:setDefaultYfit:InvalidSetting'));
                end
                if strcmpi(yfit,'')
                    bagger.DefaultYfit = '';
                    bagger.DefaultScore = NaN(1,length(bagger.ClassNames));
                else
                    [~,i] = max(bagger.ClassProb);
                    bagger.DefaultYfit = bagger.ClassNames{i};
                    bagger.DefaultScore = bagger.ClassProb;
                end
            elseif bagger.Method(1)=='r'
                bagger.DefaultYfit = yfit;
                bagger.DefaultScore = yfit;
            end
        end
        
        function [pd,varargout] = partialDependence(bagger,features,varargin)
        %PARTIALDEPENDENCE Partial Dependence for one and two variables
        %   PD = partialDependence(MODEL,VARS,DATA) is for a fitted regression
        %   CompactTreeBagger. It computes the partial dependence PD between the
        %   predictor variables listed in VARS and the response of the model, by
        %   averaging over the data in the matrix or table DATA. You can specify
        %   one or two variables using VARS. VARS is either column indices of
        %   predictors or names of predictors. The names must match the entries in
        %   MODEL.PredictorNames.
        %
        %   PD = partialDependence(MODEL,VARS,LABELS,DATA) is for a fitted
        %   classification CompactTreeBagger. You can specify one or multiple
        %   classes using LABELS. The values in LABELS must match and have the same
        %   data type as the class names in MODEL.ClassNames. If there are multiple
        %   values in both VARS and LABELS, then PD is a three-dimensional array.
        %   PD is the partial dependence between the predictor variables listed in
        %   VARS and the scores for the classes specified in LABELS.
        %
        %   PD = partialDependence(...,Name,Value) specifies additional options
        %   using one or more name-value pair arguments:
        % 
        %       'NumObservationsToSample':      An integer specifying the number
        %                                       of rows to sample at random from the
        %                                       dataset (either the data specified by
        %                                       argument DATA or the data used to train
        %                                       MODEL). The default is to use all rows.
        % 
        %       'QueryPoints':                  The points at which to calculate
        %                                       the partial dependence. For the case of
        %                                       a single predictor, the value of
        %                                       'QueryPoints' must be a column vector.
        %                                       For the case of two predictors, the
        %                                       value of 'QueryPointsâ€' must be a 1x2
        %                                       cell array containing a separate vector
        %                                       for each predictor. The default is to
        %                                       use 100 points equally spaced across
        %                                       the range of the predictor.
        % 
        %       'UseParallel'                   Can be true to instruct the method to
        %                                       perform the averaging calculations in
        %                                       parallel (using parfor), or false
        %                                       (default) to instruct the method not to
        %                                       use parallel computation.
        % 
        %   [PD,X,Y] = partialDependence(...) also returns the arrays X and Y
        %   containing the query point values of the first and second predictor in VARS,
        %   respectively.
        % 
        %   Examples:
        %       % Partial Dependence for bagged decision trees
        %       load census1994
        %       X = adultdata(:,{'age','workClass','education_num','marital_status','race',...
        %           'sex','capital_gain','capital_loss','hours_per_week','salary'});
        %       f = TreeBagger(10,X,'salary');
        %       f = compact(f); 
        % 
        %       pd = partialDependence(f,'age','<=50K',X);
        % 
        %       % Obtain values x of the variable 'age'
        %       [pd,x] = partialDependence(f,'age',{'<=50K','>50K'},X);
        %
        %       % plot partial dependence of the first and second classes scores on
        %       % 'age'
        %       plot(x,pd)
        % 
        %       % Obtain values x and y of variables 'age' and 'workClass'
        %       [pd,x,y] = partialDependence(f,{'age','workClass'},{'<=50K','>50K'},X);
        % 
        %       % plot partial dependence for class '<=50K' on variables
        %       % 'age' and 'workClass'
        %       surf(x,y,pd(:,:,1))
        %
        %       % plot partial dependence for class '>50K' on variables
        %       % 'age' and 'workClass'
        %       surf(x,y,pd(:,:,2))
        % 
        %       % With optional name-value pairs
        %       pd = partialDependence(f,1,'<=50K',X,'NumObservationsToSample',100);
        %       pd = partialDependence(f,1,'<=50K',X,'UseParallel',true);
        % 
        %       % Plot the Individual Conditional Expectation
        %       pd = partialDependence(f,1,'<=50K',X,'Conditional','absolute');
        % 
        %       % Provide alternative query points
        %       xi = linspace(min(X{:,'age'}),max(X{:,'age'}))';
        %       pd = partialDependence(f,1,'<=50K',X,'QueryPoints',xi);
        % 
        %       % Provide alternative query points
        %       xi = linspace(min(X{:,'age'}),max(X{:,'age'}))';
        %       pd = partialDependence(f,1,'<=50K',X,'QueryPoints',xi);
        % 
        %       xi = cell(1,2);
        %       xi{1} = linspace(min(X{:,1}),max(X{:,1}))';
        %       xi{2} = linspace(min(X{:,9}),max(X{:,9}))';
        %       pd = partialDependence(f,[1,9],'<=50K',X,'QueryPoints',xi);

        %-------Check number of inputs----
        narginchk(3,10);
        features = convertStringsToChars(features);
        [varargin{:}] = convertStringsToChars(varargin{:});
        
        % Regression tree bagger
        if strcmp(bagger.Method,'regression')
            % class label is not an input for regression models
            % For a compact regression model the third argument is Data
            X = varargin{1};
            % Exclude Data from varargin 
            varargin(1) = [];
            % Classification flag is false and class label is empty
            classif = false;
            labels = [];
        % Classification tree bagger
        else
            % class label is an input for classification models
            if nargin < 4
                error(message('MATLAB:narginchk:notEnoughInputs'));
            else
            labels = varargin{1}; 
            X = varargin{2};
            % Exclude Data from varargin
            varargin(1:2) = [];
            % classification flag is true
            classif = true;
            end
        end     
        
        % Call the function from classreg.learning.impl
        [pd,x,y] = classreg.learning.impl.partialDependence(bagger,features,...
            labels,X,classif,varargin{:});
        varargout{1} = x;
        varargout{2} = y;
        end
        
        function [AX] = plotPartialDependence(bagger,features,varargin)
        %PLOTPARTIALDEPENDENCE Partial Dependence Plot for 1-D or 2-D visualization
        %   plotPartialDependence(MODEL,VARS,DATA) is for a fitted regression
        %   CompactTreeBagger. It plots the partial dependence PD between the
        %   predictor variables listed in VARS and the response of the model, by
        %   averaging over the data in the matrix or table DATA. You can specify
        %   one or two variables using VARS. VARS is either column indices of
        %   predictors or names of predictors. The names must match the entries in
        %   MODEL.PredictorNames.
        %
        %   plotPartialDependence(MODEL,VARS,LABELS,DATA) is for a fitted
        %   classification CompactTreeBagger. You can specify one or multiple
        %   classes using LABELS. The values in LABELS must match and have the same
        %   data type as the class names in MODEL.ClassNames. The plotted value PD
        %   is the partial dependence between the predictor variables listed in
        %   VARS and the scores for the classes specified in LABELS. If there are
        %   multiple values in VARS, then LABELS must be a single label.
        % 
        %   plotPartialDependence(...,Name,Value) specifies additional options using
        %   one or more name-value pair arguments:
        % 
        %   'Conditional'               'none' (default) instructs the method to create
        %                               a partial dependence plot with no conditioning,
        %                               'absolute' to create an individual conditional
        %                               expectation (ICE) plot, and 'centered' to
        %                               create an ICE plot with centered data.
        % 
        %   'NumObservationsToSample'   An integer specifying the number
        %                               of rows to sample at random from the dataset
        %                               (either the data specified by the argument DATA
        %                               or the data used to train MODEL). The default
        %                               is to use all rows.
        % 
        %   'QueryPoints'               The points at which to calculate
        %                               the partial dependence. For the case of a
        %                               single predictor, the value of 'QueryPoints'
        %                               must be a column vector. For the case of two
        %                               predictors, the value of 'QueryPoints' must be
        %                               a 1x2 cell array containing a separate vector
        %                               for each predictor. The default is to use 100
        %                               points equally spaced across the range of the
        %                               predictor.
        % 
        %   'UseParallel'               Can be true to instruct the method to perform
        %                               the averaging calculations in parallel (using
        %                               parfor), or false (default) to instruct the
        %                               method not to use parallel computations.
        % 
        % 	'Parent'                    Instructs the method to create a partial
        %                               dependence plot using the axes with handle
        %                               specified by the value of this parameter.
        % 
        %   AX = plotPartialDependence(...) returns a handle AX to the axes of the
        %   plot.
        %  
        %   Examples:
        %       % Partial Dependence Plot for bagged decision tree
        %       load census1994
        %       X = adultdata(:,{'age','workClass','education_num','marital_status','race',...
        %           'sex','capital_gain','capital_loss','hours_per_week','salary'});
        %       f = TreeBagger(10,X,'salary');
        %       f = compact(f); 
        % 
        %       plotPartialDependence(f,'age','<=50K',X);
        %       plotPartialDependence(f,{'age','workClass'},'>50K',X);
        %       plotPartialDependence(f,'age',{'<=50K' '>50K'},X)
        % 
        %       % Obtain optional output Axes handle
        %       ax = plotPartialDependence(f,'age','<=50K',X);
        % 
        %       % With optional name-value pairs
        %       plotPartialDependence(f,1,'<=50K',X,'NumObservationsToSample',100);
        %       plotPartialDependence(f,1,'<=50K',X,'UseParallel',true);
        %       plotPartialDependence(f,1,'<=50K',X,'UseParallel',true,'Conditional','none');
        % 
        %       % Plot the Individual Conditional Expectation
        %       plotPartialDependence(f,1,'<=50K',X,'Conditional','absolute');
        % 
        %       % Provide alternative query points
        %       xi = linspace(min(X{:,'age'}),max(X{:,'age'}))';
        %       plotPartialDependence(f,1,'<=50K',X,'QueryPoints',xi);
        % 
        %       % Provide alternative query points
        %       xi = linspace(min(X{:,'age'}),max(X{:,'age'}))';
        %       plotPartialDependence(f,1,'<=50K',X,'QueryPoints',xi);
        % 
        %       xi = cell(1,2);
        %       xi{1} = linspace(min(X{:,1}),max(X{:,1}))';
        %       xi{2} = linspace(min(X{:,9}),max(X{:,9}))';
        %       plotPartialDependence(f,[1,9],'<=50K',X,'QueryPoints',xi);
        
        %-------Check number of inputs----
        narginchk(3,14);
        features = convertStringsToChars(features);
        [varargin{:}] = convertStringsToChars(varargin{:});
        
        % Regression tree bagger
        if strcmp(bagger.Method,'regression')
            % class label is not an input for regression models
            % For a compact regression model the third argument is Data
            X = varargin{1};
            % Exclude Data from varargin 
            varargin(1) = [];
            % Classification flag is false and class label is empty
            classif = false;
            labels = [];
            
        % Classification tree bagger
        else
            % class label is an input for regression models
            if nargin < 4
                error(message('MATLAB:narginchk:notEnoughInputs'));
            else
            labels = varargin{1}; 
            X = varargin{2};
            % Exclude Data from varargin 
            varargin(1:2) = [];
            % classification flag is true
            classif = true;
            end
        end
           
        % Call the function from classreg.learning.impl
        ax = classreg.learning.impl.plotPartialDependence(bagger,features,...
            labels,X,classif,varargin{:});
        if(nargout > 0)
            AX = ax;
        end
        end
    end  
    
    methods(Hidden=true,Static=true)
        function a = empty(varargin),      throwUndefinedError(); end %#ok<STOUT>
    end
    
    
    methods(Hidden=true)
        function bagger = CompactTreeBagger(trees,classNames,varNames)
            bagger.Trees = trees;
            bagger.NTrees = length(trees);
            bagger.ClassNames = classNames;
            if isempty(varNames)
                error(message('stats:CompactTreeBagger:CompactTreeBagger:InvalidInput'));
            end
            bagger.VarNames = varNames;
        end
        
        function bagger = addTrees(bagger,trees)
            if ~iscell(trees)
                error(message('stats:CompactTreeBagger:addTrees:InvalidInput'));
            end
            nTrees = length(trees);
            bagger.Trees(end+1:end+nTrees) = trees;
            bagger.NTrees = bagger.NTrees + nTrees;
        end
        
        function disp(obj)
            internal.stats.displayClassName(obj);

            % Display body
            fprintf(1,'%s\n',getString(message('stats:TreeBagger:DispHeader',obj.NTrees)));
            fprintf(1,'%20s: %20s\n','Method',obj.Method);
            fprintf(1,'%20s: %20i\n','NumPredictors',length(obj.VarNames));
            if obj.Method(1)=='c'
                fprintf(1,'%20s:','ClassNames');
                for i=1:length(obj.ClassNames)
                    fprintf(1,' %s',['''' obj.ClassNames{i} '''']);
                end
                fprintf(1,'\n');
            end
            
            internal.stats.displayMethodsProperties(obj);
        end
        
        function [scores,nodes,labels] = treeEval(bagger,treeInd,x,doclassregtree)
            % Get the tree and classes
            tree = bagger.Trees{treeInd};
            if doclassregtree
                cTreeNames = tree.classname;
            else
                if strcmp(bagger.Method(1),'c')
                    cTreeNames = tree.ClassNames;
                else
                    cTreeNames = {};
                end
            end
            Nclasses = length(bagger.ClassNames);
            
            % Empty data?
            if isempty(x)
                scores = NaN(0,max(Nclasses,1));
                nodes = zeros(0,1);
                if Nclasses==0
                    labels = scores;
                else
                    labels = repmat(bagger.ClassNames{1},0,1);
                end
                return;
            end
                
            % Compute responses
            if Nclasses==0
                % For regression, get Yfit values
                if doclassregtree
                    [scores,nodes] = eval(tree,x); %#ok<EVLC>
                else
                    [scores,nodes] = predict(tree,x);
                end
                labels = scores;
            else
                % For classification, get class probabilities
                if doclassregtree
                    [labels,nodes] = eval(tree,x); %#ok<EVLC>
                    unmapped = classprob(tree,nodes);
                else
                    [labels,~,nodes] = predict(tree,x);
                    unmapped = tree.ClassProb(nodes,:);
                end
                                
                % Map classregtree classes onto bagger classes
                cFullNames = bagger.ClassNames;
                nFullClasses = length(cFullNames);
                N = size(x,1);
                scores = zeros(N,nFullClasses);
                [~,pos] = ismember(cTreeNames,cFullNames);
                scores(:,pos) = unmapped;
            end
        end
        
        function prox = flatprox(bagger,X,varargin)
            % Process inputs
            args = {'trees' {'nprint','numprint'}};
            defs = {  'all'        0};
            [useTrees,nprint] = internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check print-out frequency
            if isempty(nprint) || ~isnumeric(nprint) || numel(nprint)~=1 || nprint<0
                warning(message('stats:CompactTreeBagger:flatprox:BadPrintoutFreq'));
            end
            
            % How many trees to use?
            nHaveTrees = bagger.NTrees;
            if strcmpi(useTrees,'all')
                useTrees = 1:nHaveTrees;
            end
            
            % Any trees specified?
            if isempty(useTrees)
                error(message('stats:CompactTreeBagger:flatprox:EmptyUseTrees'));
            end
            
            % Correct type?
            if ~isnumeric(useTrees)
                error(message('stats:CompactTreeBagger:flatprox:UseTreesNotNumeric'));
            end
                
            % Invalid tree indices?
            if any(useTrees<1 | useTrees>nHaveTrees)
                error(message('stats:CompactTreeBagger:flatprox:BadUseTreesIndices', nHaveTrees));
            end
            
            % How many trees?
            nUseTrees = length(useTrees);

            % Init
            N = size(X,1);
            prox = zeros(1,N*(N-1)/2);
            
            % Is this old classregtree object?
            doclassregtree = isa(bagger.Trees{1},'classregtree');
            
            % Loop over trees
            for it=1:nUseTrees
                % Get the tree index
                itb = useTrees(it);
                
                % Get leaves on which training instances fall
                [~,nodes] = treeEval(bagger,itb,X,doclassregtree);
                
                % Get 1-proximities for this tree
                thisprox = pdist(nodes,'hamming');
                
                % Update proximities
                prox = prox + thisprox;
                
                % Report progress
                if nprint>0 && floor(it/nprint)*nprint==it
                    fprintf(1,'%s\n',getString(message('stats:TreeBagger:TreesDone',it)));
                end
            end
            
            % Normalize by the number of trees
            prox = 1 - prox/nUseTrees;
        end
    end
    
    
    methods(Access=protected)
        function [useTrees,treeWeights,useIforT] = ...
                checkTreeArgs(bagger,N,useTrees,treeWeights,useIforT)
            % How many trees to use?
            nHaveTrees = bagger.NTrees;
            if strcmpi(useTrees,'all')
                useTrees = 1:nHaveTrees;
            end
            
            % Any trees specified?
            if isempty(useTrees)
                error(message('stats:CompactTreeBagger:checkTreeArgs:EmptyUseTrees'));
            end
            
            % Correct type?
            if ~isnumeric(useTrees)
                error(message('stats:CompactTreeBagger:checkTreeArgs:UseTreesNotNumeric'));
            end
                
            % Invalid tree indices?
            if any(useTrees<1 | useTrees>nHaveTrees)
                error(message('stats:CompactTreeBagger:checkTreeArgs:UseTreesOutOfRange', nHaveTrees));
            end
            
            % How many trees?
            nUseTrees = length(useTrees);
                
            % Check tree weights
            if isempty(treeWeights)
                treeWeights = ones(nUseTrees,1);
            else
                if length(treeWeights)~=nUseTrees
                    error(message('stats:CompactTreeBagger:checkTreeArgs:TreeWeightsIndicesSizeMismatch'));
                end
            end
            
            % Are tree-instance indices supplied?
            if strcmpi(useIforT,'all')
                useIforT = true(N,nHaveTrees);
            else
                % Filled?
                if isempty(useIforT)
                    error(message('stats:CompactTreeBagger:checkTreeArgs:EmptyUseIforT'));
                end
               
                % Logical?
                if ~islogical(useIforT)
                    error(message('stats:CompactTreeBagger:checkTreeArgs:NotLogicalUseIforT'));
                end
                
                % Check size
                [n1,n2] = size(useIforT);
                if n1~=N || n2~=nHaveTrees
                    error(message('stats:CompactTreeBagger:checkTreeArgs:UseIforTSizeMismatch'));
                end
            end
        end    end
    
    
    methods(Access=public,Hidden=true)
        function [sc,eigen] = mdsProx(bagger,data,varargin)
            % Process inputs
            args = {'data' 'colors' 'labels' {'mdscoords','mdscoordinates'}};
            defs = {'predictors'       []       []       [1 2]};
            [datatype,plotColor,Y,axes] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check status and inputs
            if ~isempty(plotColor) && ~ischar(plotColor)
                error(message('stats:CompactTreeBagger:mdsProx:BadColor'));
            end
            
            % Check the type of input data
            if ~ischar(datatype) || ...
                    ~any(strncmpi(datatype,{'predictors' 'proximity'},length(datatype)))
                error(message('stats:CompactTreeBagger:mdsProx:BadDataArg'));
            end
            
            % Get proximity and predictors
            prox = [];
            X = [];
            if strncmpi(datatype,'predictors',length(datatype))
                X = data;
                N = size(X,1);
            end
            if strncmpi(datatype,'proximity',length(datatype))
                prox = data;
                [n1,n2] = size(prox);
                if n1~=n2
                    error(message('stats:CompactTreeBagger:mdsProx:InputSizeMismatch'));
                end
                N = n1;
            end
            
            % Check input labels
            if islogical(Y)
                Y = double(Y);
            end
            if bagger.Method(1)=='c' && ~isempty(Y)
                if ~iscellstr(Y)
                    Y = cellstr(classreg.learning.internal.ClassLabel(Y));
                end
                if length(Y)~=N
                    error(message('stats:CompactTreeBagger:mdsProx:XYSizeMismatch'));
                end
            end
            
            % Check axes
            if ~isnumeric(axes) || any(axes<1)
                error(message('stats:CompactTreeBagger:mdsProx:AxesNotIntegers'));
            end
            if numel(axes)<2 || numel(axes)>3
                error(message('stats:CompactTreeBagger:mdsProx:BadNumAxes'));
            end
            
            % Compute proximities if they were not supplied
            if isempty(prox)
                prox = proximity(bagger,X);
            end
            
            % Apply classical multi-D scaling
            dist = 1 - prox;
            [sc,eigen] = cmdscale(dist);
            
            % Check axes again
            if any(axes>size(sc,2))
                error(message('stats:CompactTreeBagger:mdsProx:BadAxisIndex'));
            end
            
            % Plot the scaling coordinates
            if ~isempty(plotColor)
                if bagger.Method(1)=='c' && ~isempty(Y)
                    [tf,truepos] = ismember(Y,bagger.ClassNames);
                    if any(~tf)
                        error(message('stats:CompactTreeBagger:mdsProx:BadLabels'));
                    end
                    nclass = max(truepos);
                    nclass = min(nclass,length(plotColor));
                    C = false(N,nclass);
                    for c=1:nclass
                        C(:,c) = truepos==c;
                    end
                else
                    C = true(N,1);
                    nclass = 1;
                end
                if numel(axes)==2
                    for c=1:nclass
                        scatter(sc(C(:,c),axes(1)),sc(C(:,c),axes(2)),plotColor(c));
                        hold on;
                    end
                else
                    for c=1:nclass
                        scatter3(sc(C(:,c),axes(1)),sc(C(:,c),axes(2)),...
                            sc(C(:,c),axes(3)),plotColor(c));
                        hold on;
                    end
                end
                if nclass>1
                    legend(bagger.ClassNames(1:nclass));
                end
                hold off;
            end
        end
                
        function [labels,scores,stdevs,scoreWeights] = ...
                predictAccum(bagger,X,varargin)
            % Process inputs
            args = {{'useifort' 'useinstancefortree'} 'trees' 'scores' 'stdevs' ...
                'scoreweights' 'treeweights'};
            defs = {                            'all'   'all'       []       [] ...
                            []            []};
            [useIforT,useTrees,scores,stdevs,scoreWeights,treeWeights] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check input data
            if ~istable(X)
                if ~isnumeric(X)
                    error(message('stats:CompactTreeBagger:predictAccum:XNotNumeric'));
                end
                if size(X,2)~=length(bagger.VarNames)
                    error(message('stats:CompactTreeBagger:predictAccum:XVarNamesSizeMismatch'));
                end
            end
            
            % Get data size.
            N = size(X,1);
            
            % Check score weights.
            if isempty(scoreWeights)
                scoreWeights = zeros(N,1);
            end
            if ~isnumeric(scoreWeights)
                error(message('stats:CompactTreeBagger:predictAccum:ScoreWeightsNotNumeric'));
            end
            if N~=length(scoreWeights)
                error(message('stats:CompactTreeBagger:predictAccum:ScoreWeightsSizeMismatch'));
            end
            
            % How many classes and columns in the score matrix?
            Nclasses = length(bagger.ClassNames);
            nCols = max(Nclasses,1);
                        
            % Check scores.
            if isempty(scores)
                scores = NaN(N,nCols);
            end
            if ~isnumeric(scores)
                error(message('stats:CompactTreeBagger:predictAccum:ScoresNotNumeric'));
            end
            if size(scores,1)~=N || size(scores,2)~=nCols
                error(message('stats:CompactTreeBagger:predictAccum:ScoresSizeMismatch'));
            end

            % Check standard deviations.
            if isempty(stdevs)
                stdevs = NaN(N,nCols);
            end
            if ~isnumeric(stdevs)
                error(message('stats:CompactTreeBagger:predictAccum:StdevsNotNumeric'));
            end
            if size(stdevs,1)~=N || size(stdevs,2)~=nCols
                error(message('stats:CompactTreeBagger:predictAccum:StdevsSizeMismatch'));
            end

            % Check input tree arguments
            [useTrees,treeWeights,useIforT] = ...
                checkTreeArgs(bagger,N,useTrees,treeWeights,useIforT);
            
            % Is this an old classregtree object?
            doclassregtree = isa(bagger.Trees{1},'classregtree');
            
            % Update mean scores and cumulative standard deviations looping
            % over trees
            nUseTrees = length(useTrees);
            for iuse=1:nUseTrees
                % Get the tree index in the full list
                it = useTrees(iuse);
                
                % Which observations can be processed by this tree?
                tf = useIforT(:,it);

                % If there are NaN rows in input scores, treat them as
                % "no observations" and reset to 0
                nanscore = all(isnan(scores),2) & tf;
                scores(nanscore,:) = 0;
                if nargout>2
                    stdevs(nanscore,:) = 0;
                end
            
                % Get data, weights and scores for this tree
                thisX = X(tf,:);
                thisW = treeWeights(iuse);
                thisR = treeEval(bagger,it,thisX,doclassregtree);
                
                % Get new weights
                newScoreWeights = scoreWeights(tf) + thisW;
                
                % Update means
                delta = thisR - scores(tf,:);
                gamma = bsxfun(@rdivide,delta*thisW,newScoreWeights);
                scores(tf,:) = scores(tf,:) + gamma;
                
                % Update standard deviations
                if nargout>2
                    stdevs(tf,:) = stdevs(tf,:) + ...
                        bsxfun(@times,delta.*gamma,scoreWeights(tf));
                end
                
                % Update weights
                scoreWeights(tf) = newScoreWeights;
            end
            
            % Find class with max probability for classification
            isNaN = all(isnan(scores),2);
            if ~all(isnan(bagger.DefaultScore))
                scores(isNaN,:) = repmat(bagger.DefaultScore,sum(isNaN),1);
            end
            if bagger.Method(1)=='c' % classification
                % Init
                labels = cell(N,1);
                
                % Check if default scores are NaN's
                if all(isnan(bagger.DefaultScore))
                    notNaN = ~isNaN;
                    labels(:) = cellstr('');
                else
                    notNaN = (1:N)';
                end
                
                % Find class with max prob
                [~,classNum] = max(scores(notNaN,:),[],2);
                labels(notNaN) = bagger.ClassNames(classNum);
            else % regression
                labels = scores;
            end
        end
        
        function [tau,trainNode,testNode,TW] = ...
                quantileNode(bagger,X,Xtrain,wtrain,varargin)
        %QUANTILENODE Quantile node assignment
        %   [TAU,TRAINNODE,TESTNODE,TW] = QUANTILENODE(BAGGER,X,XTRAIN,WTRAIN)
        %       Inputs:
        %           BAGGER      - Object of type CompactTreeBagger.
        %           X           - Test data.
        %           XTRAIN      - Training data used to train this ensemble.
        %           WTRAIN      - Observation weights for the training data.
        %       Outputs:
        %           TAU         - Vector of desired quantiles. 0.5 (scalar) if no
        %                         values are passed.
        %           TRAINNODE   - An Ntrain-by-T array for Ntrain training
        %                         observations and T trees. Element (I,J) of this
        %                         array is the leaf node number for J-th tree on
        %                         which I-th observation landed. If I-th
        %                         observation was not used for training J-th tree
        %                         (was out of bag), this element is NaN. See
        %                         the 'InBagIndices' argument.
        %           TESTNODE    - An N-by-T array for N test observations. Element
        %                         (I,J) of this array is the leaf node number for
        %                         J-th tree on which I-th observation landed. If
        %                         I-th observation should not be used for J-th
        %                         tree, this elements is NaN. See the 'UseIforT'
        %                         argument.
        %           TW          - An Ntrain-by-T array for Ntrain training
        %                         observations and T trees. Element (I,J) of this
        %                         array is the weight of observation I for tree J.
        %                         This weight is computed by taking the weight for
        %                         observation I passed into TreeBagger and
        %                         multiplying it by the number of times observation
        %                         I was used for tree J. If observation I was not
        %                         used for training tree J, this element is zero.
        %
        %       Optional arguments:
        %           'Quantiles'      - Vector with elements between 0 and 1.
        %           'UseIforT'       - String 'all' or logical array of size N-by-T
        %                              for N test observations and T trees.
        %           'InBagIndices'   - Array of size Ntrain-by-T for Ntrain
        %                              training observations in Xtrain and T trees.
        %                              Column J of this array are indices of
        %                              observations in Xtrain used for training
        %                              tree J.
        %           'Trees'          - String 'all' or vector of tree indices to be
        %                              used.
        %           'TreeWeights'    - Vector of tree weights. If empty, all trees
        %                              get equal weights.

            % You can get quantile estimates only for regression.
            if ~strcmpi(bagger.Method,'regression')
                error(message('stats:CompactTreeBagger:quantileNode:UseForRegression'));
            end
            
            % Check input data
            if ~istable(X)
                if ~isnumeric(X)
                    error(message('stats:CompactTreeBagger:predictAccum:XNotNumeric'));
                end
                if size(X,2)~=length(bagger.VarNames)
                    error(message('stats:CompactTreeBagger:predictAccum:XVarNamesSizeMismatch'));
                end
            end
            
            % Process args
            args = {'quantiles' {'useifort' 'useinstancefortree'} 'inbagindices' ...
                'trees' 'treeweights'};
            defs = {        0.5                             'all'             [] ...
                  'all'            []};
            [tau,useIforT,ibIdx,useTrees,treeWeights] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
                        
            % Check quantiles
            internal.stats.checkSupportedNumeric('Quantiles',tau);
            if ~isvector(tau) || any(tau<0) || any(tau>1)
                error(message('stats:CompactTreeBagger:quantileNode:BadQuantiles'));
            end

            % Get data size
            N = size(X,1);
            
            % Check input tree arguments
            [useTrees,treeWeights,useIforT] = ...
                checkTreeArgs(bagger,N,useTrees,treeWeights,useIforT);
            
            trees = bagger.Trees(useTrees);
            nHaveTrees = bagger.NTrees;
                        
            % In-bag indices must be passed.
            if isempty(ibIdx)
                error(message('stats:CompactTreeBagger:quantileNode:NoInBagIdx'));
            end            
            if size(ibIdx,2)~=nHaveTrees
                error(message('stats:CompactTreeBagger:quantileNode:BadIdxSize',...
                    nHaveTrees));
            end
            
            useIforT = useIforT(:,useTrees);
            ibIdx = ibIdx(:,useTrees);
            
            [trainNode,TW] = localGetTrainingNodes(trees,Xtrain,wtrain,ibIdx);
                        
            TW = bsxfun(@times,TW,treeWeights(:)');
            
            testNode = localPredictNode(trees,X,useIforT);
        end
        
        function [Yhat,YW] = quantilePredictCompact(bagger,X,Xtrain,Ytrain,Wtrain,varargin)
            %QUANTILEPREDICTCOMPACT Quantile prediction by bagged decision trees.
            %   YHAT=QUANTILEPREDICTCOMPACT(ENS,X,XTRAIN,YTRAIN,WTRAIN) returns
            %   expected median of the response distribution for an ensemble of bagged
            %   regression decision trees ENS and predictor matrix X.
            %   QUANTILEPREDICTCOMPACT assumes that the ensemble has been trained using
            %   predictor matrix XTRAIN, response vector YTRAIN and observation weights
            %   WTRAIN. If you pass X as an N-by-P matrix for N observations and P
            %   predictors, YHAT is an N-by-1 vector of conditional medians at the
            %   respective rows of X.
            %
            %   [YHAT,YW]=QUANTILEPREDICTCOMPACT(ENS,X,XTRAIN,YTRAIN,WTRAIN) also
            %   returns an M-by-N sparse matrix YW specifying weights of M response
            %   values YTRAIN for N test observations in X.
            %
            %   [...]=QUANTILEPREDICTCOMPACT(ENS,X,XTRAIN,YTRAIN,WTRAIN,PARAM1',val1,'PARAM2',val2,...)
            %   specifies optional parameter name/value pairs:
            %       'InBagIndices' - Numeric matrix of size Nobs-by-NTrees for Nobs
            %                        observations in the training data and NTrees trees
            %                        in the ensemble holding integer indices of in-bag
            %                        observations for trees.
            %       'Quantile'  - Vector of quantiles between 0 and 1. If you pass a
            %                     vector of Q quantiles, QUANTILEPREDICTCOMPACT returns an
            %                     N-by-Q matrix YHAT.
            %       'Trees'     - Array of tree indices to use for computation of
            %                     responses.  Default is 'all'.
            %       'TreeWeights' - Array of NTrees weights for weighting votes from the
            %                     specified trees.
            %       'UseIforT'  - Logical matrix of size Nobs-by-NTrees indicating
            %                     which trees should be used to make predictions for
            %                     each observation.  By default all trees are used for
            %                     all observations. If no trees are used for an
            %                     observation (the row of the matrix has no true
            %                     values), QUANTILEPREDICTCOMPACT returns the respective
            %                     quantile for ENS.Y.

            [tau,trainNode,testNode,TW] = ...
                quantileNode(bagger,X,Xtrain,Wtrain,varargin{:});

            YW = [];
            doyw = nargout > 1;
            if doyw
                [Yhat,YW] = localPredictResponse('ensemble',[],...
                    tau,Ytrain,trainNode,testNode,TW);
            else
                Yhat = localPredictResponse('ensemble',[],...
                    tau,Ytrain,trainNode,testNode,TW);
            end
            
            [Yhat,YW] = qrfixnans(Yhat,YW,Ytrain,tau,Wtrain);
        end
        
        function err = quantileErrorCompact(bagger,X,Y,Xtrain,Ytrain,Wtrain,varargin)            
            %QUANTILEERRORCOMPACT Quantile regression error by bagged decision trees.
            %   ERR=QUANTILEERRORCOMPACT(ENS,X,Y,XTRAIN,YTRAIN,WTRAIN) returns
            %   quantile prediction error for a TreeBagger ensemble of bagged
            %   regression decision trees ENS, predictor matrix X and response vector
            %   Y. QUANTILEERRORCOMPACT assumes that the ensemble has been trained
            %   using predictor matrix XTRAIN, response vector YTRAIN and observation
            %   weights WTRAIN. If you pass X as an N-by-P matrix for N observations
            %   and P predictors, YHAT is an N-by-1 vector of conditional medians at
            %   the respective rows of X.
            %
            %   QUANTILEERRORCOMPACT accepts the same optional arguments as
            %   QUANTILEPREDICTCOMPACT.

            args = {    'mode' 'weights'};
            defs = {'ensemble'        []};
            [mode,w,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Process mode
            if ~ischar(mode) || ...
                    ~any(strncmpi(mode,{'cumulative' 'individual' 'ensemble'},length(mode)))
                error(message('stats:CompactTreeBagger:error:BadMode'));
            end
            
            [tau,trainNode,testNode,TW] = ...
                quantileNode(bagger,X,Xtrain,Wtrain,extraArgs{:});

            Q = numel(tau);
            [N,T] = size(testNode);
            
            % Process weights
            if isempty(w)
                w = ones(N,1)/N;
            end
            if ~isnumeric(w) || length(w)~=N || any(w<0) || all(w==0)
                error(message('stats:CompactTreeBagger:predictAccum:BadW', N));
            end
            
            if     any(strncmpi(mode,{'individual' 'cumulative'},length(mode)))
                err = zeros(T,Q);
                for t=1:T
                    Yhat = localPredictResponse(...
                        mode,t,tau,Ytrain,trainNode,testNode,TW);
                    Yhat = qrfixnans(Yhat,[],Ytrain,tau,Wtrain);
                    err(t,:) = qrloss(Yhat,Y,tau,w);
                end
            elseif strncmpi(mode,'ensemble',length(mode))
                Yhat = localPredictResponse(...
                    mode,[],tau,Ytrain,trainNode,testNode,TW);
                Yhat = qrfixnans(Yhat,[],Ytrain,tau,Wtrain);
                err = qrloss(Yhat,Y,tau,w);
            end
        end
        
        function [Y,args] = inferResponse(bagger,X,varargin)
            % Infer whether the first element of varargin should be taken as the
            % response, or if the response appears as part of the table X.
            if isa(bagger.Trees{1},'classregtree')
                % Response must be part of varargin
                [Y,args] = classreg.learning.internal.inferResponse('',[],varargin{:});
            else
                % Response may be part of varargin or may be taken from table
                rname = bagger.Trees{1}.ResponseName;
                [Y,args] = classreg.learning.internal.inferResponse(rname,X,varargin{:});
            end
            
            % Char Y is not accepted by TreeBagger
            if size(Y,1)~=size(X,1)
                error(message('stats:CompactTreeBagger:inferResponse:BadY',size(X,1)));
            end
        end        
    end
    
    
    %%%Methods for tall array
      methods(Access=protected)
             
         function [Y,args] = inferResponseTall(this,X,varargin)
        %    This utility function replicates the same behavior as
        %    TreeBagger/inferResponse but for tall inputs.
             if isa(this.Trees{1},'classregtree')
                [Y,args] = classreg.learning.internal.inferResponse('',[],varargin{:});
             else
                [Y,args] = internal.stats.bigdata.ClassRegModelTallAdapter.inferResponse(...
                    this.Trees{1}.ResponseName,X,varargin{:});
             end
         end %method inferResponse
         
         function mar= meanMarginTall(this,X,varargin)
            %MEANMARGIN average classification margin for TreeBagger TALL
            %  Replicates the same behavior as TreeBagger.meanMargin 
            %  but for tall inputs.  
             if ~istall(X)
                error(message('stats:tall:classreg:PredictorsNotTall'))
             end            
         
            [Y,varargin] = inferResponseTall(this,X,varargin{:});
            if ~istall(Y)
                error(message('stats:tall:classreg:ResponseNotTall'))
            end
            
            %disallow trees and treeweights
            args = { 'trees' 'treeweights' 'weights'};
            defs = {     []            []   []};
            [useTrees,treeWeights, W, ~, extraArgs] = ...
            internal.stats.parseArgs(args,defs,varargin{:});
            if ~isempty(useTrees)||~isempty(treeWeights) || ~gather(isempty(W))
                     error(message('stats:tall:TreeBagger:NonSupportOpt'))
            end
   
%             args = {'weights'};
%             defs = {       []};
%             [W,os,extraArgs] = ...
%                 internal.stats.parseArgs(args,defs,extraArgs{:});
%             
%             if strcmp(tall.getClass(X),'table') && internal.stats.isString(W)
%                 W = X.(W); 
%             elseif isempty(W) %&&~istall(W) 
%                 os.weights = false; 
%                 W =hSlicefun(@(x) ones(size(x,1),1),Y);
%             elseif ~istall(W) 
%                 error(message('stats:tall:classreg:WeightsNotTall'))
%             end
           
            % Compute mean margins
            m0 = margin(this,X,Y,extraArgs{:});
                      
            W =hSlicefun(@(x) ones(size(x,1),1),Y);
            %If there is NaN values in m0
            sumW = sum(~isnan(m0).*W,1,'omitnan');  % sum of weights with NaNs treated independently
            weightedSum_m0 = sum(m0.*W,1,'omitnan');        % weighted mean times sum of weights
            mar = weightedSum_m0 ./sumW;
%             
%            if os.weights
%              tf =  size(X,1)==size(W,1) & size(W,2)==1 & all(W>=0 | isnan(W)); 
%             %                                                        ^
%             % Note: non-tall behavior propagates NaN's so that the mean
%             % margin will be all NaNs if there is NaN in W.

%             % Checks on weights lazily added to the code branch that returns mar:
%              [mar,~] = lazyValidate(mar,tf,{@(e,tf) tf,'stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:BadW'});
%            end
          
         end %meanMarginTall
         
         function err = errorTall(this,X,varargin)
            if ~istall(X)
                error(message('stats:tall:classreg:PredictorsNotTall'))
             end            
         
            [Y,varargin] = inferResponseTall(this,X,varargin{:});
            if ~istall(Y)
                error(message('stats:tall:classreg:ResponseNotTall'))
            end
            
            args = {'mode'          {'useifort' 'useinstancefortree'} 'trees' 'treeweights' 'weights'};
            defs = {'cumulative'                            []    []            []     []};
            [mode,useIforT,useTrees,treeWeights,W] = ...
            internal.stats.parseArgs(args,defs,varargin{:});
            if ~isempty(useIforT)||~isempty(useTrees)||~isempty(treeWeights) 
                     error(message('stats:tall:TreeBagger:NonSupportOpt'))
            end
            
            if ~istall(W) && isempty(W)
                W = hSlicefun(@(x) ones(size(x,1),1),Y);
            elseif istall(W)
                [X,W] = validateSameTallSize(X,W);
            elseif ~isscalar(W) 
                error(message('stats:tall:classreg:WeightsNotTall'))
            else
                W = hSlicefun(@(x) ones(size(x,1),1).*W,Y);
            end
           
%             if os.weights
%                  tf =  size(X,1)==size(W,1) & size(W,2)==1 & all(W>=0);%| isnan(W)); 
%                  
%                 %                                                        ^
%                 % Note: non-tall behavior propagates NaN's so that the mean
%                 % margin will be all NaNs if there is NaN in W. We don't add check
%                 % here. Instead, we ingore the observations with NaN 
%                 
%                 % Checks on weight lazily added to the code branch that returns err:
%                [err,~] = lazyValidate(err,tf,{@(e,tf) tf,'stats:classreg:learning:FullClassificationRegressionModel:prepareDataCR:BadW'});
%            
%              end
%        
             % Based on comments above, remove observations with NaNs or zeros in W
             lidx_W = ~(isnan(W)|W==0);
             [X,W] = hFilterslices(lidx_W,X,W);

             [err,~] = hAggregatefun(@(chunkX,chunkY,W)inmemory_error(this,chunkX,chunkY,W,mode),...
                          @reduceError,X,Y,W);
             err=gather(err)';
         end
        
%         function adapter = makeTreeBaggerAdapter(obj,varargin)
%         %Decide not to have a TreeBaggerAdapter to avoid make another 
          %copy of TreeBagger model that could takes a lot memory.
%           if  any(cellfun(@(x) istall(x),varargin))
%              adapter = internal.stats.bigdata.TreeBaggerTallAdapter(obj);
%           else
%              adapter = [];
%           end
%  
%         end %makeTreeBaggerAdapter
    end
    %%%end Methods for tall array
end

function throwUndefinedError()
error(message('stats:CompactTreeBagger:UndefinedFunction'));
end

% Quantile response for one set of weights. These weights are determined by
% counting coincident nodes in the training and test data and cannot contain
% zeros.
function yhat = qrw(y,tau,w)
idx = w>0;
if any(idx)
    %yhat = ksdensity(y(idx),tau,'function','icdf','weights',w(idx));

    % Use linear interpolation because it is much faster than ksdensity.

    y = y(idx);
    w = w(idx);
    [y,sorted] = sort(y);
    w = w(sorted);
    
    wsum = sum(w);
    w = w/wsum;
    
    W = ( cumsum(w) + (1-flipud(cumsum(flipud(w)))) ) /2;
    
    yhat = interp1([0-eps;W;1+eps],[y(1);y;y(end)],tau);
else
    yhat = NaN(1,numel(tau));
end
end

% Loss for quantile regression. ERR is an 1-by-Q vector, where
% Q=numel(tau).
function err = qrloss(Yhat,y,tau,w)
tau = tau(:)';
DY = bsxfun(@minus,y,Yhat);
err = tau.*classreg.learning.internal.wnanmean(DY.*(DY>0),w) ...
    - (1-tau).*classreg.learning.internal.wnanmean(DY.*(DY<0),w);
end

% Replace NaN responses for observations without trees with responses
% obtained from the entire dataset. This function is used for quantile
% regression only.
function [Yhat,YW] = qrfixnans(Yhat,YW,y,tau,w)
badrows = all(isnan(Yhat),2);
if any(badrows)
    yhat = qrw(y,tau,w);
    nbad = sum(badrows);
    
    Yhat(badrows,:) = repmat(yhat(:)',nbad,1);
    
    if ~isempty(YW) && nargout>1
        YW(:,badrows) = repmat(w(:),1,nbad);
    end
end
end

function [node,TW] = localGetTrainingNodes(trees,X,w,ibIdx)

T = numel(trees);
N = size(X,1);

node = NaN(N,T);
TW = NaN(N,T);

for t=1:T
    tree = trees{t};
    
    idxtrain = ibIdx(:,t);
    nidx = accumarray(idxtrain,1);
    nidx(end+1:N) = 0;
    tw = w.*nidx;
    
    ibtf = tw>0;
    [~,tnode] = predict(tree,X(ibtf,:));
    
    S = sparse(1:sum(ibtf),tnode,tw(ibtf));
    sumPerNode = sum(S,1);
    goodNodes = sumPerNode>0;
    tw(ibtf) = full(sum(bsxfun(@rdivide,S(:,goodNodes),sumPerNode(goodNodes)),2));
    
    node(ibtf,t) = tnode;
    TW(:,t) = tw;
end
end


function node = localPredictNode(trees,X,useIforT)

T = numel(trees);
N = size(X,1);
node = NaN(N,T);

for t=1:T
    idx = useIforT(:,t);
    [~,node(idx,t)] = predict(trees{t},X(idx,:));
end
end

% Predict response for one of the 3 modes: 'individual', 'cumulative' or
% 'ensemble'. If mode is 'individual' or 'cumulative', predict for tree
% indexed by t (2nd input). If mode is 'ensemble', set t to []. Yhat is
% matrix of size Ntest-by-Q, where Q is the number of quantiles.
%
% The second output YW is a matrix of size Ntrain-by-Ntest specifying
% response distribution weights for each test observation.
function [Yhat,YW] = localPredictResponse(mode,t,tau,Ytrain,trainNode,testNode,TW)

N = size(testNode,1);

Yhat = NaN(N,numel(tau));

doyw = nargout > 1;
if doyw
    YW = sparse(numel(Ytrain),N);
end

% qw is a vector of weights with numel(Ytrain) elements specifying the
% distribution of response values.

if     strncmpi(mode,'individual',length(mode))
    maxnode = max(testNode(:,t));
    if isnan(maxnode)
        return
    end
    done = zeros(maxnode,1);
    for n=1:N
        idx = done(testNode(n,t));
        if idx > 0 % already did this, fetch results
            Yhat(n,:) = Yhat(idx,:);
            if doyw
                YW(:,n) = YW(:,idx); %#ok<SPRIX>
            end
        else   % need to compute this for first time
            qw = (trainNode(:,t)==testNode(n,t)).*TW(:,t);
            Yhat(n,:) = qrw(Ytrain,tau,qw);
            if doyw
                YW(:,n) = qw; %#ok<SPRIX>
            end
            done(testNode(n,t)) = n;
        end
    end
elseif strncmpi(mode,'cumulative',length(mode))
    for n=1:N
        qw = sum(bsxfun(@eq,trainNode(:,1:t),testNode(n,1:t)).*TW(:,1:t),2);
        Yhat(n,:) = qrw(Ytrain,tau,qw);
        if doyw
            YW(:,n) = qw; %#ok<SPRIX>
        end
    end
elseif strncmpi(mode,'ensemble',length(mode))
    for n=1:N
        qw = sum(bsxfun(@eq,trainNode,testNode(n,:)).*TW,2);
        Yhat(n,:) = qrw(Ytrain,tau,qw);
        if doyw
            YW(:,n) = qw; %#ok<SPRIX>
        end
    end
end

% Normalize YW
if doyw
    yw = sum(YW,1);
    good = yw>0;
    if any(good)
        YW(:,good) = bsxfun(@rdivide,YW(:,good),yw(good));
    end
end

end

function [err,sumW]=inmemory_error(bag,X,Y,Weight,mode)
if isempty(X)
    [err, sumW] = deal(zeros(0, bag.NumTrees));
else
  [err, sumW] = error(bag,X,Y,'weights',Weight,'mode',mode);
  err = err';
  sumW = sumW';
end
end

function [err,sumW]= reduceError(err,sumW)
if isempty(err)
    [err, sumW] = deal(zeros(0, 1));
    return;
else
   err  = sum(err.*sumW,1,'omitnan');
   sumW = sum(sumW,1,'omitnan');
   idx = sumW > 0; %filter the trees with zero total observation weights
   err = err(idx)./sumW(idx);
end
end

function m = localMargin(bagger,x,varargin)
if isempty(x)
    m = [];
else
    m = bagger.margin(x,varargin{:});
end
end
