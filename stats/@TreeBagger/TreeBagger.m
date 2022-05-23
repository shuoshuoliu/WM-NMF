classdef TreeBagger < classreg.learning.internal.DisallowVectorOps
%TREEBAGGER Bootstrap aggregation for an ensemble of decision trees.
%   TreeBagger bags an ensemble of decision trees for either classification
%   or regression. 'Bagging' stands for 'bootstrap aggregation'. Every tree
%   in the ensemble is grown on an independently-drawn bootstrap replica of
%   input data. Observations not included in this replica are "out of bag"
%   for this tree. To compute prediction of an ensemble of trees for unseen
%   data, TREEBAGGER takes an average of predictions from individual trees.
%   
%   To estimate the prediction error of the bagged ensemble, you can compute
%   predictions for each tree on its out-of-bag observations, average these
%   predictions over the entire ensemble for each observation and then
%   compare the predicted out-of-bag response with the true value at this
%   observation.
%
%   TreeBagger relies on functionality of ClassificationTree and
%   RegressionTree for growing individual trees. In particular,
%   ClassificationTree and RegressionTree accept the number of features
%   selected at random for each decision split as an optional input
%   argument.
%
%   The 'compact' method returns an object of another class,
%   CompactTreeBagger, with sufficient information to make predictions
%   using new data. This information includes the tree ensemble, variable
%   names, and class names (for classification).  CompactTreeBagger
%   requires less memory than TreeBagger, but only TreeBagger has methods
%   for growing more trees for the ensemble.  Once you grow an ensemble of
%   trees using TreeBagger, if you no longer need access to the training
%   data, you can opt to work with the compact version of the trained
%   ensemble from then on.
%
%   TreeBagger properties:
%      X                     - Matrix or table of predictors used to create the ensemble.
%      Y                     - Y data used to create the ensemble.
%      W                     - Weights of observations used to create the ensemble.
%      InBagFraction         - Fraction of in-bag observations.
%      SampleWithReplacement - Flag to sample with replacement.
%      TreeArguments  		 - Cell array of arguments for ClassificationTree or RegressionTree.
%      ComputeOOBPrediction  - Flag to compute out-of-bag predictions.
%      ComputeOOBPredictorImportance - Flag to compute out-of-bag predictor importance.
%      Prune          		 - Flag to prune trees.
%      MergeLeaves   		 - Flag to merge leaves that do not improve risk.
%      NumPredictorsToSample - Number of predictors for random feature selection.
%      MinLeafSize           - Minimum number of observations per tree leaf.
%      OOBIndices            - Indicator matrix for out-of-bag observations.
%      Trees                 - Decision trees in the ensemble.
%      NumTrees              - Number of decision trees in the ensemble.
%      ClassNames            - Names of classes.
%      Prior                 - Prior class probabilities.
%      Cost           		 - Misclassification costs.
%      PredictorNames		 - Predictor names.
%      Method         		 - Method used by trees (classification or regression).
%      OOBInstanceWeight     - Count of out-of-bag trees for each observation.
%      OOBPermutedPredictorDeltaError       - Predictor importance for classification error.
%      OOBPermutedPredictorDeltaMeanMargin  - Predictor importance for classification margin.
%      OOBPermutedPredictorCountRaiseMargin - Predictor importance for raising margin.
%      DeltaCriterionDecisionSplit          - Split criterion contributions for each predictor.
%      NumPredictorSplit      - Number of decision splits on each predictor.
%      SurrogateAssociation   - Predictor associations.
%      Proximity      		  - Proximity matrix for observations.
%      OutlierMeasure 		  - Measure for determining outliers.
%      DefaultYfit   		  - Default value returned by PREDICT and OOBPREDICT.
%
%   TreeBagger methods:
%      TreeBagger/TreeBagger - Create an ensemble of bagged decision trees.
%      append           - Append new trees to ensemble.
%      compact          - Compact ensemble of decision trees.
%      error            - Error (misclassification probability or MSE).
%      fillprox         - Fill proximity matrix for training data.
%      growTrees        - Train additional trees and add to ensemble.
%      margin           - Classification margin.
%      mdsprox          - Multidimensional scaling of proximity matrix.
%      meanMargin       - Mean classification margin per tree.
%      oobError         - Out-of-bag error.
%      oobMargin        - Out-of-bag margins.
%      oobMeanMargin    - Out-of-bag mean margins.
%      oobPredict       - Ensemble predictions for out-of-bag observations.
%      oobQuantileError - Out-of-bag quantile loss.
%      oobQuantilePredict - Predicted quantiles for out-of-bag observations.
%      predict          - Predict response.
%      quantileError    - Quantile loss.
%      quantilePredict  - Predicted quantiles.
%
%   Example:
%      % Using matrices
%      load fisheriris
%      b = TreeBagger(50,meas,species,'OOBPrediction','on')
%      plot(oobError(b))
%      xlabel('number of grown trees')
%      ylabel('out-of-bag classification error')
%
%      % Using a table
%      tbl = readtable('fisheriris.csv','format','%f%f%f%f%C');
%      b2 = TreeBagger(50,tbl,'Species')
% 
%    See also TreeBagger/TreeBagger, CompactTreeBagger, fitctree, fitrtree,
%    classreg.learning.classif.CompactClassificationTree,
%    classreg.learning.regr.CompactRegressionTree.
    
%   Copyright 2008-2020 The MathWorks, Inc.

    properties(SetAccess=protected,GetAccess=public,Dependent)
        %X X data used to create the ensemble.
        %   The X property contains the predictor values. It is a table for an
        %   ensemble trained on a table, or a numeric matrix for a model trained on
        %   a matrix. It has size Nobs-by-Nvars, where Nobs is the number of
        %   observations (rows) and Nvars is the number of predictor variables or
        %   columns in the training data.  This matrix contains the predictor (or
        %   feature) values.
        %
        %   See also TREEBAGGER.
        X;
    end
        
    properties(SetAccess=protected,GetAccess=public)
        %Y Y data used to create the ensemble.
        %   The Y property is an array of true class labels for classification, or
        %   response values for regression. Y is a numeric column vector for
        %   regression and a cell array of strings for classification.
        %
        %   See also TREEBAGGER.
        Y = [];
        
        %W Weights of observations used to create the ensemble.
        %   The W property is a numeric vector of size Nobs, where Nobs is the
        %   number of observations. These weights are used for growing every
        %   decision tree in the ensemble.
        %
        %   See also TREEBAGGER.
        W = [];
        
        %SAMPLEWITHREPLACEMENT Flag to sample with replacement.
        %   The SampleWithReplacement property is a logical flag specifying if data
        %   are sampled for each decision tree with replacement. True if data are
        %   sampled with replacement and false otherwise. True by default.
        %
        %   See also TREEBAGGER.
        SampleWithReplacement = true;
        
        %COMPUTEOOBPREDICTION Flag to compute out-of-bag predictions.
        %   The ComputeOOBPrediction property is a logical flag specifying whether out-of-bag
        %   predictions for training observations should be computed.  The default is
        %   false.
        %
        %   If this flag is true, the following properties are available:
        %       OOBIndices, OOBInstanceWeight.
        %   If this flag is true, the following methods can be called:
        %       oobError, oobMargin, oobMeanMargin.
        %
        %   See also TREEBAGGER, OOBERROR, OOBMARGIN, OOBMEANMARGIN, OOBINDICES,
        %   OOBINSTANCEWEIGHT.
        ComputeOOBPrediction = false;
        
        %PRUNE Flag to prune trees.
        %   The Prune property is true if decision trees are pruned and false if
        %   they are not. Pruning decision trees is not recommended for ensembles.
        %   The default value is false.
        %
        %   See also TREEBAGGER, ClassificationTree/prune, RegressionTree/prune.
        Prune = false;
        
        %MERGELEAVES Flag to merge leaves that do not improve risk.
        %   The MergeLeaves property is true if decision trees have their leaves
        %   with the same parent merged for splits that do not decrease the total
        %   risk, and false otherwise. The default value is false.
        %
        % See also TREEBAGGER.
        MergeLeaves = false;
        
        %OOBINDICES Indicator matrix for out-of-bag observations.
        %   The OOBIndices property is a logical array of size Nobs-by-NumTrees where
        %   Nobs is the number of observations in the training data and NumTrees is
        %   the number of trees in the ensemble.  The (I,J) element is true if
        %   observation I is out-of-bag for tree J and false otherwise.  In other
        %   words, a true value means observation I was not selected for the
        %   training data used to grow tree J.
        %
        %   See also TREEBAGGER.
        OOBIndices = [];
    end
    
    properties(SetAccess=public,GetAccess=public,Dependent=true,Hidden=true)
        DefaultScore;
    end
    
    properties(SetAccess=protected,GetAccess=public,Hidden=true)
        InBagIndices = [];
    end
    
    properties(SetAccess=protected,GetAccess=protected)
        Compact = [];
        PrivOOBPermutedVarDeltaError = [];% NumTrees-by-Nvars
        PrivOOBPermutedVarDeltaMeanMargin = [];
        PrivOOBPermutedVarCountRaiseMargin = [];
        PrivProx = [];% Proximity matrix in the form of a 1D array
        DataSummary = [];
        ClassSummary = [];
        ModelParams = [];
        PrivateX = [];
        Version = 1; % 1 for loaded objects, constructor sets higher value
    end
    
    properties(SetAccess=protected,GetAccess=public,Dependent=true)
        %INBAGFRACTION Fraction of in-bag observations.
        %   The InBagFraction property is the fraction of observations to be randomly
        %   selected with replacement for each bootstrap replica. The size of each
        %   replica is given by N*InBagFraction, where N is the number of observations in
        %   the training set.  The default value is 1.
        %
        %   See also TREEBAGGER.
        InBagFraction;
        
        %TREEARGUMENTS Cell array of arguments for FITCTREE or FITRTREE.
        %   The TreeArguments property is a cell array of arguments for the FITCTREE or
        %   FITRTREE functions. These arguments are used by TreeBagger in growing
        %   new trees for the ensemble.
        %
        %   See also TreeBagger, fitctree, fitrtree.
        TreeArguments;
        
        %COMPUTEOOBPREDICTORIMPORTANCE Flag to compute out-of-bag variable importance.
        %   The ComputeOOBPredictorImportance property is a logical flag specifying
        %   whether out-of-bag estimates of variable importance should be computed.
        %   The default is false. If this flag is true, COMPUTEOOBPREDICTION is
        %   true as well.
        %
        %   If this flag is true, the following properties are available:
        %      OOBPermutedPredictorDeltaError, OOBPermutedPredictorDeltaMeanMargin,
        %      OOBPermutedPredictorCountRaiseMargin
        %
        %   See also TREEBAGGER, COMPUTEOOBPREDICTION, OOBPERMUTEDVARDELTAERROR,
        %   OOBPERMUTEDVARDELTAMEANMARGIN, OOBPERMUTEDVARCOUNTRAISEMARGIN.
        ComputeOOBPredictorImportance;
        
        %NUMPREDICTORSTOSAMPLE Number of variables for random feature selection.
        %   The NumPredictorsToSample property specifies the number of predictor or
        %   feature variables to select at random for each decision split. By
        %   default, it is set to the square root of the total number of variables
        %   for classification and one third of the total number of variables for
        %   regression.
        %
        %   See also TREEBAGGER.
        NumPredictorsToSample;
        
        %MINLEAFSIZE Minimum number of observations per tree leaf.
        %   The MinLeafSize property specifies the minimum number of observations
        %   per tree leaf. By default it is 1 for classification and 5 for
        %   regression. For decision tree training, the 'MinParentSize' value is
        %   set to 2*MinLeafSize.
        %
        %   See also TREEBAGGER.
        MinLeafSize;
        
        %TREES Decision trees in the ensemble.
        %   The Trees property is a cell array of size NumTrees-by-1 containing the
        %   trees in the ensemble.
        %
        %   See also TREEBAGGER, NUMTREES.
        Trees;
        
        %NUMTREES Number of decision trees in the ensemble.
        %   The NumTrees property is a scalar equal to the number of decision trees
        %   in the ensemble.
        %
        %   See also TREEBAGGER, TREES.
        NumTrees;
        
        %CLASSNAMES Names of classes.
        %   The ClassNames property is a cell array containing the class names for
        %   the response variable Y.  This property is empty for regression trees.
        %
        %   See also TREEBAGGER.
        ClassNames;
        
        %PRIOR Prior class probabilities.
        %   The Prior property is a vector with prior probabilities for
        %   classes.  This property is empty for ensembles of regression trees.
        %
        %   See also TREEBAGGER.
        Prior = [];
        
        %COST Misclassification costs.
        %   The Cost property is a matrix with misclassification costs.
        %   This property is empty for ensembles of regression trees.
        %
        %   See also TREEBAGGER.
        Cost = [];
        
        %PREDICTORNAMES Variable names.
        %   The PredictorNames property is a cell array containing the names of the
        %   predictor variables (features). For an ensemble trained on a table,
        %   these are taken from the variable names in the table. For an ensemble
        %   trained on a matrix, these names are taken from the optional 'names'
        %   parameter, and the default names are 'x1', 'x2', etc. 
        %
        %   See also TREEBAGGER.
        PredictorNames;
        
        %METHOD Method used by trees.
        %   The Method property is 'classification' for classification ensembles
        %   and 'regression' for regression ensembles.
        %
        %   See also TREEBAGGER.
        Method;
        
        %OOBINSTANCEWEIGHT Count of out-of-bag trees for each observation.
        %   The OOBInstanceWeight property is a numeric array of size Nobs-by-1
        %   containing the number of trees used for computing OOB response for each
        %   observation.  Nobs is the number of observations in the training data
        %   used to create the ensemble.
        %
        %   See also TREEBAGGER.
        OOBInstanceWeight;
                
        %OOBPERMUTEDPREDICTORDELTAERROR Variable importance for prediction error.
        %   The OOBPermutedPredictorDeltaError property is a numeric array of size
        %   1-by-Nvars containing a measure of importance for each predictor
        %   variable (feature).  For any variable, the measure is the increase in
        %   prediction error if the values of that variable are permuted across the
        %   out-of-bag observations. This measure is computed for every tree, then
        %   averaged over the entire ensemble and divided by the standard deviation
        %   over the entire ensemble.
        %
        %   See also TREEBAGGER.
        OOBPermutedPredictorDeltaError;
        
        %OOBPERMUTEDPREDICTORDELTAMEANMARGIN Variable importance for classification
        %   margin. The OOBPermutedPredictorDeltaMeanMargin property is a numeric array
        %   of size 1-by-Nvars containing a measure of importance for each
        %   predictor variable (feature).  For any variable, the measure is the
        %   decrease in the classification margin if the values of that variable
        %   are permuted across the out-of-bag observations. This measure is
        %   computed for every tree, then averaged over the entire ensemble and
        %   divided by the standard deviation over the entire ensemble. This
        %   property is empty for regression trees.
        %
        %   See also TREEBAGGER.
        OOBPermutedPredictorDeltaMeanMargin;
        
        %OOBPERMUTEDPREDICTORCOUNTRAISEMARGIN Variable importance for raising margin.
        %   The OOBPermutedPredictorCountRaiseMargin property is a numeric array of size
        %   1-by-Nvars containing a measure of variable importance for each
        %   predictor.  For any variable, the measure is the difference between the
        %   number of raised margins and the number of lowered margins if the
        %   values of that variable are permuted across the out-of-bag
        %   observations. This measure is computed for every tree, then averaged
        %   over the entire ensemble and divided by the standard deviation over the
        %   entire ensemble. This property is empty for regression trees.
        %
        %   See also TREEBAGGER.
        OOBPermutedPredictorCountRaiseMargin;
        
        %DELTACRITERIONDECISIONSPLIT Split criterion contributions for each predictor.
        %   The DeltaCriterionDecisionSplit property is a numeric array of size
        %   1-by-Nvars of changes in the split criterion summed over splits on each
        %   variable, averaged across the entire ensemble of grown trees.
        %
        % See also TREEBAGGER, COMPACTTREEBAGGER/DELTACRITERIONDECISIONSPLIT,
        % ClassificationTree/predictorImportance, RegressionTree/predictorImportance.
        DeltaCriterionDecisionSplit;
        
        %NUMPREDICTORSPLIT Number of decision splits for each predictor.
        %   The NumPredictorSplit property is a numeric array of size 1-by-Nvars, where
        %   every element gives a number of splits on this predictor summed over
        %   all trees.
        %
        %   See also TREEBAGGER, COMPACTTREEBAGGER/NUMPREDICTORSPLIT.
        NumPredictorSplit;

        %SURROGATEASSOCIATION Predictor associations.
        %   The SurrogateAssociation  property is a matrix of size Nvars-by-Nvars with
        %   predictive measures of variable association, averaged across the entire
        %   ensemble of grown trees. If you grew the ensemble setting 'surrogate'
        %   to 'on', this matrix for each tree is filled with predictive measures
        %   of association averaged over the surrogate splits. If you grew the
        %   ensemble setting 'surrogate' to 'off' (default), SurrogateAssociation  is diagonal.
        %
        % See also COMPACTTREEBAGGER, COMPACTTREEBAGGER/SURROGATEASSOCIATION,
        % ClassificationTree/surrogateAssociation, RegressionTree/surrogateAssociation.
        SurrogateAssociation ;

        %PROXIMITY Proximity matrix for observations.
        %   The Proximity property is a numeric matrix of size Nobs-by-Nobs, where
        %   Nobs is the number of observations in the training data, containing
        %   measures of the proximity between observations.  For any two
        %   observations, their proximity is defined as the fraction of trees for
        %   which these observations land on the same leaf.  This is a symmetric
        %   matrix with 1's on the diagonal and off-diagonal elements ranging from
        %   0 to 1.
        %
        %   See also TREEBAGGER, COMPACTTREEBAGGER/PROXIMITY.
        Proximity;
        
        %OUTLIERMEASURE Measure for determining outliers.
        %   The OutlierMeasure property is a numeric array of size Nobs-by-1, where
        %   Nobs is the number of observations in the training data, containing
        %   outlier measures for each observation.
        %
        %   See also TREEBAGGER, COMPACTTREEBAGGER/OUTLIERMEASURE.
        OutlierMeasure;
    end

        
    properties(SetAccess=public,GetAccess=public,Dependent=true)
        %DEFAULTYFIT Default value returned by PREDICT and OOBPREDICT.
        %   The DefaultYfit property controls what predicted value is returned when
        %   no prediction is possible, for example when the OOBPREDICT method needs
        %   to predict for an observation that is in-bag for all trees in the
        %   ensemble.  For classification, you can set this property to either ''
        %   or 'MostPopular'.  If you choose 'MostPopular' (default), the property
        %   value becomes the name of the most probable class in the training data.
        %   For regression, you can set this property to any numeric scalar.  The
        %   default is the mean of the response for the training data.  If you set
        %   this property to '' for classification or NaN for regression, the
        %   in-bag observations are excluded from computation of the out-of-bag
        %   error and margin.
        %
        %   See also TREEBAGGER, OOBPREDICT, PREDICT, OOBINDICES,
        %   COMPACTTREEBAGGER/DEFAULTYFIT.
        DefaultYfit;
    end
    
    properties(SetAccess=protected,GetAccess=public,Hidden=true)        
        FBoot = 1;
        TreeArgs = {};
        ComputeOOBVarImp = false;
        NVarToSample = [];
        MinLeaf = [];
    end
    
    properties(SetAccess=protected,GetAccess=public,Dependent=true,Hidden=true)
        NTrees;
        VarNames;
        OOBPermutedVarDeltaError;
        OOBPermutedVarDeltaMeanMargin;
        OOBPermutedVarCountRaiseMargin;
        DeltaCritDecisionSplit;
        NVarSplit
        VarAssoc;
    end

    % Properties related to table support, and their access methods
    properties(Dependent=true,GetAccess=public,SetAccess=public,Hidden=true)
        VariableRange = {};
        TableInput = false;
    end

    methods
        function vr = get.VariableRange(this)
            try
                vr = this.DataSummary.VariableRange;
            catch
                vr = {}; % field may be missing for old saved object
            end
        end
        function this = set.VariableRange(this,vr)
            this.DataSummary.VariableRange = vr;
        end
        function ti = get.TableInput(this)
            try
                ti = this.DataSummary.TableInput;
            catch
                ti = false; % field may be missing for old saved object
            end
        end
        function this = set.TableInput(this,t)
            this.DataSummary.TableInput = t;
        end
    end


    methods
        function bagger = TreeBagger(NumTrees,X,Y,varargin)
            %TREEBAGGER Create an ensemble of bagged decision trees.
            %   B = TREEBAGGER(NUMTREES,TBL,Y) creates an ensemble B of NUMTREES decision
            %   trees for predicting response Y from data in the table TBL. Y can be
            %   any of the following:
            % 
            %      1. An array of true class labels for classification, or numeric
            %         values for regression. Class labels can be a categorical array,
            %         logical vector, numeric vector, or cell array of strings.
            %      2. The name of a variable in TBL. This variable is used as the
            %         response Y, and the remaining variables in TBL are used as
            %         predictors.
            %      3. A formula string such as 'y ~ x1 + x2 + x3' specifying that the
            %         variable y is to be used as the response, and the other variables
            %         in the formula are predictors. Any table variables not listed in
            %         the formula are not used.
            %
            %   TREEBAGGER converts Y to a cell array of strings for classification.
            %
            %   B = TREEBAGGER(NUMTREES,X,Y) is an alternative syntax that accepts X as
            %   an N-by-P numeric matrix of predictors with one row per observation and
            %   one column per predictor. Y is the response and is an array of N class
            %   labels for classification, or numeric values for regression.
            %
            %   Use of a matrix X rather than a table TBL saves both memory and
            %   execution time.
            %
            %   By default, TREEBAGGER builds an ensemble of classification trees. The
            %   function can build an ensemble of regression trees by setting the
            %   optional input argument 'Method' to 'regression'.
            %
            %   B = TREEBAGGER(NTREES,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies
            %   optional parameter name/value pairs:
            %
            %     'Cost'        Square matrix, where COST(I,J) is the cost of
            %                   classifying a point into class J if its true class is
            %                   I. Alternatively, COST can be a structure S having two
            %                   fields: S.ClassificationCosts containing the cost
            %                   matrix C, and S.ClassNames containing the class names
            %                   and defining the ordering of classes used for the rows
            %                   and columns of the cost matrix. Default: COST(I,J)=1 if
            %                   I~=J, and COST(I,J)=0 if I=J.
            %                NOTE: TreeBagger generates bootstrap replicas (in-bag
            %                      samples) by oversampling classes with large
            %                      misclassification costs and undersampling classes
            %                      with small misclassification costs. Out-of-bag (OOB)
            %                      samples consequently have fewer observations from
            %                      classes with large misclassification costs and more
            %                      observations from classes with small
            %                      misclassification costs. For small datasets and
            %                      highly skewed costs, the number of OOB observations
            %                      per class can be very low. The OOB error then has
            %                      large variance and may be hard to interpret.
            %     'InBagFraction'       Fraction of input data to sample with replacement
            %                   from the input data for growing each new tree.
            %                   By default set to 1.
            %     'SampleWithReplacement'   'on' to sample with replacement or 'off' to
            %                   sample without replacement. If you sample without 
            %                   replacement, you need to set 'InBagFraction' to a value 
            %                   less than one.  Default is 'on'.
            %     'OOBPrediction'     'on' to store info on what observations are out of bag
            %                   for each tree. This info can be used by OOBPREDICT to
            %                   compute the out-of-bag predicted class probabilities
            %                   for each tree in the ensemble.  Default is 'off'.
            %     'OOBPredictorImportance'   'on' to store out-of-bag estimates of feature
            %                   importance in the ensemble.  Default is 'off'.
            %                   Specifying 'on' also sets the 'OOBPrediction' value to 'on'.
            %     'Method'      'classification' or 'regression'. Regression requires a
            %                   numeric Y.
            %     'NumPredictorsToSample' Number of variables to select at random for each
            %                   decision split. Default is the square root of the
            %                   number of variables for classification and one third of
            %                   the number of variables for regression. Valid values
            %                   are 'all' or a positive integer. Setting this argument
            %                   to any valid value but 'all' invokes Breiman's 'random
            %                   forest' algorithm.
            %     'NumPrint'    Number of training cycles (grown trees) after which
            %                   TREEBAGGER displays a diagnostic message showing training
            %                   progress. Default is no diagnostic messages.
            %     'MinLeafSize'     Minimum number of observations per tree leaf. Default
            %                   is 1 for classification and 5 for regression.
            %     'Options'     A struct that contains options specifying whether to use
            %                   parallel computation when growing the ensemble of decision
            %                   trees, and options specifying how to use random numbers
            %                   when drawing replicates of the original data. This argument
            %                   can be created by a call to STATSET. TREEBAGGER uses the
            %                   following fields:
            %                    'UseParallel'
            %                    'UseSubstreams'
            %                    'Streams'
            %                   For information on these fields see PARALLELSTATS.
            %                NOTE: If 'UseParallel' is TRUE and 'UseSubstreams' is
            %                      FALSE, then the length of 'Streams' must equal the
            %                      number of workers used by TREEBAGGER. If a parallel
            %                      pool is already open, this will be the size of the
            %                      parallel pool. If a parallel pool is not already
            %                      open, then MATLAB may try to open a pool for you
            %                      (depending on your installation and preferences). To
            %                      ensure more predictable results, it is best to use
            %                      the PARPOOL command and explicitly create a parallel
            %                      pool prior to invoking TREEBAGGER with 'UseParallel'
            %                      set to TRUE.
            %     'Prior'       Prior probabilities for each class.
            %                     Specify as one of:
            %                       * A string:
            %                         - 'empirical' determines class probabilities from
            %                            class frequencies in Y
            %                         - 'uniform' sets all class probabilities equal
            %                       * A vector (one scalar value for each class)
            %                       * A structure S with two fields: S.ClassProbs
            %                         containing a vector of class probabilities, and
            %                         S.ClassNames containing the class names and
            %                         defining the ordering of classes used for the
            %                         elements of this vector. If you pass numeric
            %                         values, TreeBagger normalizes them to add up to
            %                         one. Default: 'empirical'
            %                NOTE: TreeBagger generates bootstrap replicas (in-bag
            %                      samples) by oversampling classes with large prior
            %                      probabilities and undersampling classes with small
            %                      prior probabilities. Out-of-bag (OOB) samples
            %                      consequently have fewer observations from classes
            %                      with large prior probabilities and more observations
            %                      from classes with small prior probabilities. For
            %                      small datasets and highly skewed prior
            %                      probabilities, the number of OOB observations per
            %                      class can be very low. The OOB error then has large
            %                      variance and may be hard to interpret.
            %
            %   In addition to the optional arguments above, this method accepts some
            %   optional arguments for FITCTREE or FITRTREE. For more details, see the
            %   <a href="matlab:doc TreeBagger/TreeBagger">TreeBagger reference page</a>.
            %
            %   See also TreeBagger, oobPredict, fitctree, fitrtree, groupingvariable,
            %   RandStream, statset, statget, parfor, parallelstats.
            if nargin > 2
                Y = convertStringsToChars(Y);
            end
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            % Process inputs for the bagger
            growTreeArgs = {{'nprint','numprint'}             'options'};
            growTreeDefs = {                    0 statset('TreeBagger')};
            [nprint,parallelOptions,~,makeArgs] ...
                = internal.stats.parseArgs(growTreeArgs,growTreeDefs,varargin{:});

            % Make an empty bagger
            bagger = init(bagger,X,Y,makeArgs{:});
            
            % Add trees
            bagger = growTrees(bagger,NumTrees,'Options', parallelOptions, 'nprint', nprint);
            
            % Version 3 introduced in 16b. Adds the InBagIndices property
            % and quantile forest.
            bagger.Version = 3;
        end
        
        function bagger = growTrees(bagger,NumTrees,varargin)
            %GROWTREES Train additional trees and add to ensemble.
            %   B = GROWTREES(B,NUMTREES) grows NUMTREES new trees and appends them to
            %   those already stored in the ensemble B.
            %
            %   B = GROWTREES(B,NUMTREES,'PARAM1',val1,'PARAM2',val2,...) specifies
            %   optional parameter name/value pairs:
            %
            %     'NumPrint'      Specifies that a diagnostic message showing training
            %                   progress should display after every VALUE training
            %                   cycles (grown trees). Default is no diagnostic messages.
            %
            %     'Options'     A struct that contains options specifying whether to use
            %                   parallel computation when growing the ensemble of decision
            %                   trees, and options specifying how to use random numbers 
            %                   when drawing replicates of the original data. This argument 
            %                   can be created by a call to STATSET. TREEBAGGER uses the 
            %                   following options:
            %                    'UseParallel'
            %                    'UseSubstreams'
            %                    'Streams'
            %                   For information on these options see PARALLELSTATS.
            %      NOTE: If 'UseParallel' is TRUE and 'UseSubstreams' is FALSE, then the
            %                   length of 'Streams' must equal the number of workers used by
            %                   TREEBAGGER.  If a parallel pool is already open, this will 
            %                   be the size of the parallel pool.  If a parallel pool is 
            %                   not already open, then MATLAB may try to open a pool 
            %                   for you (depending on your installation and preferences).
            %                   To ensure more predictable results, it is best to use the PARPOOL
            %                   command and explicitly create a parallel pool prior to invoking 
            %                   TREEBAGGER with 'UseParallel' set to TRUE.
            %
            %   See also TREEBAGGER, ClassificationTree.fit, RegressionTree.fit,
            %   RANDSTREAM, STATSET, STATGET, PARFOR, PARPOOL, PARALLELSTATS.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            % Process inputs for the bagger
            baggerArgs = {{'nprint','numprint'}         'options'};
            baggerDefs = {       0                statset('TreeBagger')};
            [nprint,parallelOptions] ...
                = internal.stats.parseArgs(baggerArgs,baggerDefs,varargin{:});
            
            % Check print-out frequency
            if isempty(nprint) || ~isnumeric(nprint) || numel(nprint)~=1 || nprint<0
                warning(message('stats:TreeBagger:growTrees:BadPrintoutFreq'));
            end
            
            % Check number of trees
            if ~isnumeric(NumTrees) || NumTrees<0
                error(message('stats:TreeBagger:growTrees:BadTreeNum'));
            end
            if NumTrees==0
                return;
            end
            
            % Process options for parallel execution and random number control
            [useParallel, RNGscheme, poolsz] = ...
                internal.stats.parallel.processParallelAndStreamOptions(parallelOptions,true);
            usePool = useParallel && poolsz>0;
            
            % Determine data size and number of vars
            [N,Nvars] = size(bagger.PrivateX);
            
            % How many trees have been grown, so far?
            NTreesBefore = bagger.NTrees;
            
            % Resume growing old classregtrees?
            doclassregtree = false;
            if NTreesBefore>0 && isa(bagger.Trees{1},'classregtree')
                doclassregtree = true;
            end
            
            % Prepare oob arrays
            if bagger.ComputeOOBPrediction
                % Expand array of oob indices
                if isempty(bagger.OOBIndices)
                    bagger.OOBIndices = false(N,NumTrees);
                else
                    bagger.OOBIndices(1:N,end+1:end+NumTrees) = false;
                end
                
                % Expand array with delta error and delta margin
                if bagger.ComputeOOBVarImp
                    bagger.PrivOOBPermutedVarDeltaError(end+1:end+NumTrees,1:Nvars) ...
                        = zeros(NumTrees,Nvars);
                    if bagger.Method(1)=='c'
                        bagger.PrivOOBPermutedVarDeltaMeanMargin(end+1:end+NumTrees,1:Nvars) ...
                            = zeros(NumTrees,Nvars);
                        bagger.PrivOOBPermutedVarCountRaiseMargin(end+1:end+NumTrees,1:Nvars) ...
                            = zeros(NumTrees,Nvars);
                    end
                end
            end
            
            if bagger.Method(1)=='r'
                fsample = bagger.InBagFraction;
                Nsample = ceil(N*fsample);
                if isempty(bagger.InBagIndices)
                    bagger.InBagIndices = zeros(Nsample,NumTrees);
                else
                    bagger.InBagIndices(1:Nsample,end+1:end+NumTrees) = zeros(Nsample,NumTrees);
                end
            end
                        
            % Dump all bagger properties needed by loopBody
            fboot = bagger.FBoot;
            sampleWithReplacement = bagger.SampleWithReplacement;
            x = bagger.PrivateX;
            y = bagger.Y;
            w = bagger.W;
            method = bagger.Method;
            ds = bagger.DataSummary;
            cs = bagger.ClassSummary;
            mp = bagger.ModelParams;
            computeOOBPrediction = bagger.ComputeOOBPrediction;
            computeOOBVarImp = bagger.ComputeOOBVarImp;
            if computeOOBVarImp
                xorig = bagger.X;
            else
                xorig = [];
            end
            args = {doclassregtree fboot sampleWithReplacement x y w method ...
                ds cs mp computeOOBPrediction computeOOBVarImp nprint usePool xorig};

            % Set up running counts for optional progress reports
            if nprint>0
                if usePool
                    % We are asked for a periodic tally of progress in growing
                    % the ensemble. On each worker, we maintain a running count
                    % of the number of trees created during this function
                    % invocation. This count is on an individual worker basis,
                    % not a cumulative count across all the workers.
                    % Here we initialize the count on each worker.
                    %
                    % Since usePool is true, there is a parallel pool open,
                    % and the unguarded call to get pool size should
                    % never error.
                    parfor i=1:internal.stats.parallel.getParallelPoolSize
                        internal.stats.parallel.statParallelStore('mylabindex', i);
                        internal.stats.parallel.statParallelStore('ntreesGrown',0);
                    end
                else
                    % Initialize running count variables on the client
                    internal.stats.parallel.statParallelStore('mylabindex', 1);
                    internal.stats.parallel.statParallelStore('ntreesGrown',0);
                end
            end
            
            % Grow trees
            [trees, ...
                slicedOOBIndices, ...
                slicedInBagIndices, ...
                slicedPrivOOBPermutedVarDeltaError, ...
                slicedPrivOOBPermutedVarDeltaMeanMargin, ...
                slicedPrivOOBPermutedVarCountRaiseMargin] = ...
                localGrowTrees(NumTrees, useParallel, RNGscheme, args);
                            
            % Append scores for the newly grown trees to the accumulated list
            if computeOOBPrediction
                bagger.OOBIndices(:,NTreesBefore+1:NTreesBefore+NumTrees) = slicedOOBIndices;
                
                if computeOOBVarImp
                    bagger.PrivOOBPermutedVarDeltaError(NTreesBefore+1:NTreesBefore+NumTrees,:) = ...
                        slicedPrivOOBPermutedVarDeltaError;
                    if bagger.Method(1)=='c'
                        bagger.PrivOOBPermutedVarDeltaMeanMargin(NTreesBefore+1:NTreesBefore+NumTrees,:) = ...
                            slicedPrivOOBPermutedVarDeltaMeanMargin;
                        bagger.PrivOOBPermutedVarCountRaiseMargin(NTreesBefore+1:NTreesBefore+NumTrees,:) = ...
                            slicedPrivOOBPermutedVarCountRaiseMargin;
                    end
                end
            end
            
            if bagger.Method(1)=='r'
                bagger.InBagIndices(:,NTreesBefore+1:NTreesBefore+NumTrees) = slicedInBagIndices;
            end
            
            % Store grown trees in the Compact object
            bagger.Compact = addTrees(bagger.Compact,trees);                        
        end %-growTrees

        function cmp = compact(bagger)
            %COMPACT Compact ensemble of decision trees.
            %   Return an object of class CompactTreeBagger holding the structure of
            %   the trained ensemble. The class is more compact than the full
            %   TreeBagger class because it does not contain information for growing
            %   more trees for the ensemble. In particular it does not contain X and Y
            %   used for training. 
            %
            % See also COMPACTTREEBAGGER, Method.

            cmp = bagger.Compact;
        end
        
        function bagger = fillprox(bagger,varargin)
            %FILLPROX Proximity matrix for training data.
            %   B = FILLPROX(B) computes a proximity matrix for the training data
            %   and stores it in the Properties field of B.
            %
            %   B = FILLPROX(B,'PARAM1',val1,'PARAM2',val2,...) specifies
            %   optional parameter name/value pairs:
            %
            %     'Trees'       Either 'all' or a vector of indices of the trees in the
            %                   ensemble to be used in computing the proximity matrix.
            %                   Default is 'all'.
            %     'NumPrint'      Number of training cycles (grown trees) after which
            %                   TREEBAGGER displays a diagnostic message showing training
            %                   progress. Default is no diagnostic messages.
            %
            % See also TREEBAGGER, COMPACTTREEBAGGER/PROXIMITY,
            % COMPACTTREEBAGGER/OUTLIERMEASURE.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            bagger = fillProximities(bagger,varargin{:});
        end
        
        function bagger = append(bagger,other)
            %APPEND Append new trees to ensemble.
            %   B = APPEND(B,OTHER) appends the trees from the OTHER
            %   ensemble to those in B.  This method checks for consistency
            %   of the X and Y properties of the two ensembles, as well as
            %   consistency of their compact objects and out-of-bag
            %   indices, before appending the trees.  The output ensemble B
            %   takes the training parameters such as InBagFraction, Prior,
            %   Cost, and other from the B input. There is no attempt to
            %   check if these training parameters are consistent between
            %   the two objects.
            %
            % See also TREEBAGGER, COMPACTTREEBAGGER/COMBINE.
            
            % Check X
            if ~isequaln(bagger.PrivateX,other.PrivateX)
                error(message('stats:TreeBagger:append:IncompatibleX'));
            end
            
            % Check Y
            if ~isequaln(bagger.Y,other.Y)
                error(message('stats:TreeBagger:append:IncompatibleY'));
            end
            
            % Check ensemble type
            if ~strcmpi(bagger.Method,other.Method)
                error(message('stats:TreeBagger:append:IncompatibleMethod'));
            end
            
            % Check OOB info
            if bagger.ComputeOOBPrediction~=other.ComputeOOBPrediction ...
                    || bagger.ComputeOOBVarImp~=other.ComputeOOBVarImp
                error(message('stats:TreeBagger:append:IncompatibleOOB'));
            end
            
            % Combine Compact objects
            bagger.Compact = combine(bagger.Compact,other.Compact);
            
            % Combine OOB indices
            nTrees = other.NTrees;
            if nTrees>0 && other.ComputeOOBPrediction
                if size(other.OOBIndices,1)~=size(bagger.OOBIndices,1) || ...
                        size(other.OOBIndices,2)~=nTrees
                    error(message('stats:TreeBagger:append:BadOOBIndices'));
                end
                bagger.OOBIndices(:,end+1:end+nTrees) = other.OOBIndices;
            end
            
            % Combine in-bag indices
            nTrees = other.NTrees;
            if nTrees>0 && bagger.Method(1)=='r'
                if size(other.InBagIndices,1)~=size(bagger.InBagIndices,1) || ...
                        size(other.InBagIndices,2)~=nTrees
                    error(message('stats:TreeBagger:append:BadInBagIndices'));
                end
                bagger.InBagIndices(:,end+1:end+nTrees) = other.InBagIndices;
            end
            
            % Combine OOB permuted info
            if nTrees>0 && other.ComputeOOBVarImp
                % Error
                if size(bagger.PrivOOBPermutedVarDeltaError,2) ...
                        ~=size(other.PrivOOBPermutedVarDeltaError,2) || ...
                        size(other.PrivOOBPermutedVarDeltaError,1)~=nTrees
                    error(message('stats:TreeBagger:append:BadOOBError'));
                end
                bagger.PrivOOBPermutedVarDeltaError(end+1:end+nTrees,:) = ...
                    other.PrivOOBPermutedVarDeltaError;

                % Mean margin
                if size(bagger.PrivOOBPermutedVarDeltaMeanMargin,2) ...
                        ~=size(other.PrivOOBPermutedVarDeltaMeanMargin,2) || ...
                        size(other.PrivOOBPermutedVarDeltaMeanMargin,1)~=nTrees
                    error(message('stats:TreeBagger:append:BadOOBMargin'));
                end
                bagger.PrivOOBPermutedVarDeltaMeanMargin(end+1:end+nTrees,:) = ...
                    other.PrivOOBPermutedVarDeltaMeanMargin;

                % Raised margins
                if size(bagger.PrivOOBPermutedVarCountRaiseMargin,2) ...
                        ~=size(other.PrivOOBPermutedVarCountRaiseMargin,2) || ...
                        size(other.PrivOOBPermutedVarCountRaiseMargin,1)~=nTrees
                    error(message('stats:TreeBagger:append:BadOOBCountMargin'));
                end
                bagger.PrivOOBPermutedVarCountRaiseMargin(end+1:end+nTrees,:) = ...
                    other.PrivOOBPermutedVarCountRaiseMargin;
            end
            
            % Combine proximities
            if isempty(other.PrivProx)
                bagger.PrivProx = [];
            else
                if ~isempty(bagger.PrivProx)
                    if numel(bagger.PrivProx)~=numel(other.PrivProx)
                        error(message('stats:TreeBagger:append:BadProx'));
                    end
                    bagger.PrivProx = ...
                        (bagger.NTrees*bagger.PrivProx + nTrees*other.PrivProx) / ...
                        (bagger.NTrees + nTrees);
                end
            end
        end
    end
    
    methods(Static)
        function this = loadobj(obj)
            if ~isstruct(obj) && obj.Version>=2
                % 16a or later
                this = obj;
            elseif ~isempty(obj.Compact.Trees) ...
                    && isa(obj.Compact.Trees{1},'classreg.learning.Predictor')
                % 13a or later
                this = obj;
                this.PrivateX = this.X;
            else % 12b or earlier  
                sampleWithReplacement = 'off';
                if obj.SampleWithReplacement
                    sampleWithReplacement = 'on';
                end
                computeOOBPrediction = 'off';
                if obj.ComputeOOBPrediction
                    computeOOBPrediction = 'on';
                end
                computeOOBVarImp = 'off';
                if obj.ComputeOOBVarImp
                    computeOOBVarImp = 'on';
                end
                prune = 'off';
                if obj.Prune
                    prune = 'on';
                end
                mergeLeaves = 'off';
                if obj.MergeLeaves
                    mergeLeaves = 'on';
                end
                
                weightsSet = false;
                w = ones(size(obj.X,1),1);
                if isfield(obj,'W')
                    w = obj.W;
                    weightsSet = true;
                end
                
                this = TreeBagger(0,obj.X,obj.Y,...
                    'fboot',obj.FBoot,'samplewithreplacement',sampleWithReplacement,...
                    'oobpred',computeOOBPrediction,'oobvarimp',computeOOBVarImp,...
                    'prune',prune,'mergeleaves',mergeLeaves,...
                    'nvartosample',obj.NVarToSample,'minleaf',obj.MinLeaf,...
                    'weights',w,obj.TreeArgs{:});
                this.Compact = obj.Compact;
                this.OOBIndices = obj.OOBIndices;
                this.PrivOOBPermutedVarDeltaError = obj.PrivOOBPermutedVarDeltaError;
                this.PrivOOBPermutedVarDeltaMeanMargin = obj.PrivOOBPermutedVarDeltaMeanMargin;
                this.PrivOOBPermutedVarCountRaiseMargin = obj.PrivOOBPermutedVarCountRaiseMargin;
                
                if ~isempty(obj.PriorStruct)
                    newgrp = cellstr(this.ClassSummary.NonzeroProbClasses);
                    oldgrp = obj.PriorStruct.group;
                    [~,loc] = ismember(newgrp,oldgrp);
                    this.ClassSummary.Prior = obj.PriorStruct.prob(loc);
                end
                
                if ~isempty(obj.CostStruct)
                    newgrp = cellstr(this.ClassSummary.NonzeroProbClasses);
                    oldgrp = obj.CostStruct.group;
                    [~,loc] = ismember(newgrp,oldgrp);
                    this.ClassSummary.Cost = obj.CostStruct.cost(loc,loc);
                end

                % Need to adjust weights for prior probabilities and
                % misclassification costs only if the weights were not
                % filled in the older object and if prior or cost was set
                % to a non-default value.
                if ~weightsSet && (~isempty(obj.PriorStruct) || ~isempty(obj.CostStruct))
                    [~,this.W] = classreg.learning.internal.adjustPrior(...
                        this.ClassSummary,classreg.learning.internal.ClassLabel(this.Y),w);
                end
            end
        end
    end

    
    methods(Hidden=true,Static=true)
        function a = empty(varargin),      throwUndefinedError(); end %#ok<STOUT>
    end
    
    
    methods(Hidden=true)
        function disp(obj)
            internal.stats.displayClassName(obj);

            % Display body
            fprintf(1,'%s\n',getString(message('stats:TreeBagger:DispHeader',obj.NTrees)));
            sx = ['[' num2str(size(obj.PrivateX,1)) 'x' num2str(size(obj.PrivateX,2)) ']'];
            sy = ['[' num2str(size(obj.Y,1)) 'x' num2str(size(obj.Y,2)) ']'];
            fprintf(1,'%30s: %20s\n','Training X',sx);
            fprintf(1,'%30s: %20s\n','Training Y',sy);
            fprintf(1,'%30s: %20s\n','Method',obj.Method);
            fprintf(1,'%30s: %20i\n','NumPredictors',length(obj.VarNames));
            fprintf(1,'%30s: %20s\n','NumPredictorsToSample',num2str(obj.NVarToSample));
            fprintf(1,'%30s: %20i\n','MinLeafSize',obj.MinLeaf);
            fprintf(1,'%30s: %20g\n','InBagFraction',obj.FBoot);
            fprintf(1,'%30s: %20i\n','SampleWithReplacement',obj.SampleWithReplacement);
            fprintf(1,'%30s: %20i\n','ComputeOOBPrediction',obj.ComputeOOBPrediction);
            fprintf(1,'%30s: %20i\n','ComputeOOBPredictorImportance',obj.ComputeOOBVarImp);
            if ~isempty(obj.PrivProx)
                sprox = ['[' num2str(size(obj.PrivateX,1)) 'x' num2str(size(obj.PrivateX,1)) ']'];
            else
                sprox = '[]';
            end
            fprintf(1,'%30s: %20s\n','Proximity',sprox);
            if obj.Method(1)=='c'
                sform = ' %s';
                if ~isempty(obj.Prior) || ~isempty(obj.Cost)
                    sform = ' %15s';
                end
                fprintf(1,'%30s:','ClassNames');
                for i=1:length(obj.ClassNames)
                    fprintf(1,sform,['''' obj.ClassNames{i} '''']);
                end
                fprintf(1,'\n');
            end
            
            internal.stats.displayMethodsProperties(obj);
        end
    end
    
    
    methods
        function x = get.X(bagger)
            x = bagger.PrivateX;
            % Decode categorical variables into original values
            if bagger.TableInput
                t = array2table(x,'VariableNames',bagger.PredictorNames);
                for j=1:size(x,2)
                    vrj = bagger.VariableRange{j};
                    newx = decodeX(x(:,j),vrj);
                    t.(bagger.PredictorNames{j}) = newx;
                end
                x = t;
            end
        end
        function bagger = set.X(bagger,x)
            bagger.PrivateX = x;
        end
        
        function fboot = get.InBagFraction(bagger)
            fboot = bagger.FBoot;
        end
        
        function treeargs = get.TreeArguments(bagger)
            treeargs = bagger.TreeArgs;
        end
        
        function computeOOBVarImp = get.ComputeOOBPredictorImportance(bagger)
            computeOOBVarImp = bagger.ComputeOOBVarImp;
        end
        
        function nvartosample = get.NumPredictorsToSample(bagger)
            nvartosample = bagger.NVarToSample;
        end
        
        function minleaf = get.MinLeafSize(bagger)
            minleaf = bagger.MinLeaf;
        end
        
        function trees = get.Trees(bagger)
            trees = bagger.Compact.Trees;
        end
        
        function n = get.NTrees(bagger)
            n = length(bagger.Trees);
        end
        
        function n = get.NumTrees(bagger)
            n = bagger.NTrees;
        end
        
        function cnames = get.ClassNames(bagger)
            cnames = bagger.Compact.ClassNames;
        end

        function prior = get.Prior(this)
            K = length(this.ClassSummary.ClassNames);
            prior = zeros(1,K);
            [~,pos] = ismember(this.ClassSummary.NonzeroProbClasses,...
                this.ClassSummary.ClassNames);
            prior(pos) = this.ClassSummary.Prior;            
        end
        
        function cost = get.Cost(this)
            K = length(this.ClassSummary.ClassNames);
            if isempty(this.ClassSummary.Cost)
                cost = ones(K) - eye(K);
            else
                cost = zeros(K);
                [~,pos] = ismember(this.ClassSummary.NonzeroProbClasses,...
                    this.ClassSummary.ClassNames);
                cost(pos,pos) = this.ClassSummary.Cost;
                unmatched = 1:K;
                unmatched(pos) = [];
                cost(:,unmatched) = NaN;
                cost(1:K+1:end) = 0;
            end
        end
        
        function vnames = get.VarNames(bagger)
            vnames = bagger.Compact.VarNames;
        end
        
        function vnames = get.PredictorNames(bagger)
            vnames = bagger.VarNames;
        end
        
        function meth = get.Method(bagger)
            meth = bagger.Compact.Method;
        end
        
        function weights = get.OOBInstanceWeight(bagger)
            % Check if OOB info was filled
            if ~bagger.ComputeOOBPrediction
                error(message('stats:TreeBagger:OOBInstanceWeight:InvalidProperty'));
            end
            
            % Get weights
            weights = sum(bagger.OOBIndices,2);
        end
               
        function deltacrit = get.DeltaCritDecisionSplit(bagger)
            deltacrit = bagger.Compact.DeltaCritDecisionSplit;
        end
        
        function deltacrit = get.DeltaCriterionDecisionSplit(bagger)
            deltacrit = bagger.DeltaCritDecisionSplit;
        end
         
        function nsplit = get.NVarSplit(bagger)
            nsplit = bagger.Compact.NVarSplit;
        end
        
        function nsplit = get.NumPredictorSplit(bagger)
            nsplit = bagger.NVarSplit;
        end
        
        function assoc = get.VarAssoc(bagger)
            assoc = bagger.Compact.VarAssoc;
        end
        
        function assoc = get.SurrogateAssociation(bagger)
            assoc = bagger.VarAssoc;
        end

        function delta = get.OOBPermutedVarDeltaError(bagger)
            % Check if permutation info was filled
            if ~bagger.ComputeOOBVarImp
                error(message('stats:TreeBagger:OOBPermutedVarDeltaError:InvalidProperty'));
            end
            
            % Get error shifts due to variable permutation
            delta = normalizeMean1(bagger.PrivOOBPermutedVarDeltaError);
        end
        
        function delta = get.OOBPermutedPredictorDeltaError(bagger)
            delta = bagger.OOBPermutedVarDeltaError;
        end
        
        function delta = get.OOBPermutedVarDeltaMeanMargin(bagger)
            % Check if permutation info was filled
            if ~bagger.ComputeOOBVarImp
                error(message('stats:TreeBagger:OOBPermutedVarDeltaMeanMargin:InvalidProperty'));
            end
            
            % Get shifts in mean margin due to variable permutations
            delta = normalizeMean1(bagger.PrivOOBPermutedVarDeltaMeanMargin);
        end
        
        function delta = get.OOBPermutedPredictorDeltaMeanMargin(bagger)
            delta = bagger.OOBPermutedVarDeltaMeanMargin;
        end
        
        function delta = get.OOBPermutedVarCountRaiseMargin(bagger)
            % Check if permutation info was filled
            if ~bagger.ComputeOOBVarImp
                error(message('stats:TreeBagger:OOBPermutedVarCountRaiseMargin:InvalidProperty'));
            end
            
            % Get shifts in raised-lowered margin counts due to variable
            % permutations
            delta = normalizeMean1(bagger.PrivOOBPermutedVarCountRaiseMargin);
        end
        
        function delta = get.OOBPermutedPredictorCountRaiseMargin(bagger)
            delta = bagger.OOBPermutedVarCountRaiseMargin;
        end
        
        function prox = get.Proximity(bagger)
            % Check if proximities were computed
            if isempty(bagger.PrivProx)
                error(message('stats:TreeBagger:Proximity:InvalidProperty'));
            end
            
            % Get proximities
            prox = squareform(bagger.PrivProx);
            N = size(bagger.PrivateX,1);
            prox(1:N+1:end) = 1;
        end
        
        function outlier = get.OutlierMeasure(bagger)
            % Check if proximities have been computed
            if isempty(bagger.PrivProx)
                error(message('stats:TreeBagger:OutlierMeasure:InvalidProperty'));
            end
            
            % Compute outliers
            if bagger.Method(1)=='c'
                outlier = outlierMeasure(bagger.Compact,bagger.Proximity,...
                    'data','proximity','labels',bagger.Y);
            else
                outlier = outlierMeasure(bagger.Compact,bagger.Proximity,...
                    'data','proximity');
            end
        end
        
        function yfit = get.DefaultYfit(bagger)
            yfit = bagger.Compact.DefaultYfit;
        end
        
        function bagger = set.DefaultYfit(bagger,yfit)
            bagger.Compact = setDefaultYfit(bagger.Compact,yfit);
        end
        
        function sc = get.DefaultScore(bagger)
            sc = bagger.Compact.DefaultScore;
        end
        
        function bagger = set.DefaultScore(bagger,score)
            bagger.Compact.DefaultScore = score;
        end
    end
   
    methods(Access=protected)
        function bagger = init(bagger,x,y,varargin)
            % Process inputs for the bagger
            baggerArgs = {{'fboot','inbagfraction'} 'samplewithreplacement' ...
                          {'oobpred','oobprediction'} {'oobvarimp','oobpredictorimportance'} ...
                          'method' 'names' 'splitmin'};
            baggerDefs = {      1                    'on' ...
                              'off'       'off' ...
                  'classification'     {}          []};
            [bagger.FBoot,...
                samplemethod,oobpred,oobvarimp, ...
                method,varnames,splitmin,~,bagger.TreeArgs] ...
                = internal.stats.parseArgs(baggerArgs,baggerDefs,varargin{:});
            
            % Check status and inputs
            checkOnOff = ...
                @(x) ischar(x) && (strcmpi(x,'off') || strcmpi(x,'on'));

            if ~isnumeric(bagger.FBoot) ...
                    || bagger.FBoot<=0 || bagger.FBoot>1
                error(message('stats:TreeBagger:init:BadFBoot'));
            end
                        
            if ~checkOnOff(samplemethod)
                error(message('stats:TreeBagger:init:BadSampleWithReplacement'));
            end
            bagger.SampleWithReplacement = strcmpi(samplemethod,'on');
            
            if ~checkOnOff(oobpred)
                error(message('stats:TreeBagger:init:BadOOBPred'));
            end
            bagger.ComputeOOBPrediction = strcmpi(oobpred,'on');

            if ~checkOnOff(oobvarimp)
                error(message('stats:TreeBagger:init:BadOOBVarImp'));
            end
            bagger.ComputeOOBVarImp = strcmpi(oobvarimp,'on');

            allowedVals = {'classification' 'regression'};
            tf = strncmpi(method,allowedVals,length(method));
            if isempty(method) || ~ischar(method) || ~any(tf)
                error(message('stats:TreeBagger:init:BadMethod'));
            end
            method = allowedVals{tf};
            
            [bagger.ModelParams,extraArgs] = classreg.learning.modelparams.TreeParams.make(...
                method,bagger.TreeArgs{:});

            if islogical(y)
                y = double(y);
            end
            
            if strcmp(method(1),'c')
                % Prepare data. Require that y be treated as a string
                % label, though it is returned as a ClassLabel.
                [bagger.PrivateX,y,bagger.W,bagger.DataSummary,classSummary] = ...
                    ClassificationTree.prepareData(x,y,...
                    'predictornames',varnames,'ResponseString',true,extraArgs{:});
                bagger.Y = cellstr(y);
                [classSummary,bagger.W] = ...
                    classreg.learning.internal.adjustPrior(classSummary,y,bagger.W);
                bagger.ClassSummary = classSummary;
            else
                [bagger.PrivateX,bagger.Y,bagger.W,bagger.DataSummary] = ...
                    RegressionTree.prepareData(...
                    x,y,'predictornames',varnames,extraArgs{:});
            end
            
            if isempty(bagger.ModelParams.Prune)
                bagger.ModelParams.Prune = 'off';
            end
            bagger.Prune = strcmpi(bagger.ModelParams.Prune,'on');
            if bagger.Prune
                warning(message('stats:TreeBagger:init:BadPruneValue'));
            end
            
            if isempty(bagger.ModelParams.MergeLeaves)
                bagger.ModelParams.MergeLeaves = 'off';
            end
            bagger.MergeLeaves = strcmpi(bagger.ModelParams.MergeLeaves,'on');
            if bagger.MergeLeaves
                warning(message('stats:TreeBagger:init:BadMergeLeavesValue'));
            end
                                    
            % Figure out logic for input params
            bagger.ComputeOOBPrediction = ...
                bagger.ComputeOOBPrediction || bagger.ComputeOOBVarImp;

            % Have enough observations for OOB error?
            N = size(bagger.PrivateX,1);
            if bagger.ComputeOOBPrediction && ~bagger.SampleWithReplacement ...
                    && N*(1-bagger.FBoot)<1
                error(message('stats:TreeBagger:init:NotEnoughOOBobservations'));
            end
                
            % Make Compact object
            if strcmp(method(1),'c')
                classnames = cellstr(bagger.ClassSummary.ClassNames);
            else
                classnames = {};
            end
            varnames = bagger.DataSummary.PredictorNames;
            if isnumeric(varnames)
                varnames = classreg.learning.internal.defaultPredictorNames(varnames);
            end
            bagger.Compact = CompactTreeBagger({},classnames,varnames);

            % Check number of variables to sample at random
            nvartosample = bagger.ModelParams.NVarToSample;
            if isempty(nvartosample)
                if     method(1)=='c'
                    nvartosample = ceil(sqrt(length(varnames)));
                elseif method(1)=='r'
                    nvartosample = ceil(length(varnames)/3);
                end
            end
            bagger.NVarToSample = nvartosample;
            bagger.ModelParams.NVarToSample = nvartosample;
            
            % Check tree leaf size
            minleaf = bagger.ModelParams.MinLeaf;
            if isempty(minleaf)
                if     method(1)=='c'
                    minleaf = 1;
                elseif method(1)=='r'
                    minleaf = 5;
                end
            end
            bagger.MinLeaf = minleaf;
            bagger.ModelParams.MinLeaf = minleaf;
            
            % Catch user attempts to set leaf size through decision tree
            % implementation
            minparent = bagger.ModelParams.MinParent;
            if ~isempty(minparent) || ~isempty(splitmin)
                error(message('stats:TreeBagger:init:MinparentSplitminNotAllowed'));
            end
            
            bagger.ModelParams = fillDefaultParams(bagger.ModelParams,...
                bagger.PrivateX,bagger.Y,bagger.W,bagger.DataSummary,bagger.ClassSummary);
 
            % Compute default scores
            if method(1)=='c'
                bagger.Compact.ClassProb = bagger.Prior;
                bagger.DefaultYfit = 'mostpopular';
            else
                bagger.DefaultYfit = dot(bagger.W,bagger.Y)/sum(bagger.W);
            end
        end

        function x = getXForPrediction(bagger,x)
            vrange = bagger.DataSummary.VariableRange;
            x = classreg.learning.internal.table2PredictMatrix(x,[],[],...
                vrange,bagger.DataSummary.CategoricalPredictors,bagger.PredictorNames);
        end
    end    
    
    
    methods
        function [varargout] = predict(bagger,X,varargin)
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
            %      'Trees'     Array of tree indices to use for computation of
            %                  responses.  Default is 'all'.
            %      'TreeWeights'  Array of NumTrees weights for weighting votes from the
            %                  specified trees.
            %      'UseInstanceForTree'  Logical matrix of size N-by-NumTrees
            %                  indicating which trees should be used to make
            %                  predictions for each observation. By default all trees
            %                  are used for all observations. If no trees are used for
            %                  an observation (the row of the matrix has no true
            %                  values), PREDICT returns the default prediction, most
            %                  popular class for classification or sample mean for
            %                  regression. The default prediction can be changed using
            %                  the DefaultYfit property.
            %
            %   See also TREEBAGGER, CLASSNAMES, ClassificationTree/predict,
            %   RegressionTree/predict, DEFAULTYFIT, SETDEFAULTYFIT, Method.

            [varargout{1:nargout}] = predict(bagger.Compact,X,varargin{:});
        end
        
        
        function [varargout] = quantilePredict(bagger,X,varargin)
            %QUANTILEPREDICT Predict response quantile.
            %   YFIT = QUANTILEPREDICT(B,X) computes medians of the predicted
            %   response of the trained ensemble B for predictors X. X should be a
            %   table if B was originally trained on a table, or a matrix if B was
            %   originally trained on a matrix. The numeric array YFIT has one
            %   prediction for each row of X.
            %
            %   [YFIT,YW] = QUANTILEPREDICT(B,X) also returns an Ntrain-by-N sparse
            %   matrix of response weights YW for Ntrain elements in B.Y and N
            %   observations (rows) in X. You can use this matrix for alternative
            %   estimation of the predicted quantiles.
            %
            %   QUANTILEPREDICT errors if B is a classification ensemble.
            %
            %   [...] = QUANTILEPREDICT(B,X,'PARAM1',val1,'PARAM2',val2,...) specifies
            %   optional parameter name/value pairs:
            %      'Quantile'  Vector of quantiles between 0 and 1. If you pass a
            %                  vector of Q quantiles, QUANTILEPREDICT returns an N-by-Q
            %                  matrix YFIT for N observations (rows) in X. Default: 0.5
            %      'Trees'     Array of tree indices to use for computation of
            %                  responses.  Default is 'all'.
            %      'TreeWeights'  Array of NumTrees weights for weighting votes from the
            %                  specified trees.
            %      'UseInstanceForTree'  Logical matrix of size N-by-NumTrees
            %                  indicating which trees should be used to make
            %                  predictions for each observation. By default all trees
            %                  are used for all observations. If no trees are used for
            %                  an observation (the row of the matrix has no true
            %                  values), QUANTILEPREDICT returns quantiles derived from
            %                  B.Y.
            %
            %   See also TREEBAGGER, RegressionTree/predict, Method.

            [varargin{:}] = convertStringsToChars(varargin{:});
            [varargout{1:nargout}] = quantilePredictCompact(bagger.Compact,X,...
                bagger.X,bagger.Y,bagger.W,...
                'inbagindices',bagger.InBagIndices,varargin{:});
        end
        
        function [varargout] = oobPredict(bagger,varargin)
            %OOBPREDICT Ensemble predictions for out-of-bag observations.
            %   Y = OOBPREDICT(B) computes predicted responses using the trained bagger
            %   B for out-of-bag observations in the training data.  The output has one
            %   prediction for each observation in the training data. The returned Y is
            %   a cell array of strings for classification and a numeric array for
            %   regression. For observations that are in bag for all trees, OOBPREDICT
            %   returns the default prediction, most popular class for classification
            %   or sample mean for regression.
            %
            %   Y = OOBPREDICT(B,'PARAM1',val1,'PARAM2',val2,...) specifies optional
            %   parameter name/value pairs:
            %      'Trees'     Array of tree indices to be used for computation of
            %                  responses.  Default is 'all'.
            %      'TreeWeights'  Array of NumTrees weights for weighting votes from the
            %                  specified trees.
            %
            % See also TREEBAGGER, COMPACTTREEBAGGER/PREDICT, OOBINDICES.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            % Check if OOB info was filled
            if ~bagger.ComputeOOBPrediction
                error(message('stats:TreeBagger:oobPredict:InvalidOperation'));
            end
            
            % Get OOB predictions
            [varargout{1:nargout}] = predict(bagger.Compact,bagger.X,...
                'useifort',bagger.OOBIndices,varargin{:});
        end
                
        function [varargout] = oobQuantilePredict(bagger,varargin)
            %OOBQUANTILEPREDICT Quantile predictions for out-of-bag observations.
            %   YFIT = OOBQUANTILEPREDICT(B) computes medians of the predicted
            %   response using the trained bagger B for out-of-bag observations in the
            %   training data. The numeric array YFIT has one prediction for each
            %   observation (row) in the training data.
            %
            %   [YFIT,YW] = OOBQUANTILEPREDICT(B) also returns an N-by-N sparse matrix
            %   of response weights YW for N observations in the training data. I-th
            %   column of this matrix specifies weights for I-th observation (row) in
            %   B.X. You can use this matrix for alternative estimation of the
            %   predicted quantiles.
            %
            %   OOBQUANTILEPREDICT errors if B is a classification ensemble.
            %
            %   [...] = OOBQUANTILEPREDICT(B,'PARAM1',val1,'PARAM2',val2,...) specifies
            %   optional parameter name/value pairs:
            %      'Quantile'  Vector of quantiles between 0 and 1. If you pass a
            %                  vector of Q quantiles, OOBQUANTILEPREDICT returns an N-by-Q
            %                  matrix YFIT for N observations (rows) in X. Default: 0.5
            %      'Trees'     Array of tree indices to use for computation of
            %                  responses.  Default is 'all'.
            %      'TreeWeights'  Array of NumTrees weights for weighting votes from the
            %                  specified trees.
            %
            %   See also TREEBAGGER, RegressionTree/predict, Method, predict, OOBIndices.

            [varargin{:}] = convertStringsToChars(varargin{:});
            % Check if OOB info was filled
            if ~bagger.ComputeOOBPrediction
                error(message('stats:TreeBagger:oobPredict:InvalidOperation'));
            end
            
            % Get OOB predictions
            [varargout{1:nargout}] = quantilePredict(bagger,bagger.X,...
                'useifort',bagger.OOBIndices,varargin{:});
        end        
        
        function err = oobError(bagger,varargin)
            %OOBERROR Out-of-bag error (misclassification probability or MSE).
            %   ERR = OOBERROR(B,) computes the misclassification probability (for
            %   classification trees) or mean squared error (MSE, for regression trees)
            %   for out-of-bag observations in the training data, using the trained
            %   bagger B. ERR is a vector with one error measure for each of the
            %   NumTrees trees in the ensemble B.
            %
            %   ERR = OOBERROR(B,'PARAM1',val1,'PARAM2',val2,...) specifies optional
            %   parameter name/value pairs:.
            %     'Mode'     String indicating how errors are computed. If set to
            %                'cumulative' (default), cumulative errors are computed and
            %                ERR is a vector of length NumTrees, where the 1st element
            %                gives error from tree 1, 2nd element gives error from
            %                trees 1:2 etc, up to 1:NumTrees. If set to 'individual',
            %                ERR is a vector of length NumTrees, where each element is
            %                an error from each tree in the ensemble. If set to
            %                'ensemble', ERR is a scalar showing the cumulative error
            %                for the entire ensemble.
            %     'Trees'    Vector of indices indicating what trees to include
            %                in this calculation. By default, this argument is set to
            %                'all' and the method uses all trees. If 'Trees' is a numeric
            %                vector, the method returns a vector of length NumTrees for
            %                'cumulative' and 'individual' modes, where NumTrees is the
            %                number of elements in the input vector, and a scalar for
            %                'ensemble' mode. For example, in the 'cumulative' mode,
            %                the first element gives error from trees(1), the second
            %                element gives error from trees(1:2) etc.
            %     'TreeWeights' Vector of tree weights. This vector must have the same
            %                length as the 'Trees' vector. OOBERROR uses these weights to
            %                combine output from the specified trees by taking a
            %                weighted average instead of the simple non-weighted
            %                majority vote. You cannot use this argument in the
            %                'individual' mode.
            %
            %   See also TREEBAGGER, COMPACTTREEBAGGER/ERROR.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            % Check if OOB info was filled
            if ~bagger.ComputeOOBPrediction
                error(message('stats:TreeBagger:oobError:InvalidOperation'));
            end

            % Get errors
            err = error(bagger.Compact,bagger.X,bagger.Y,...
                'weights',bagger.W,'useifort',bagger.OOBIndices,varargin{:});
        end
        
        function err = oobQuantileError(bagger,varargin)
            %OOBQUANTILEERROR Quantile loss.
            %   ERR = OOBQUANTILEERROR(B) returns half of the mean absolute deviation
            %   for out-of-bag observations in the training data, using the trained
            %   bagger B. ERR is a numeric scalar.
            %
            %   ERR = OOBQUANTILEERROR(B,'PARAM1',val1,'PARAM2',val2,...) specifies
            %   optional parameter name/value pairs.
            %     'Mode'     String indicating how errors are computed. If set to
            %                'cumulative', cumulative errors are computed and ERR is a
            %                vector of length NumTrees, where the 1st element gives
            %                error from tree 1, 2nd element gives error from trees 1:2
            %                etc, up to 1:NumTrees. If set to 'individual', ERR is a
            %                vector of length NumTrees, where each element is an error
            %                from each tree in the ensemble. If set to 'ensemble'
            %                (default), ERR is a scalar showing the error for the
            %                entire ensemble.
            %     'Quantile' Vector of quantiles between 0 and 1. If you pass a
            %                vector of Q quantiles, OOBQUANTILEERROR returns a
            %                row-vector ERR with Q elements for the 'ensemble' mode. If
            %                you set the mode to 'individual' or 'cumulative',
            %                OOBQUANTILEERROR returns a T-by-Q matrix for T trees and Q
            %                quantiles. Default: 0.5
            %     'Trees'    Vector of indices indicating what trees to include
            %                in this calculation. By default, this argument is set to
            %                'all' and OOBQUANTILEERROR uses all trees. If 'Trees' is a
            %                numeric vector, the method returns a vector of length
            %                NumTrees for 'cumulative' and 'individual' modes, where
            %                NumTrees is the number of elements in the input vector,
            %                and a scalar for 'ensemble' mode. For example, in the
            %                'cumulative' mode, the 1st element gives an error from
            %                trees(1), the 2nd element gives an error from trees(1:2)
            %                etc.
            %     'TreeWeights' Vector of tree weights. This vector must have the same
            %                length as the 'Trees' vector. OOBQUANTILEERROR uses these
            %                weights to combine output from the specified trees. You
            %                cannot use this argument in the 'individual' mode.
            %
            % See also TREEBAGGER, oobQuantilePredict, quantileError, NumTrees.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            % Check if OOB info was filled
            if ~bagger.ComputeOOBPrediction
                error(message('stats:TreeBagger:oobError:InvalidOperation'));
            end
            
            % Get errors
            err = quantileError(bagger,bagger.X,bagger.Y,...
                'weights',bagger.W,'useifort',bagger.OOBIndices,...
                varargin{:});           
        end
        
        function mar = oobMargin(bagger,varargin)
            %OOBMARGIN Out-of-bag margins.
            %   MAR = OOBMARGIN(B) computes an Nobs-by-NumTrees matrix of classification
            %   margins for out-of-bag observations in the training data, using the
            %   trained bagger B.
            %
            %   MAR = OOBMARGIN(B,'PARAM1',val1,'PARAM2',val2,...) specifies optional
            %   parameter name/value pairs.
            %
            %     'Mode'     String indicating how OOBMARGIN computes errors. If set to
            %                'cumulative' (default), the method computes cumulative
            %                margins and MAR is an Nobs-by-NumTrees matrix, where the 1st
            %                column gives margins from tree 1, 2nd column gives margins
            %                from trees 1:2 etc, up to 1:NumTrees. If set to 'individual',
            %                MAR is an Nobs-by-NumTrees matrix, where each column gives
            %                margins from each tree in the ensemble. If set to
            %                'ensemble', MAR is a single column of length Nobs showing
            %                the cumulative margins for the entire ensemble.
            %     'Trees'    Vector of indices indicating what trees to include
            %                in this calculation. By default, this argument is set to
            %                'all' and the method uses all trees. If 'Trees' is a numeric
            %                vector, the method returns an Nobs-by-NumTrees matrix for
            %                'cumulative' and 'individual' modes, where NumTrees is the
            %                number of elements in the input vector, and a single
            %                column for 'ensemble' mode. For example, in the
            %                'cumulative' mode, the first column gives margins from
            %                trees(1), the second column gives margins from
            %                trees(1:2) etc.
            %     'TreeWeights' Vector of tree weights. This vector must have the same
            %                length as the 'Trees' vector. OOBMARGIN uses these weights
            %                to combine output from the specified trees by taking a
            %                weighted average instead of the simple non-weighted
            %                majority vote. You cannot use this argument in the
            %                'individual' mode.
            %
            %   See also TREEBAGGER, COMPACTTREEBAGGER/MARGIN.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            % Check if OOB info was filled
            if ~bagger.ComputeOOBPrediction
                error(message('stats:TreeBagger:oobMargin:InvalidOperation'));
            end
            
            % Get margins
            mar = margin(bagger.Compact,bagger.X,bagger.Y,...
                'useifort',bagger.OOBIndices,varargin{:});
        end
        
        function mar = oobMeanMargin(bagger,varargin)
            %OOBMEANMARGIN Out-of-bag mean margins.
            %   MAR = OOBMEANMARGIN(B) computes average classification margins for
            %   out-of-bag observations in the training data, using the trained bagger
            %   B. OOBMEANMARGIN averages the margins over all out-of-bag observations.
            %   MAR is a row-vector of length NumTrees, where NumTrees is the number of
            %   trees in the ensemble.
            %
            %   MAR = OOBMEANMARGIN(B,'PARAM1',val1,'PARAM2',val2,...) specifies
            %   optional parameter name/value pairs:
            %
            %     'Mode'     String indicating how the method computes
            %                errors. If set to 'cumulative' (default), OOBMEANMARGIN
            %                computes cumulative errors and MAR is a vector of length
            %                NumTrees, where the 1st element gives mean margin from tree
            %                1, 2nd element gives mean margin from trees 1:2 etc, up to
            %                1:NumTrees. If set to 'individual', MAR is a vector of
            %                length NumTrees, where each element is a mean margin from
            %                each tree in the ensemble. If set to 'ensemble', MAR is a
            %                scalar showing the cumulative mean margin for the entire
            %                ensemble.
            %     'Trees'    Vector of indices indicating what trees to
            %                include in this calculation. By default, this argument is
            %                set to 'all' and the method uses all trees. If 'Trees' is
            %                a numeric vector, the method returns a vector of length
            %                NumTrees for 'cumulative' and 'individual' modes, where
            %                NumTrees is the number of elements in the input vector, and
            %                a scalar for 'ensemble' mode. For example, in the
            %                'cumulative' mode, the first element gives mean margin
            %                from trees(1), the second element gives mean margin from
            %                trees(1:2) etc.
            %     'TreeWeights' Vector of tree weights. This vector must
            %                have the same length as the 'Trees' vector. OOBMEANMARGIN
            %                uses these weights to combine output from the specified
            %                trees by taking a weighted average instead of the simple
            %                non-weighted majority vote. You cannot use this argument
            %                in the 'individual' mode.
            %
            % See also TREEBAGGER, COMPACTTREEBAGGER/MEANMARGIN.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            % Check if OOB info was filled
            if ~bagger.ComputeOOBPrediction
                error(message('stats:TreeBagger:oobMeanMargin:InvalidOperation'));
            end
            
            % Get margins
            mar = meanMargin(bagger.Compact,bagger.X,bagger.Y,...
                'weights',bagger.W,'useifort',bagger.OOBIndices,varargin{:});
        end
        
        function err = error(bagger,X,varargin)
            %ERROR Error (misclassification probability or MSE).
            %   ERR = ERROR(B,X,Y) computes the misclassification probability (for
            %   classification trees) or mean squared error (MSE, for regression trees)
            %   for predictors X given true response Y. X must be a table if B was
            %   originally trained on a table, or a matrix if B was originally trained
            %   on a matrix. For classification, Y can be either a numeric vector, cell
            %   array of strings, categorical vector, or logical vector. For
            %   regression, Y must be a numeric vector. Y can be omitted if X is a
            %   table that includes the response variable. For classification and MSE
            %   regression, ERR is a vector with one error measure for each of the
            %   NumTrees trees in the ensemble B.
            %
            %   ERR = ERROR(B,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
            %   parameter name/value pairs.
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
            % See also TREEBAGGER, COMPACTTREEBAGGER/ERROR, NumTrees.

            % Get Y from input or table as required 
            [y,varargin] = inferResponse(bagger.Compact,X,varargin{:});
            
            % Get errors
            err = error(bagger.Compact,X,y,varargin{:});
        end
        
        function err = quantileError(bagger,X,varargin)
            %QUANTILEERROR Quantile loss.
            %   ERR = QUANTILEERROR(B,X,Y) returns half of the mean absolute deviation
            %   for predictors X given true response Y. X must be a table if B was
            %   originally trained on a table, or a matrix if B was originally trained
            %   on a matrix. Y must be a numeric vector. Y can be omitted if X is a
            %   table that includes the response variable. ERR is a numeric scalar.
            %
            %   ERR = QUANTILEERROR(B,X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies
            %   optional parameter name/value pairs.
            %     'Mode'     String indicating how errors are computed. If set to
            %                'cumulative', cumulative errors are computed and ERR is a
            %                vector of length NumTrees, where the 1st element gives
            %                error from tree 1, 2nd element gives error from trees 1:2
            %                etc, up to 1:NumTrees. If set to 'individual', ERR is a
            %                vector of length NumTrees, where each element is an error
            %                from each tree in the ensemble. If set to 'ensemble'
            %                (default), ERR is a scalar showing the cumulative error
            %                for the entire ensemble.
            %     'Weights'  Vector of observation weights to use for error
            %                averaging. By default the weight of every observation
            %                is set to 1. The length of this vector must be equal
            %                to the number of rows in X.
            %     'Quantile' Vector of quantiles between 0 and 1. If you pass a
            %                vector of Q quantiles, QUANTILEERROR returns a row-vector
            %                ERR with Q elements for the 'ensemble' mode. If you set
            %                the mode to 'individual' or 'cumulative', QUANTILEERROR
            %                returns a T-by-Q matrix for T trees and Q quantiles.
            %                Default: 0.5
            %     'Trees'    Vector of indices indicating what trees to include
            %                in this calculation. By default, this argument is set to
            %                'all' and QUANTILEERROR uses all trees. If 'Trees' is a
            %                numeric vector, the method returns a vector of length
            %                NumTrees for 'cumulative' and 'individual' modes, where
            %                NumTrees is the number of elements in the input vector,
            %                and a scalar for 'ensemble' mode. For example, in the
            %                'cumulative' mode, the 1st element gives an error from
            %                trees(1), the 2nd element gives an error from trees(1:2)
            %                etc.
            %     'TreeWeights' Vector of tree weights. This vector must have the same
            %                length as the 'Trees' vector. QUANTILEERROR uses these
            %                weights to combine output from the specified trees. You
            %                cannot use this argument in the 'individual' mode.
            %     'UseInstanceForTree' Logical matrix of size N-by-NumTrees
            %                indicating which trees QUANTILEERROR should use to make
            %                predictions for each observation. By default the method
            %                uses all trees for all observations.
            %
            % See also TREEBAGGER, quantilePredict, NumTrees.

            [varargin{:}] = convertStringsToChars(varargin{:});
            % Get Y from input or table as required 
            [y,varargin] = inferResponse(bagger.Compact,X,varargin{:});

            % Compute error
            err = quantileErrorCompact(bagger.Compact,X,y,...
                bagger.X,bagger.Y,bagger.W,...
                'inbagindices',bagger.InBagIndices,varargin{:});
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
            %   parameter name/value pairs:
            %
            %     'Mode'     String indicating how the method computes errors. If set
            %                to 'cumulative' (default), the method computes cumulative
            %                margins and MAR is an Nobs-by-NumTrees matrix, where the 1st
            %                column gives margins from tree 1, 2nd column gives margins
            %                from trees 1:2 etc, up to 1:NumTrees. If set to
            %                'individual', MAR is an Nobs-by-NumTrees matrix, where each
            %                column gives margins from each tree in the ensemble. If
            %                set to 'ensemble', MAR is a single column of length Nobs
            %                showing the cumulative margins for the entire ensemble.
            %     'Trees'    Vector of indices indicating what trees to include
            %                in this calculation. By default, this argument is set to
            %                'all' and the method uses all trees. If 'Trees' is a
            %                numeric vector, the method returns a vector of length
            %                NumTrees for 'cumulative' and 'individual' modes, where
            %                NumTrees is the number of elements in the input vector, and
            %                a scalar for 'ensemble' mode. For example, in the
            %                'cumulative' mode, the first element gives error from tree
            %                trees(1), the second element gives error from trees
            %                trees(1:2) etc.
            %     'TreeWeights' Vector of tree weights. This vector must have the same
            %                length as the 'Trees' vector. The method uses these
            %                weights to combine output from the specified trees by
            %                taking a weighted average instead of the simple
            %                non-weighted majority vote. You cannot use this argument
            %                in the 'individual' mode.
            %     'UseInstanceForTree' Logical matrix of size Nobs-by-NumTrees 
            %                indicating which trees to use to make
            %                predictions for each observation.  By default
            %                the method uses all trees for all
            %                observations.
            %
            % See also TREEBAGGER, COMPACTTREEBAGGER/MARGIN.

            % Get Y from input or table as required 
            [y,varargin] = inferResponse(bagger.Compact,X,varargin{:});
            
            % Get margins
            mar = margin(bagger.Compact,X,y,varargin{:});
        end
        
        function [varargout] = meanMargin(bagger,X,varargin)
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
            %   optional parameter name/value pairs:
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
            %                'all' and the method uses all trees. If 'Trees' is a
            %                numeric vector, the method returns a vector of length
            %                NumTrees for 'cumulative' and 'individual' modes, where
            %                NumTrees is the number of elements in the input vector, and
            %                a scalar for 'ensemble' mode. For example, in the
            %                'cumulative' mode, the first element gives error from tree
            %                trees(1), the second element gives error from trees
            %                trees(1:2) etc.
            %     'TreeWeights' Vector of tree weights. This vector must have the same
            %                length as the 'Trees' vector. The method uses these
            %                weights to combine output from the specified trees by
            %                taking a weighted average instead of the simple
            %                non-weighted majority vote. You cannot use this argument
            %                in the 'individual' mode.
            %     'UseInstanceForTree' Logical matrix of size Nobs-by-NumTrees 
            %                indicating which trees to use to make
            %                predictions for each observation.  By default
            %                the method uses all trees for all
            %                observations.
            %
            %   See also TREEBAGGER, COMPACTTREEBAGGER/MEANMARGIN.

            % Get Y from input or table as required 
            [y,varargin] = inferResponse(bagger.Compact,X,varargin{:});
            
            [varargout{1:nargout}] = meanMargin(bagger.Compact,X,y,varargin{:});
        end
        
        function [varargout] = mdsprox(bagger,varargin)
            %MDSPROX Multidimensional scaling of proximity matrix.
            %   [S,E] = MDSPROX(B) returns scaled coordinates S and eigenvalues E for
            %   the proximity matrix in the ensemble B.  The proximity matrix must be
            %   created by an earlier call to FILLPROX(B).
            %
            %   [S,E] = MDSPROX(B,'PARAM1',val1,'PARAM2',val2,...) specifies optional
            %   parameter name/value pairs:
            %
            %    'Keep'   Array of indices of observations in the training data to
            %             use for multidimensional scaling. By default, this argument
            %             is set to 'all'. If you provide numeric or logical indices,
            %             the method uses only the subset of the training data specified
            %             by these indices to compute the scaled coordinates and
            %             eigenvalues.
            %    'Colors' If you supply this argument, MDSPROX makes overlaid scatter
            %             plots of two scaled coordinates using specified colors for
            %             different classes. You must supply the colors as a string
            %             with one character for each color.  If there are more classes
            %             in the data than characters in the supplied string, MDSPROX
            %             only plots the first C classes, where C is the length of the
            %             string. For regression or if you do not provide the vector of
            %             true class labels, the method uses the first color for all
            %             observations in X.
            %    'MDSCoordinates' Indices of the two or three scaled coordinates to be
            %             plotted. By default, MDSPROX makes a scatter plot of the 1st
            %             and 2nd scaled coordinates which correspond to the two
            %             largest eigenvalues.  You can specify any other two or three
            %             indices not exceeding the dimensionality of the scaled data.
            %             This argument has no effect unless you also supply the
            %             'Colors' argument.
            %
            %   See also TREEBAGGER, COMPACTTREEBAGGER/MDSPROX, CMDSCALE,
            %   FILLPROX.
            
            [varargin{:}] = convertStringsToChars(varargin{:});
            [varargout{1:nargout}] = mdsProx(bagger,varargin{:});
        end
        
        
        function [pd,varargout] = partialDependence(bagger,features,varargin)
        %PARTIALDEPENDENCE Partial Dependence for one and two variables
        %   PD = partialDependence(MODEL,VARS) is for a fitted regression
        %   TreeBagger. It computes the partial dependence PD between the predictor
        %   variables listed in VARS and the response of the model, by averaging
        %   over the data used in training the model. You can specify one or two
        %   variables using VARS. VARS is either column indices of predictors or
        %   names of predictors. The names must match the entries in
        %   MODEL.PredictorNames.
        %
        %   PD = partialDependence(MODEL,VARS,LABELS) is for a fitted classification
        %   TreeBagger. You can specify one or multiple classes using LABELS. The
        %   values in LABELS must match and have the same data type as the class
        %   names in MODEL.ClassNames. If there are multiple values in both VARS
        %   and LABELS, then PD is a three-dimensional array. PD is the partial
        %   dependence between the predictor variables listed in VARS and the
        %   scores for the classes specified in LABELS.
        %
        %   PD = partialDependence(...,DATA) specifies the data to be used for 
        %   averaging. DATA is a matrix or table of data to be used in place of the
        %   data used in fitting the model.
        %
        %   PD = partialDependence(...,Name,Value) specifies additional options
        %   using one or more name-value pair arguments:
        %
        %      'NumObservationsToSample':   An integer specifying the number
        %                                   of rows to sample at random from the
        %                                   dataset (either the data specified by
        %                                   argument DATA or the data used to train
        %                                   MODEL). The default is to use all rows.
        %
        %      'QueryPoints':               The points at which to calculate
        %                                   the partial dependence. For the case of
        %                                   a single predictor, the value of
        %                                   'QueryPoints' must be a column vector.
        %                                   For the case of two predictors, the
        %                                   value of 'QueryPoints' must be a 1x2
        %                                   cell array containing a separate vector
        %                                   for each predictor. The default is to
        %                                   use 100 points equally spaced across
        %                                   the range of the predictor.
        %
        %      'UseParallel'                Can be true to instruct the method to
        %                                   perform the averaging calculations in
        %                                   parallel (using parfor), or false
        %                                   (default) to instruct the method not to
        %                                   use parallel computation.
        %
        %   [PD,X,Y] = partialDependence(...) also returns the arrays X and Y
        %   containing the query point values of the first and second predictor in
        %   VARS, respectively.
        % 
        %   Examples:
        %       % Partial Dependence for bagged decision trees
        %       load census1994
        %       X = adultdata(:,{'age','workClass','education_num','marital_status','race',...
        %           'sex','capital_gain','capital_loss','hours_per_week','salary'});
        %       f = TreeBagger(50,X,'salary');
        % 
        %       pd = partialDependence(f,'age','<=50K');
        % 
        %       % Obtain values x of the variable 'age'
        %       [pd,x] = partialDependence(f,'age',{'<=50K','>50K'});
        %
        %       % plot partial dependence of the first and second classes scores on
        %       % 'age'
        %       plot(x,pd)
        % 
        %       % Obtain values x and y of variables 'age' and 'workClass'
        %       [pd,x,y] = partialDependence(f,{'age','workClass'},{'<=50K','>50K'});
        % 
        %       % plot partial dependence for class '<=50K' on variables
        %       % 'age' and 'workClass'
        %       surf(x,y,pd(:,:,1))
        %
        %       % plot partial dependence for class '>50K' on variables
        %       % 'age' and 'workClass'
        %       surf(x,y,pd(:,:,2))
        % 
        %       % With additional Data
        %       pd = partialDependence(f,'age','<=50K',X);
        %       pd = partialDependence(f,1,'<=50K',X);
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
        %       pd = partialDependence(f,1,'<=50K','QueryPoints',xi);
        % 
        %       % Provide alternative query points
        %       xi = linspace(min(X{:,'age'}),max(X{:,'age'}))';
        %       pd = partialDependence(f,1,'<=50K','QueryPoints',xi);
        % 
        %       xi = cell(1,2);
        %       xi{1} = linspace(min(X{:,1}),max(X{:,1}))';
        %       xi{2} = linspace(min(X{:,9}),max(X{:,9}))';
        %       pd = partialDependence(f,[1,9],'<=50K','QueryPoints',xi);
        
        %-------Check number of inputs----
        narginchk(2,10);
        features = convertStringsToChars(features);
        [varargin{:}] = convertStringsToChars(varargin{:});
               
        if strcmp(bagger.Method,'regression')
            % class label is not an input for regression models
            % Check inputs with inputParser. This step ensures that the third
            % argument is either a Name-Value pair or data, no other strings/char
            % array allowed.
            p = inputParser;
            addRequired(p,'Model');
            addRequired(p,'Var');
            addOptional(p,'Data',bagger.X); % Default - training data
            addParameter(p,'NumObservationsToSample',0);
            addParameter(p,'QueryPoints',[]);
            addParameter(p,'UseParallel',false); 
            parse(p,bagger,features,varargin{:});
            data = p.Results.Data;
            
            %------Parse Data-----------------
            if(nargin>2 && ~ischar(varargin{1}))
                % Pass everything but the first argument(Data)to compact method
                varargin = varargin(2:end);
            end  
            % Classification flag is false and class label is empty
            classif = false;
            label = [];
         
        % Classification tree bagger
        else
            % class label is an input for regression models 
            if nargin < 3
                error(message('MATLAB:narginchk:notEnoughInputs'));
            else
            label = varargin{1};
           
            % Check inputs with inputParser. This step ensures that the
            % fourth argument is either a Name-Value pair or data, no other strings/char
            % array allowed.
            p = inputParser;        
            addRequired(p,'Model');
            addRequired(p,'Var');
            addRequired(p,'Label');
            addOptional(p,'Data',bagger.X); % Default - training data
            addParameter(p,'NumObservationsToSample',0);
            addParameter(p,'QueryPoints',[]);
            addParameter(p,'UseParallel',false); 
            parse(p,bagger,features,varargin{:});
            data = p.Results.Data;
            
            %------Parse Data-----------------
            varargin(1) = [];  % remove label
            if (~isempty(varargin) && ~ischar(varargin{1}))
                % Pass everything but the first argument(Data)to compact method
                varargin(1) = [];  % remove data
            end  
            
            % Classification flag is true
            classif = true;
            end
        end
        
        % Call the function from classreg.learning.impl
        [pd,x,y] = classreg.learning.impl.partialDependence(bagger,features,...
            label,data,classif,varargin{:});
        varargout{1} = x;
        varargout{2} = y;
        end
        
        function [AX] = plotPartialDependence(bagger,features,varargin)
        %PLOTPARTIALDEPENDENCE Partial Dependence Plot for 1-D or 2-D visualization
        %   plotPartialDependence(MODEL,VARS) is for a fitted regression
        %   TreeBagger. It plots the partial dependence PD between the predictor
        %   variables listed in VARS and the response of the model, by averaging
        %   over the data used in training the model. You can specify one or two
        %   variables using VARS. VARS is either column indices of predictors or
        %   names of predictors. The names must match the entries in
        %   MODEL.PredictorNames.
        %
        %   plotPartialDependence(MODEL,VARS,LABELS) is for a fitted classification
        %   TreeBagger. You can specify one or multiple classes using LABELS. The
        %   values in LABELS must match and have the same data type as the class
        %   names in MODEL.ClassNames. The plotted value PD is the partial
        %   dependence between the predictor variables listed in VARS and the
        %   scores for the classes specified in LABELS. If there are multiple
        %   values in VARS, then LABELS must be a single label.
        % 
        %   plotPartialDependence(...,DATA) specifies the data to be used for
        %   averaging. DATA is a matrix or table of data to be used in place of the
        %   data used in fitting the model.
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
        %   'NumObservationsToSample'	An integer specifying the number
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
        %   'Parent'                    Instructs the method to create a partial
        %                               dependence plot using the axes with handle
        %                               specified by the value of this parameter.
        %
        %   AX = plotPartialDependence(...) returns a handle AX to the axes of the
        %   plot.
        %
        %   Examples:
        %       % Partial Dependence Plot for bagged decision trees
        %       load census1994
        %       X = adultdata(:,{'age','workClass','education_num','marital_status','race',...
        %           'sex','capital_gain','capital_loss','hours_per_week','salary'});
        %       f = TreeBagger(50,X,'salary');
        %
        %       plotPartialDependence(f,'age','<=50K');
        %       plotPartialDependence(f,{'age','workClass'},'>50K');
        %       plotPartialDependence(f,'age',help t)
        %
        %       % Retrieve optional output Axes handle
        %       ax = plotPartialDependence(f,'age','<=50K');
        %
        %       % Specify alternative data
        %       plotPartialDependence(f,'age','<=50K',X);
        %       plotPartialDependence(f,1,'<=50K',X);
        %
        %       % Specify optional name-value pair arguments
        %       plotPartialDependence(f,1,'<=50K',X,'NumObservationsToSample',100);
        %       plotPartialDependence(f,1,'<=50K',X,'UseParallel',true);
        %       plotPartialDependence(f,1,'<=50K',X,'UseParallel',true,'Conditional','none');
        %
        %       % Plot the Individual Conditional Expectation
        %       plotPartialDependence(f,1,'<=50K',X,'Conditional','absolute');
        %
        %       % Provide alternative query points
        %       xi = linspace(min(X{:,'age'}),min(X{:,'age'}))';
        %       plotPartialDependence(f,1,'<=50K','QueryPoints',xi);
        %
        %       % Provide alternative query points for one predictor
        %       xi = linspace(min(X{:,'age'}),max(X{:,'age'}))';
        %       plotPartialDependence(f,1,'<=50K','QueryPoints',xi);
        %       % Provide alternative query points for two predictors
        %
        %       xi = cell(1,2);
        %       xi{1} = linspace(min(X{:,1}),max(X{:,1}))';
        %       xi{2} = linspace(min(X{:,9}),max(X{:,9}))';
        %       plotPartialDependence(f,[1,9],'<=50K','QueryPoints',xi);
        
        %-------Check number of inputs----
        narginchk(2,14);
        features = convertStringsToChars(features);
        [varargin{:}] = convertStringsToChars(varargin{:});
               
        if strcmp(bagger.Method,'regression')
            % class label is not an input for regression models
            if nargin > 13
                error(message('MATLAB:narginchk:tooManyInputs'));
            else
            
            % Check inputs with inputParser. This step ensures that the third
            % argument is either a Name-Value pair or data, no other strings/char
            % array allowed.
            p = inputParser;
            addRequired(p,'Model');
            addRequired(p,'Var');
            addOptional(p,'Data',bagger.X); % Default - training data
            addParameter(p,'Conditional',{'none','absolute','centered'});
            addParameter(p,'NumObservationsToSample',0);
            addParameter(p,'QueryPoints',[]);
            addParameter(p,'UseParallel',false);
            addParameter(p,'ParentAxisHandle',[]);
            parse(p,bagger,features,varargin{:});
            data = p.Results.Data;
            
            %------Parse Data-----------------
            if(nargin>2 && ~ischar(varargin{1}))
                % Pass everything but the first argument(Data)to compact method
                varargin = varargin(2:end);
            end  
            % Classification flag is false and class label is empty
            classif = false;
            labels = [];
            end
        
        % Classification tree bagger
        else
            % class label is an input for regression models 
            if nargin < 3
                error(message('MATLAB:narginchk:notEnoughInputs'));
            else
            labels = varargin{1};
            
            % Check inputs with inputParser. This step ensures that the
            % fourth argument is either a Name-Value pair or data, no other strings/char
            % array allowed.
            p = inputParser;        
            addRequired(p,'Model');
            addRequired(p,'Var');
            addRequired(p,'Label');
            addOptional(p,'Data',bagger.X); % Default - training data
            addParameter(p,'Conditional',{'none','absolute','centered'});
            addParameter(p,'NumObservationsToSample',0);
            addParameter(p,'QueryPoints',[]);
            addParameter(p,'UseParallel',false); 
            addParameter(p,'ParentAxisHandle',[]);
            parse(p,bagger,features,varargin{:});
            data = p.Results.Data;
            
            varargin = varargin(2:end);
            %------Parse Data-----------------
            if(nargin>3 && ~ischar(varargin{1}))
                % Pass everything but the first argument(Data)to compact method
                varargin = varargin(2:end);
            end  
            
            % Classification flag is true
            classif = true;
            end
        end
        
        % Call the function from classreg.learning.impl
        ax = classreg.learning.impl.plotPartialDependence(bagger,features,...
            labels,data,classif,varargin{:});
        if(nargout > 0)
            AX = ax;
        end 
        end
    end
    
    methods(Hidden=true)
        function [varargout] = mdsProx(bagger,varargin)
            % Check if proximities have been computed
            if isempty(bagger.PrivProx)
                error(message('stats:TreeBagger:mdsProx:InvalidProperty'));
            end
            
            % Process inputs
            args = {'keep'};
            defs = { 'all'};
            [keep,~,compactArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            % Check indices to be kept
            N = size(bagger.PrivateX,1);
            if strcmpi(keep,'all')
                keep = 1:N;
            end
            if ~isnumeric(keep) && (~islogical(keep) || length(keep)~=N)
                error(message('stats:TreeBagger:mdsProx:InvalidInput'));
            end
            
            % Get scaled coordinates
            [varargout{1:nargout}] = mdsProx(bagger.Compact,...
                bagger.Proximity(keep,keep),'labels',bagger.Y(keep),...
                'data','proximity',compactArgs{:});
        end
        
        function bagger = fillProximities(bagger,varargin)
            bagger.PrivProx = flatprox(bagger.Compact,bagger.X,varargin{:});
        end
    end
end

function throwUndefinedError()
error(message('stats:TreeBagger:UndefinedFunction'));
end

function nm = normalizeMean1(A)
% Init
nm = zeros(1,size(A,2));

% Get mean and stdev
m = mean(A,1);
s = std(A,1,1);

% Get deltas
above0 = s>0 | m>0;
nm(above0) = m(above0)./s(above0);
end

% function idx = plainSample(s,N,fsample,replace)
% % Get number of instances to draw
% Nsample = ceil(N*fsample);
% if Nsample==0
%     error(message('stats:TreeBagger:plainSample:NoData'))
% end
% 
% % Sample
% if isempty(s)
%     idx = randsample(N,Nsample,replace);
% else
%     idx = randsample(s,N,Nsample,replace);
% end
% end

function idx = weightedSample(s,w,fsample,replace)
% Get number of instances to draw
N = numel(w);
Nsample = ceil(N*fsample);
% Nsample will be always >=1 because fsample is required to be within (0,1]
% and N must be >=1

% Sample
if isempty(s)
    idx = datasample((1:N)',Nsample,'replace',replace,'weights',w);
else
    idx = datasample(s,(1:N)',Nsample,'replace',replace,'weights',w);
end
end

function [trees, ...
    slicedOOBIndices, ...
    slicedInBagIndices, ...
    slicedPrivOOBPermutedVarDeltaError, ...
    slicedPrivOOBPermutedVarDeltaMeanMargin, ...
    slicedPrivOOBPermutedVarCountRaiseMargin] = ...
    localGrowTrees(NTrees, useParallel, RNGscheme, args)

doclassregtree            = args{1};            
fboot                     = args{2};
sampleWithReplacement     = args{3};
x                         = args{4};
y                         = args{5};
w                         = args{6};
method                    = args{7};
dataSummary               = args{8};
classSummary              = args{9};
modelParams               = args{10};
computeOOBPrediction      = args{11};
computeOOBVarImp          = args{12};
nprint                    = args{13};
usePool                   = args{14};
xorig                     = args{15};

% Retrieve parameters
prune = modelParams.Prune;
nvartosample = modelParams.NVarToSample;
minleaf = modelParams.MinLeaf;
maxsplits = modelParams.MaxSplits;
mergeleaves = modelParams.MergeLeaves;
nsurrogate = modelParams.NSurrogate;
qetoler = modelParams.QEToler;
splitcrit = modelParams.SplitCriterion;
catpred = dataSummary.CategoricalPredictors;
maxcat = modelParams.MaxCat;
algcat = modelParams.AlgCat;
prunecrit = modelParams.PruneCriterion;
predictorsel = modelParams.PredictorSelection;
usechisq = modelParams.UseChisqTest;

if strcmp(method(1),'c')
    localY = grp2idx(classreg.learning.internal.ClassLabel(y));
end

% Preallocate storage for temporary variables
N = size(x,1);

% Initialize output
slicedOOBIndices                              = [];
slicedInBagIndices                            = [];
slicedPrivOOBPermutedVarDeltaError            = [];
slicedPrivOOBPermutedVarDeltaMeanMargin       = [];
slicedPrivOOBPermutedVarCountRaiseMargin      = [];

% Grow the trees
if computeOOBPrediction
    if computeOOBVarImp
        if method(1)=='c'
            [trees, ...
                slicedOOBIndices, ...
                slicedInBagIndices, ...
                slicedPrivOOBPermutedVarDeltaError, ...
                slicedPrivOOBPermutedVarDeltaMeanMargin, ...
                slicedPrivOOBPermutedVarCountRaiseMargin] = ...
                internal.stats.parallel.smartForSliceout(NTrees, @loopBody, useParallel, RNGscheme);
        else
            [trees, ...
                slicedOOBIndices, ...
                slicedInBagIndices, ...
                slicedPrivOOBPermutedVarDeltaError] = ...
                internal.stats.parallel.smartForSliceout(NTrees, @loopBody, useParallel, RNGscheme);
        end
    else
        if method(1)=='r'
            [trees, slicedOOBIndices, slicedInBagIndices] = ...
                internal.stats.parallel.smartForSliceout(NTrees, @loopBody, useParallel, RNGscheme);
        else
            [trees, slicedOOBIndices] = ...
                internal.stats.parallel.smartForSliceout(NTrees, @loopBody, useParallel, RNGscheme);
        end
    end    
elseif method(1)=='r'
    [trees, ~, slicedInBagIndices] = ...
        internal.stats.parallel.smartForSliceout(NTrees, @loopBody, useParallel, RNGscheme);
else
    trees = internal.stats.parallel.smartForSliceout(NTrees, @loopBody, useParallel, RNGscheme);
end

% Make sure vectors have correct orientation
if N==1
    slicedOOBIndices = slicedOOBIndices(:)';
end

if NTrees==1
    slicedPrivOOBPermutedVarDeltaError = slicedPrivOOBPermutedVarDeltaError(:)';
    slicedPrivOOBPermutedVarDeltaMeanMargin = slicedPrivOOBPermutedVarDeltaMeanMargin(:)';
    slicedPrivOOBPermutedVarCountRaiseMargin = slicedPrivOOBPermutedVarCountRaiseMargin(:)';
end


    % --------- Nested within localGrowTrees ----------
    %
    % Grow one tree and optional metrics. This is the body of iterative loops.
    function [slicedTree, ...
            slicedOOBIndices, ...
            slicedInBagIndices, ...
            slicedPrivOOBPermutedVarDeltaError, ...
            slicedPrivOOBPermutedVarDeltaMeanMargin, ...
            slicedPrivOOBPermutedVarCountRaiseMargin] = loopBody(~,s)
        
        slicedOOBIndices                              = [];
        slicedInBagIndices                            = [];
        slicedPrivOOBPermutedVarDeltaError            = [];
        slicedPrivOOBPermutedVarDeltaMeanMargin       = [];
        slicedPrivOOBPermutedVarCountRaiseMargin      = [];
        
        if isempty(s)
            s = RandStream.getGlobalStream;
        end
        
        % Draw instances for training
        %idxtrain = plainSample(s,N,fboot,sampleWithReplacement);
        idxtrain = weightedSample(s,w,fboot,sampleWithReplacement);
        
        if isempty(classSummary)
            cost = [];
            classNames = {};
        else
            % The cost matrix is reset to ones(K)-eye(K) by
            % classreg.learning.internal.adjustPrior. It is passed to
            % TreeImpl because TreeImpl requires a cost matrix at input.
            cost = classSummary.Cost;
            classNames = cellstr(classSummary.ClassNames);
        end
        
        varnames = dataSummary.PredictorNames;
        if isnumeric(varnames)
            varnames = classreg.learning.internal.defaultPredictorNames(varnames);
        end

        % Grow the tree
        if doclassregtree
            % classregtree has been removed. Use fitctree or fitrtree
            % instead.
            error(message('stats:obsolete:ReplaceThisWithObjectUsingOptions','classregtree','TreeBagger trained on','fitctree','fitrtree'));
        else
             if strcmp(method(1),'c')
                 impl = classreg.learning.impl.TreeImpl.makeFromData(...
                     x,localY,ones(N,1)/N,idxtrain,true,...
                     catpred,splitcrit,minleaf,2*minleaf,maxsplits,nvartosample,...
                     nsurrogate,maxcat,algcat,cost,0,...
                     predictorsel,usechisq,s);
                 if strcmp(prune,'on')
                     impl = impl.prune('cost',cost,'criterion',prunecrit);
                 end
                 if strcmp(mergeleaves,'on')
                     impl = impl.prune('cost',cost,'criterion',prunecrit,'level',0);
                 end
                 tree = classreg.learning.classif.CompactClassificationTree(...
                     dataSummary,classSummary,...
                     @classreg.learning.transform.identity,'probability');
                 tree.Impl = impl;
             else
                 impl = classreg.learning.impl.TreeImpl.makeFromData(...
                     x,y,ones(N,1)/N,idxtrain,false,...
                     catpred,splitcrit,minleaf,2*minleaf,maxsplits,nvartosample,...
                     nsurrogate,0,'',[],qetoler,...
                     predictorsel,usechisq,s);
                 if strcmp(prune,'on')
                     % Only MSE is allowed for regression.
                     % computePruneInfo uses it by default =>
                     % no need to pass prunecrit in.
                     impl = impl.prune();
                 end
                 if strcmp(mergeleaves,'on')
                     impl = impl.prune('level',0);
                 end
                 tree = classreg.learning.regr.CompactRegressionTree(...
                     dataSummary,@classreg.learning.transform.identity);
                 tree.Impl = impl;
            end
        end
        
        %
        % Optional stuff
        %
        
        % Out of bag estimates
        if computeOOBPrediction
            % Record indices of OOB instances
            oobtf = true(N,1);
            oobtf(idxtrain) = false;
            slicedOOBIndices = oobtf;  % assign to bagger.OOBIndices
            % outside the parfor loop
            
            % Make a compact object with a single tree for predictions
            localcompact = CompactTreeBagger({tree},classNames,varnames);
            
            % Update scores for feature importance.
            % Here, we pass observation weights to oobPermVarUpdate, even
            % in the presence of weighted sampling. The rationale is this:
            % If we passed uniform weights, OOB estimate of feature
            % importance would be driven by improbable observations
            % (observation with small weights) which are not representative
            % of the sampling distribution.
            if computeOOBVarImp
                % Collect the quantities in the loop.
                % Assign to TreeBagger class variables (ie, Properties)
                % on exit from the loop.
                if method(1)=='c'
                    [slicedPrivOOBPermutedVarDeltaError, ...
                        slicedPrivOOBPermutedVarDeltaMeanMargin, ...
                        slicedPrivOOBPermutedVarCountRaiseMargin ] = ...
                        oobPermVarUpdate(xorig,y,w,classNames,...
                        localcompact,1,oobtf,doclassregtree,s);
                else
                    slicedPrivOOBPermutedVarDeltaError = ...
                        oobPermVarUpdate(xorig,y,w,classNames,...
                        localcompact,1,oobtf,doclassregtree,s);
                end
            end
        end
        
        if method(1)=='r'
            slicedInBagIndices = idxtrain;
        end
        
        % Report progress
        if nprint>0
            ntreesGrown = internal.stats.parallel.statParallelStore('ntreesGrown') + 1;
            internal.stats.parallel.statParallelStore('ntreesGrown',ntreesGrown);
            if floor(ntreesGrown/nprint)*nprint==ntreesGrown
                if usePool
                    fprintf(1,'%s\n',getString(message('stats:TreeBagger:TreesDoneOnWorker', ...
                        ntreesGrown, ...
                        internal.stats.parallel.statParallelStore('mylabindex'))));
                else
                    fprintf(1,'%s\n',getString(message('stats:TreeBagger:TreesDone',ntreesGrown)));
                end
            end
        end
        
        slicedTree{1} = tree;
        
    end %-nested function loopBody
end


function [slicedPrivOOBPermutedVarDeltaError, ...
    slicedPrivOOBPermutedVarDeltaMeanMargin, ...
    slicedPrivOOBPermutedVarCountRaiseMargin] = ...
    oobPermVarUpdate(x,y,w,classNames,compact,treeInd,oobtf,doclassregtree,s)

% oobPermVarUpdate:
% The output arguments correspond to TreeBagger Properties
% PrivOOBPermutedVarDeltaError, PrivOOBPermutedVarDeltaMeanMargin,
% and PrivOOBPermutedVarCountRaiseMargin, respectively.
% They are supplied as return values because in situ assignments
% to class properties cannot be done in a parfor context.

% Permute values across each input variable
% and estimate decrease in margin due to permutation
Nvars = size(x,2);

% Preallocate the output arguments
slicedPrivOOBPermutedVarDeltaError = zeros(1,Nvars);
slicedPrivOOBPermutedVarDeltaMeanMargin = zeros(1,Nvars);
slicedPrivOOBPermutedVarCountRaiseMargin = zeros(1,Nvars);

% Get oob data
Xoob = x(oobtf,:);

% Get size of oob data
Noob = size(Xoob,1);
if Noob<=1
    return;
end

% Get weights
Woob = w(oobtf);
Wtot = sum(Woob);
if Wtot<=0
    return;
end

% Get non-permuted scores and labels
[sfit,~,yfit] = treeEval(compact,treeInd,Xoob,doclassregtree);

% Get error
doclass = ~isempty(classNames);
    if iscategorical(yfit)
        yfit = cellstr(yfit);
    end
if doclass
    err = dot(Woob,~strcmp(y(oobtf),yfit))/Wtot;
else
    err = dot(Woob,(y(oobtf)-yfit).^2)/Wtot;
end

% Get margins; for classification only
if doclass
    C = classreg.learning.internal.classCount(classNames,y(oobtf));
    margin = classreg.learning.loss.classifmargin(C,sfit);
end

if doclassregtree
    % classregtree has been removed. Use fitctree or fitrtree
    % instead.
    error(message('stats:obsolete:ReplaceThisWithObjectUsingOptions','classregtree','TreeBagger trained on','fitctree','fitrtree'));
else
    usedVars = predictorImportance(compact.Trees{treeInd}) > 0;
end

for ivar=1:Nvars
    % Is this a useful variable?
    if ~usedVars(ivar) 
        continue;
    end
    
    % Get permuted scores and labels
    permuted = randsample(s,Noob,Noob);
    xperm = Xoob;
    xperm(:,ivar) = xperm(permuted,ivar);
    wperm = Woob(permuted);
    [sfitvar,~,yfitvar] = ...
        treeEval(compact,treeInd,xperm,doclassregtree);
    
    % Get the change in error
    if iscategorical(yfitvar)
        yfitvar = cellstr(yfitvar);
    end
    if doclass
        permErr = dot(wperm,~strcmp(y(oobtf),yfitvar))/Wtot;
    else
        permErr = dot(wperm,(y(oobtf)-yfitvar).^2)/Wtot;
    end
    slicedPrivOOBPermutedVarDeltaError(ivar) = permErr-err;
    
    % Get shifts in margins; for classification only
    if doclass
        permMargin = classreg.learning.loss.classifmargin(C,sfitvar);
        deltaMargin = margin-permMargin;
        slicedPrivOOBPermutedVarDeltaMeanMargin(ivar) = ...
            dot(wperm,deltaMargin)/Wtot;
        slicedPrivOOBPermutedVarCountRaiseMargin(ivar) = ...
            sum(deltaMargin>0) - sum(deltaMargin<0);
    end
end
end

function newx = decodeX(oldx,vr)
% Decode a column of X using the variable range and group numbers
if isempty(vr)
    newx = oldx;
else
    ok = oldx>0 & ~isnan(oldx);
    if all(ok) 
        newx = vr(oldx);
    else
        newx = vr(ones(length(ok),1));
        newx(ok,:) = vr(oldx(ok));
        if iscategorical(newx)
            missing = '<undefined>';
        elseif isfloat(newx)
            missing = NaN;
        elseif iscell(newx) % cellstr is the only cell-valued predictor allowed
            missing = {''};
        elseif ischar(newx)
            missing = ' ';
        end
        newx(~ok,:) = missing;
    end
end
end
