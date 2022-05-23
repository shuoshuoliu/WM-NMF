%PARALLELSTATS Parallel computation in the Statistics and Machine Learning Toolbox
%  Several commands in the Statistics and Machine Learning Toolbox can make
%  use of multiple processors to run computations in parallel. If you
%  request parallel computation, a Statistics and Machine Learning Toolbox
%  function will perform iterative calculations in parallel if you have a
%  license for the Parallel Computing Toolbox (PCT) and you have a MATLAB
%  pool open. Parallel computation incurs extra overhead, compared to
%  serial computation, but it can speed up costly iterative calculations.
%
%  You can use three optional parameters to control Statistics and Machine
%  Learning Toolbox functions that support parallel computation. You supply
%  the parameters as fields of a struct, which you may construct using the
%  STATSET function. If you are new to parallel computation in MATLAB, you
%  may wish to consult the Glossary of Terms at the end of this text, prior
%  to reading the parameter descriptions. The parameter names and a
%  description of their values are:
%
%  'UseParallel'     Indicates whether eligible functions should use
%                    capabilities of the Parallel Computing Toolbox (PCT),
%                    if the capabilities are available -- that is, if the
%                    PCT is installed and a license is available. Valid
%                    values are FALSE (the default), to indicate that you
%                    want serial computation, and TRUE, to choose
%                    parallel computation.
%
%                    It is valid to set 'UseParallel' to TRUE, even if the
%                    PCT is not available.  In this case, the Statistics
%                    and Machine Learning Toolbox command will run in
%                    serial mode on the client. This behavior allows
%                    "parallelized" code to be portable, and it facilitates
%                    more convenient debugging.  However, when debugging
%                    you should take into account that iterations may be
%                    performed in different order when in parallel mode
%                    than they would with 'UseParallel' equal to FALSE.
%
%  'UseSubstreams'   Many parallelized toolbox commands generate random
%                    numbers as part of their calculations.  This parameter
%                    indicates whether the command should use a separate
%                    substream of the random number generator for each
%                    iteration in loops that are parallelized.  By doing
%                    so, the command generates results in such a way that
%                    you can reproduce them, either in serial or parallel.
%                    See below for further discussion and examples. Valid
%                    values are TRUE, or do use a separate substream
%                    for each iteration, or FALSE, do not use a separate
%                    substream for each iteration. The value TRUE can
%                    only be used with random stream types that support
%                    substreams. Defaults to FALSE.
%
%  'Streams'         A random number stream (i.e., an object of the
%                    RandStream class), or a cell array of streams. If no
%                    value is supplied for this parameter, the default
%                    random stream on each MATLAB session will be used to
%                    generate random numbers. Otherwise 'Streams' supplies
%                    the stream or streams to be used by the MATLAB
%                    sessions involved in the computation.
%
%                    The rules governing the 'Streams' parameter depend on
%                    the 'UseSubstreams' parameter. If 'UseSubstreams' has
%                    the value TRUE, then 'Streams' must be a single
%                    random number stream, or a scalar cell array
%                    containing a single stream.
%
%                    If 'UseSubstreams' is FALSE, and the computation is
%                    in serial mode, then 'Streams' must be a single random
%                    number stream, or a scalar cell array containing a
%                    single stream.
%
%                    If 'UseSubstreams' is FALSE, and the computation is in
%                    parallel, then the number of random streams supplied
%                    by 'Streams' varies from function to function in the
%                    Statistics and Machine Learning Toolbox. There are two
%                    basic patterns. Some functions only generate random
%                    numbers on the client, even though multiple parallel
%                    workers are used elsewhere during the calculations.
%                    For these commands, 'Streams' must supply a single
%                    stream. Other functions generate random numbers on
%                    every parallel worker involved in the parallel
%                    computations.  For these functions, 'Streams' must
%                    have the same length as the size of the parallel pool.
%                    Consult each parallelized Statistics and Machine
%                    Learning Toolbox function for its particular
%                    requirements.
%
%                    When 'UseParallel' is true, even if a parallel pool is
%                    not open, then depending on your installation and
%                    preferences, the PCT may open a parallel pool for you.
%                    For any Statistics and Machine Learning Toolbox
%                    function that requires 'Streams' to have the same
%                    length as the size of the parallel pool, it is best to
%                    open a pool explicitly, using the GCP or PARPOOL
%                    command, before calling the function with
%                    'UseParallel' equal to TRUE.  See Example 3 below,
%                    Case (C), for further explanation.
%
%  RANDOM STREAMS IN PARALLEL
%
%  When you open a pool of MATLAB workers, each worker creates its own
%  default random number stream. The default streams are initialized in
%  such a way that you can treat all the different streams as being
%  statistically independent. Thus, the result of a computation is
%  statistically equivalent to the same computation performed in serial
%  mode.  Serial and parallel computations may generate different random
%  values, but both will produce equally valid results. You do not need to
%  make any special initialization of the parallel workers' default random
%  number streams in order to use one of the Statistics and Machine
%  Learning Toolbox functions that supports parallel computation.
%
%  Instead of using the parallel workers' default random number streams,
%  you can use the 'Streams' parameter to supply your own streams to
%  Statistics and Machine Learning Toolbox functions that support parallel
%  computation. The streams you supply are used within that function call,
%  and they do not replace the workers' default streams.
%
%  Note that whenever you modify the workers' default random number
%  streams, or replace them for a calculation, you must take care that the
%  changes you make do not invalidate the statistical independence of the
%  streams in the MATLAB pool. A safe and convenient way to create random
%  streams that work correctly in a MATLAB pool is to use one of the MATLAB
%  random stream types that support multiple streams. See the section on
%  "Multiple streams" in the MATLAB User Guide.
%
%  SUBSTREAMS AND REPRODUCIBILITY
%  
%  Most often, you want computations that use random numbers to exhibit
%  "ongoing" randomness. That is, you expect the results to continue to
%  vary randomly if you repeat the calculation. However, in some cases it
%  is important to be able to reproduce the same result over and over, by
%  recreating a previous set of initial conditions, and then generating
%  exactly the same results during a repeat of the computation.  For
%  Statistics and Machine Learning Toolbox functions that support parallel
%  computation, this kind of reproducibility requires that you use the
%  'UseSubstreams' parameter. When operating in parallel mode, iterative
%  calculations may be assigned to workers in an unpredictable fashion, and
%  iterations may occur in an unpredictable order. Therefore, simply
%  putting the workers' random streams in a known state before calling a
%  Statistics and Machine Learning Toolbox function does not guarantee
%  predictable results. If 'UseSubstreams' is TRUE, the Statistics and
%  Machine Learning Toolbox function associates each loop iteration with a
%  precalculated substream index, or starting point in the random stream
%  being used. Therefore it does not matter which worker performs the
%  iteration, or whether the computation is in serial or parallel mode.
%
%  FOR FURTHER INFORMATION
%
%  The MATLAB User Guide provides general information about creating and
%  managing random number streams, and the Statistics and Machine Learning
%  Toolbox User Guide has additional information about managing random
%  streams in a MATLAB pool.
%
%  GLOSSARY OF TERMS
%
%  MATLAB Session         A running instance of the MATLAB software; a
%                         MATLAB executable.
%  parallel pool          A collection of MATLAB sessions running on
%                         separate processors, with intercommunication
%                         established between processors so that they can
%                         collaborate jointly on a task.
%  serial computation     Computation performed on a single MATLAB session.
%                         This is also called "serial mode".
%  parallel computation   Computation that takes place using all the MATLAB
%                         Sessions in a parallel pool. Also called "parallel
%                         mode".
%  client                 The primary MATLAB session, which gives you the
%                         command window. You always interact with the
%                         client. When there is a parallel pool open, the
%                         client interacts on your behalf with the workers
%                         in the parallel pool. The client is also the MATLAB
%                         session that performs serial mode computation.
%  worker                 A MATLAB session in a parallel pool. The worker
%                         does not have a graphical display and does not
%                         interact with you directly.
%
%  EXAMPLES:
%
%  EXAMPLE 1
%
%  The following example demonstrates parallel computation with default
%  random number streams. This is the simplest and probably the most common
%  type of use of parallelized Statistics and Machine Learning Toolbox
%  functions. The "parallel" form of TreeBagger works whether you are
%  licensed for PCT or not. If you have a parallel pool open (this requires
%  the PCT), TreeBagger will build fifty decision trees using multiple
%  processors in parallel. If you do not have a parallel pool open, MATLAB
%  may try to open a pool for you (depending on your preferences and
%  whether you have the PCT). If MATLAB succeeds in opening a parallel
%  pool, TreeBagger will again build the trees using parallel computation.
%  Otherwise, the computations will run in serial mode.
%  
%    load fisheriris
%    opt = statset('UseParallel',true);
%    b = TreeBagger(50, meas, species, 'oobpred','on', 'opt',opt);
%    plot(oobError(b))
%    xlabel('number of grown trees')
%    ylabel('out-of-bag classification error')
%
%  EXAMPLE 2
%
%  This example demonstrates how to use the 'UseSubstreams' parameter to
%  generate results that you can replicate on a subsequent function call,
%  either in serial or parallel mode. Note that you can run this example
%  only if you have a license for PCT, otherwise the "gcp" command is
%  invalid.
%
%    % We will start with an open MATLAB pool.
%    try
%        p = gcp();
%        if isempty(p)
%            nWorkers = 0;
%        else
%            nWorkers = p.NumWorkers;
%        end
%    catch ME
%        nWorkers = 0;
%    end
% 
%    % If nWorkers == 0, MATLAB could not open a parallel pool.
%    % You need to determine why before proceeding with the example.
%
%    % You must use a random number stream that supports substreams
%    stream = RandStream('mrg32k3a');
%    optpar = statset('UseParallel',true, ...
%                     'UseSubstreams',true, ...
%                     'Streams',stream);
%
%    % Run bootstrp() in parallel mode with the 'UseSubstreams' parameter
%    y = exprnd(5,100,1);
%    bootstatPar = bootstrp(100, @mean, y, 'options', optpar);
%
%    % Set the random number stream to the state it had before you ran the
%    % parallel bootstrp.
%    reset(stream)
%
%    % Now run the same bootstrp command, again in parallel
%    bootstatPar2 = bootstrp(100, @mean, y, 'options', optpar);
%
%    % The results should be identical
%    isequal(bootstatPar, bootstatPar2)
%
%    % Reset the stream again 
%    reset(stream)
%
%    % Rerun the bootstrp() command, this time in serial mode. To do so,
%    % you can create a serial-mode Options parameter, as below. This will
%    % tell the bootstrp command to run in serial mode on the client, even
%    % though a parallel pool remains open (you also get serial computation
%    % if you simply do not specify 'UseParallel', since 'false' is the
%    % default value).
%
%    optser = statset('UseParallel',false, ...
%                     'UseSubstreams',true, ...
%                     'Streams',stream);
%    bootstatSer = bootstrp(100, @mean, y, 'options', optser);
%
%    % Again the results should be identical
%    isequal(bootstatPar, bootstatSer)
%
%  EXAMPLE 3
%
%  This example demonstrates use of the 'Streams' parameter when your goal is
%  basic randomness, not reproducibility. The example can only be run if you
%  have the PCT installed and licensed. The example will cover three cases:
%  (A) A MATLAB pool is open, and the Statistics and Machine Learning
%      Toolbox function uses random numbers on the workers.
%  (B) A MATLAB pool is open, and the Statistics and Machine Learning
%      Toolbox function uses random numbers on the client (though it uses
%      the workers for other parts of the computation).
%  (C) The Statistics and Machine Learning Toolbox function would use
%      random numbers on the workers if there were a MATLAB pool open, but a
%      MATLAB pool is not currently open.
%
%    % CASE (A)
%    %
%    % TreeBagger generates random numbers on all of the workers. 
%    % If you are supplying your own streams, you must supply one for each 
%    % worker in the MATLAB pool.
%
%    % We will start with an open MATLAB pool
%    try
%        p = gcp();
%        if isempty(p)
%            nWorkers = 0;
%        else
%            nWorkers = p.NumWorkers;
%        end
%    catch ME
%        nWorkers = 0;
%    end
% 
%    % If nWorkers == 0, MATLAB could not open a parallel pool.
%    % You need to determine why before proceeding with the example.
%
%    poolsz = nWorkers;
%    streams = RandStream.create('mlfg6331_64', ...
%                                'NumStreams', poolsz, ...
%                                'StreamIndices', 1:poolsz, ...
%                                'CellOutput', true);
%    optMany = statset('UseParallel',true,'Streams',streams);
%    load fisheriris
%    b = TreeBagger(50, meas, species, 'oobpred','on', 'options',optMany);
%    plot(oobError(b))
%    xlabel('number of grown trees')
%    ylabel('out-of-bag classification error')
%
%    % Case (B) 
%    %
%    % The crossval function only uses random numbers on the client. 
%    % Therefore it only recognizes a single stream, even if there is 
%    % an open parallel pool.
%
%    % We will start with an open MATLAB pool
%    try
%        p = gcp();
%        if isempty(p)
%            nWorkers = 0;
%        else
%            nWorkers = p.NumWorkers;
%        end
%    catch ME
%        nWorkers = 0;
%    end
% 
%    % If nWorkers == 0, MATLAB could not open a parallel pool.
%    % You need to determine why before proceeding with the example.
%
%    soloStream = RandStream('mlfg6331_64');
%    optSolo = statset('UseParallel',true, 'Streams',soloStream);
%    load('fisheriris');
%    y = meas(:,1);
%    X = [ones(size(y,1),1),meas(:,2:4)];
%    regf=@(XTRAIN,ytrain,XTEST)(XTEST*regress(ytrain,XTRAIN));
%    cvMse = crossval('mse', X, y, 'predfun',regf, 'options',optSolo); 
%
%    % Case (C)
%    %
%    % Repeat case (A) but without first opening a parallel a pool. 
%
%    % Make sure there is no pool open.
%    p = gcp('nocreate');
%    delete(p)
%
%    % Depending on your preferences and whether you have PCT available 
%    % when 'UseParallel' is true, MATLAB may try to open a parallel pool.
%    % If it succeeds, then the value of 'Streams' must match the number 
%    % of workers in the pool that was created.  Say that you know that 
%    % your default parallel pool size is 4.  Then you can issue:
%    poolsz = 4;
%    streams = RandStream.create('mlfg6331_64', ...
%                                'NumStreams', poolsz, ...
%                                'StreamIndices', 1:poolsz, ...
%                                'CellOutput', true);
%    optMany = statset('UseParallel',true,'Streams',streams);
%    load fisheriris
%    b = TreeBagger(50, meas, species, 'oobpred','on', 'options',optMany);
%
%    % The above statement will succeed if a pool was opened with number
%    % of workers equal to your expected size.  Otherwise it will fail.
%
%    % If a pool failed to open, then there is only one MATLAB session 
%    % available (the client).  Thus TreeBagger expects only a single
%    % stream, regardless of the value of 'UseParallel'.  In that case,
%    % the following command sequence would work:
%
%    soloStream = RandStream('mlfg6331_64');
%    optSolo = statset('UseParallel',true, 'Streams',soloStream);
%    b2 = TreeBagger(50, meas, species, 'oobpred','on', 'options',optSolo);
%
%    % Because it is difficult to ensure the outcome of the auto-opening
%    % of the parallel pool, it is best to use the GCP or PARPOOL command
%    % and explicitly create a parallel pool prior to using parallelized
%    % Statistics and Machine Learning Toolbox commands with 'UseParallel'
%    % set to TRUE.
%    
%  See also statset, statget, parfor, parpool, gcp, RandStream, RandStream/Substream.
%
%  For more information on parallel computing concepts and how to work 
%  with parallelized Statistics and Machine Learning Toolbox functions, see the chapter on 
%  <a href = "matlab: helpview([docroot '/toolbox/stats/stats.map'], 'ch_parallel_stats');">Parallel Statistics</a> in the User Guide.

%  Copyright 2009-2014 The MathWorks, Inc.
%

help parallelstats
