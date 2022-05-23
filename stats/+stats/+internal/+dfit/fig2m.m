function fig2m(dffig,outfilename)
%FIG2M Turn figure into generated code that can produce the figure


%   Copyright 2003-2020 The MathWorks, Inc.

dsdb = stats.internal.dfit.getdsdb;
fitdb = stats.internal.dfit.getfitdb;
if isempty(down(dsdb)) && isempty(down(fitdb))
   errordlg(getString(message('stats:dffig2m:DBRequired')),...
       getString(message('stats:dffig2m:DlgNameDBRequired')),...
       'modal');
   return
end

if nargin<1
   dffig = dfgetset('dffig'); 
end
if nargin<2
   outfilename = '';
end

% Get file name with .m suffix, and get corresponding function name
if length(outfilename)<2 || ~isequal(outfilename(end-1:end),'.m')
   outfilename = sprintf('%s.m',outfilename);
end
fcnname = outfilename(1:end-2);
k = find(fcnname(1:end-1)=='\',1,'last');
if ~isempty(k)
   fcnname = fcnname(k+1:end);
end
k = find(fcnname(1:end-1)=='/',1,'last');
if ~isempty(k)
   fcnname = fcnname(k+1:end);
end
if isempty(fcnname)
    fcnname = 'createFit';
end
   
% Set up some variables for later
allprop = {'Color' 'Marker' 'LineStyle' 'LineWidth' 'MarkerSize'};
ftype = dfgetset('ftype');
alpha = 1 - dfgetset('conflev');
showlegend = isequal(dfgetset('showlegend'),'on');
axold = findobj(dffig, 'Tag', 'dfMainAxes');
leginfo = {};
if showlegend
   legh = get(axold,'Legend');
   leginfo = stats.internal.dfit.getlegendinfo(legh);
   legloc = get(legh,'Location');
   if isequal(legloc,'none')
      oldu = get(legh,'units');
      legpos = get(legh,'position');
      legpos = hgconvertunits(dffig,legpos,oldu,'normalized',dffig);
      legloc = legpos(1:2);
   end
end

% Create arrays to receive code text
blkc = cell(0,1);    % block of comment lines
blks = cell(0,1);    % block of setup lines
blkd = cell(0,1);    % block of data-related lines
blkf = cell(0,1);    % block of fit-related lines
blke = cell(0,1);    % block of lines at end

% Examine datasets to extract the variable names that they use
[blkc,dsvarnames,exprlist,arglist] = getarglist(blkc,dsdb);

% Define some variable names, making sure they don't conflict with argument
% names
leghname = genvarname('LegHandles',arglist);
legtname = genvarname('LegText',arglist);
xlimname = genvarname('XLim',arglist);
hlegendname = genvarname('hLegend',arglist);
hlinename = genvarname('hLine',arglist);
hboundsname = genvarname('hBounds',arglist);
xgridname = genvarname('XGrid',arglist);
xdiscretename = genvarname('XDiscrete',arglist);
xcontname = genvarname('XContinuous',arglist);
incrname = genvarname('XIncr',arglist);
probdistspecname = genvarname('CustomDist',arglist);
positionname = genvarname('Position',arglist);

% We also want a series of names pd1, pd2, etc., so make sure there is
% nothing like that in the list now
pdprefix = getpdprefix(arglist);

% Write set-up code
blks{end+1} = '% Prepare figure';
blks{end+1} = 'clf;';
blks{end+1} = 'hold on;';
if showlegend
   blks{end+1} = sprintf('%s = []; %s = {};   % handles and text for legend',...
                         leghname, legtname);
end

% Preliminary pass, see how many datasets are to be plotted
ds = down(dsdb);
ndsplotted = 0;
dtype = dfgetset('dtype');
while(~isempty(ds))
   dsline = ds.line;
   if ~isempty(dsline) && ishghandle(dsline)
      ndsplotted = ndsplotted + 1;
   end
   ds = right(ds);
end
simplecase = ischar(dtype) && (ndsplotted == 1);

% Process each dataset for plotting
ds = down(dsdb);
numds = 0;
while(~isempty(ds))
   numds = numds + 1;
   [blkd,arglist,showbounds,onplot] = writedset(blkd,ds,arglist,...
       dsvarnames(numds,:),allprop,alpha,hlinename,hboundsname,simplecase);

   if onplot && showlegend
      blkd{end+1} = sprintf('%s(end+1) = %s;',leghname,hlinename); %#ok<*AGROW>
      blkd{end+1} = sprintf('%s{end+1} = ''%s'';',legtname,quotedtext(ds.name));
      if showbounds
         blkd{end+1} = sprintf('%s(end+1) = %s;',leghname,hboundsname);
         blkd{end+1} = sprintf('%s{end+1} = ''%g%% confidence bounds'';',...
                               legtname,100*(1-alpha));
      end
   end
   ds = right(ds);
end

% Force all inputs to be column vectors
blkc{end+1} = ' ';
blkc{end+1} = '% Force all inputs to be column vectors';
for j=1:length(arglist)
   blkc{end+1} = sprintf('%s = %s(:);',arglist{j},arglist{j});
end

% Set up for plotting fits
anycontinuous = false;
anydiscrete = false;
ft = down(fitdb);
while(~isempty(ft))
   if ft.iscontinuous
      anycontinuous = true;
   else
      anydiscrete = true;
   end
   ft = right(ft);
end

% Process each fit
numfit = 0;
ft = down(fitdb);
anySmoothFits = false;
while(~isempty(ft))
   numfit = numfit+1;
   fitname = ft.name;
   outname = sprintf('%s%d',pdprefix,numfit);

   % Create code to re-create this fit
   blkf{end+1} = sprintf('\n%% --- Create fit "%s"',fitname);

   % Call subfunction to generate code for each type
   if isequal(getfittype(ft),'param')
      [blkf,showbounds,onplot] = writepfit(blkf,outname,ft,alpha,allprop,...
         anycontinuous,anydiscrete,exprlist,arglist,...
         hlinename,hboundsname,xgridname,xcontname,xdiscretename,probdistspecname);
   else
      anySmoothFits = true;
      [blkf,onplot] = writenpfit(blkf,outname,ft,allprop,...
                                 anycontinuous,anydiscrete,exprlist,arglist,...
                                 hlinename,xgridname,xcontname);
      showbounds = false;
   end

   % Add legend if requested
   if onplot && showlegend
      blkf{end+1} = sprintf('%s(end+1) = %s;',leghname,hlinename);
      blkf{end+1} = sprintf('%s{end+1} = ''%s'';',legtname,quotedtext(ft.name));
      if showbounds
         blkf{end+1} = sprintf('%s(end+1) = %s;',leghname,hboundsname);
         blkf{end+1} = sprintf('%s{end+1} = ''%g%% confidence bounds'';',...
                               legtname,100*(1-alpha));
      end
   end
   ft = right(ft);
end

% In setup section, create empty axes for prob plot except in simple case
if isequal(ftype,'probplot') && ~simplecase
   if ischar(dtype)
      blks{end+1} = sprintf('probplot(''%s''); %% create empty plot of desired type', dtype);
   else
      blks{end+1} = sprintf(...
          '%s = internal.stats.getdistributions(''%s'');',...
          probdistspecname,dtype.distspec.code);
      blks{end+1} = sprintf('probplot({%s,%s});',...
                            probdistspecname,cell2text(num2cell(dtype.params)));
   end
   blks{end+1} = sprintf('title('''');');
end


% At end of data set section, set x axis limits
if ~isequal(ftype,'probplot')
    if ndsplotted==0
        % In this unusual case we have to select a range for the fit
        if isequal(ftype,'icdf')
            blkd{end+1} = sprintf('%s = [0.001, 0.999];',xlimname);
        else
            blkd{end+1} = sprintf('\n%% Get data limits to determine plotting range');
            ds = down(dsdb);
            yname = expression2name(ds.yname,exprlist,arglist);
            blkd{end+1} = sprintf('%s = [min(%s), max(%s)];',xlimname,yname,yname);
            
            
            blkd{end+1} = sprintf('%s = [min(%s), max(%s)];',xlimname,yname,yname);
            while(true)
                ds = right(ds);
                if isempty(ds)
                    break
                end
                yname = expression2name(ds.yname,exprlist,arglist);
                blkd{end+1} = sprintf('%s(1) = min(%s(1), min(%s));',xlimname,xlimname,yname);
                blkd{end+1} = sprintf('%s(2) = max(%s(2), max(%s));',xlimname,xlimname,yname);
            end
        end
    end
    
    % Most functions are plotted over a grid
    blkd{end+1} = sprintf('\n%% Create grid where function will be computed');
    if ndsplotted>0
        % This is the typical case where we plot the data and evaluate the fit
        % over the plotted range
        blkd{end+1} = sprintf('%s = get(gca,''XLim'');',xlimname);
    end
    if ~isequal(ftype,'icdf')
        if isequal(get(axold,'XScale'),'log')
            blkd{end+1} = sprintf('   %s = exp(log(%s) + [-1 1] * 0.01 * diff(log(%s)));',...
                xlimname,xlimname,xlimname);
        else
            blkd{end+1} = sprintf('   %s = %s + [-1 1] * 0.01 * diff(%s);',...
                xlimname,xlimname,xlimname);
        end
    end
    
    % Create a suitable X vector, may depend on whether it's discrete
    if ~isequal(ftype,'pdf') || ~anydiscrete
        blkd{end+1} = sprintf('%s = linspace(%s(1),%s(2),100);',xgridname,xlimname,xlimname);
    elseif ~anycontinuous
        blkd{end+1} = sprintf('%s = max(1,floor((%s(2)-%s(1))/100));',incrname,xlimname,xlimname);
        blkd{end+1} = sprintf('%s = floor(%s(1)):%s:ceil(%s(2));',xgridname,xlimname,incrname,xlimname);
    else
        blkd{end+1} = sprintf('%s = linspace(%s(1),%s(2),100);',xcontname,xlimname,xlimname);
        blkd{end+1} = sprintf('%s = max(1,floor((%s(2)-%s(1))/100));',incrname,xlimname,xlimname);
        blkd{end+1} = sprintf('%s = floor(%s(1)):%s:ceil(%s(2));',xdiscretename,xlimname,incrname,xlimname);
    end
end

% Finish up
blke{end+1} = '% Adjust figure';
blke{end+1} = sprintf('box on;');
if isequal(dfgetset('showgrid'),'on')
   blke{end+1} = sprintf('grid on;');
end
blke{end+1} = 'hold off;';

if showlegend
    if isempty(leginfo)
        legargstext = '';
    else
        legargstext = [',' cell2text(leginfo,'list')];
    end

   blke{end+1} = '';
   blke{end+1} = '% Create legend from accumulated handles and labels';
   if isnumeric(legloc)
      blke{end+1} = sprintf('%s = legend(%s,%s%s); % create and reposition legend',...
          hlegendname,leghname,legtname,legargstext);
      blke{end+1} = sprintf('set(%s,''Units'',''normalized'');',hlegendname);
      blke{end+1} = sprintf('%s = get(%s,''Position'');',positionname,hlegendname);
      blke{end+1} = sprintf('%s(1:2) = [%g,%g];',positionname,legloc);
      blke{end+1} = sprintf('set(%s,''Interpreter'',''none'',''Position'',%s);',hlegendname,positionname);
   else
      blke{end+1} = sprintf('%s = legend(%s,%s%s);  % create legend',...
          hlegendname,leghname,legtname,legargstext);
      blke{end+1} = sprintf('set(%s,''Interpreter'',''none'');',hlegendname);
   end
end

% Get text for the function output, if any
if numfit==0
    outtext = '';
    outcomment = '';
elseif numfit==1
    outlist = [pdprefix '1'];
    outtext = [outlist ' = '];
    outcomment = sprintf('%% Output fitted probablility distribution: %s', upper(outlist));
else
    outlist = [pdprefix '1'];
    for j=2:numfit
        outlist = sprintf('%s,%s%d',outlist,pdprefix,j);
    end
    outtext = ['[' outlist '] = '];
    outcomment = sprintf('%% Output fitted probablility distributions: %s', upper(outlist));
end
if ~isempty(outcomment)
    blkc = [{''} {outcomment} blkc];
end

% Get text for the function input, if any
if isempty(arglist)
   argtext = '';
else
   argtext = sprintf('%s,',arglist{:});
   argtext = sprintf('(%s)',argtext(1:end-1));
end

codestring = writecode(blkc,blks,blkd,blkf,blke,outtext,fcnname,...
                       argtext,numds,numfit,anySmoothFits,ftype);

% Write code into file (if given) or editor
if length(outfilename)<3    % no file name, just .m suffix
    editorDoc = matlab.desktop.editor.newDocument(codestring);
    editorDoc.smartIndentContents();
    editorDoc.goToPositionInLine(1,1);
else
    [fid,msg] = fopen(outfilename,'w');
    if fid==-1
        emsg = sprintf('Error trying to write to %s:\n%s',outfilename,msg);
        errordlg(emsg,'Error Saving File','modal');
        return
    end
    
    fprintf(fid,'%s',codestring);
    fclose(fid);
end

function codestring = writecode(blkc,blks,blkd,blkf,blke,outtext,fcnname,...
                                argtext,numds,numfit,anySmoothFits,ftype)

codestring = '';
codestring = sprintf('%sfunction %s%s%s\n',codestring,outtext,fcnname,argtext);
codestring = sprintf('%s%%%s    Create plot of datasets and fits\n',codestring,upper(fcnname));
codestring = sprintf('%s%%   %s%s%s\n',codestring,upper(outtext),upper(fcnname),upper(argtext));
codestring = sprintf('%s%%   Creates a plot, similar to the plot in the main distribution fitter\n',codestring);
codestring = sprintf('%s%%   window, using the data that you provide as input.  You can\n',codestring);
codestring = sprintf('%s%%   apply this function to the same data you used with distributionFitter\n',codestring);
codestring = sprintf('%s%%   or with different data.  You may want to edit the function to\n',codestring);
codestring = sprintf('%s%%   customize the code and this help message.\n',codestring);
codestring = sprintf('%s%%\n',codestring);
codestring = sprintf('%s%%   Number of datasets:  %d\n',codestring,numds);
codestring = sprintf('%s%%   Number of fits:  %d\n',codestring,numfit);
codestring = sprintf('%s%%\n',codestring);
codestring = sprintf('%s%%   See also FITDIST.\n',codestring);
codestring = sprintf('%s\n',codestring);
codestring = sprintf('%s%% This function was automatically generated on %s\n',...
            codestring,datestr(now));
for j=1:length(blkc)
   codestring = sprintf('%s%s\n',codestring,blkc{j});
end
codestring = sprintf('%s\n',codestring);
for j=1:length(blks)
   codestring = sprintf('%s%s\n',codestring,blks{j});
end
codestring = sprintf('%s\n',codestring);
for j=1:length(blkd)
   codestring = sprintf('%s%s\n',codestring,blkd{j});
end
codestring = sprintf('%s\n',codestring);
for j=1:length(blkf)
   codestring = sprintf('%s%s\n',codestring,blkf{j});
end
codestring = sprintf('%s\n',codestring);
for j=1:length(blke)
   codestring = sprintf('%s%s\n',codestring,blke{j});
end

% Create sub function to be used to support a functionline fit on a probability plot
if anySmoothFits && isequal(ftype,'probplot')
   codestring = sprintf('%s\n\n%% -----------------------------------------------\n',codestring);
   codestring = sprintf('%sfunction f=cdfnp(x,y,cens,freq,support,kernel,width)\n',codestring);
   codestring = sprintf('%s%%CDFNP Compute cdf for non-parametric fit, used in probability plot\n\n',codestring);
   codestring = sprintf('%sf = ksdensity(y,x,''cens'',cens,''weight'',freq,''function'',''cdf'',...\n',codestring);
   codestring = sprintf('%s                  ''support'',support,''kernel'',kernel,''width'',width);\n',codestring);
end


% ------------------- double up quotes in text string
function a = quotedtext(b)
if ischar(b)
   a = strrep(b,'''','''''');
else
   a = sprintf('%.13g',b);
end

% ------------------- create text to re-create cell or numeric array
function a = cell2text(b,target)

% This is not a completely general-purpose routine, but it handles the sort
% of cell arrays used here.  A cell array containing a matrix, for
% instance, would not work
if ~iscell(b)
   if ischar(b)
      a = sprintf('''%s''',quotedtext(b));
   elseif length(b)==1
      a = sprintf('%.13g',b);
   else
      numtext = num2str(b,'%.13g ');
      if size(numtext,1)>1
         numtext = [numtext repmat(';',size(numtext,1),1)]';
         numtext = numtext(:)';
         numtext = numtext(1:end-1);
      end
      a = sprintf('[%s]',numtext);
   end
   return
end

if ~isempty(b)
   bj = b{1};
   if ischar(bj)
      a = sprintf('''%s''',quotedtext(bj));
   else
      a = sprintf(' %.13g',bj);
   end
   for j=2:length(b)
      bj = b{j};
      if ischar(bj)
         a = sprintf('%s, ''%s''',a,quotedtext(bj));
      elseif isscalar(bj)
         a = sprintf('%s, %.13g',a,bj);
      else
         a = sprintf('%s, [%s]',a,sprintf(' %.13g',bj));
      end
   end
else
   a = '';
end
if nargin<2 || ~isequal(target,'list')
    % Convert to an array unless asked to preserve the list form
    a = sprintf('[%s]',a);
end


% ----------------- add censoring and frequency args to code block
function oneline = addcensfreq(oneline,censname,freqname)

if ~isempty(censname) && ~isequal(censname,'[]')
   oneline = sprintf('%s,''cens'',%s',oneline,censname);
end
if ~isempty(freqname) && ~isequal(freqname,'[]')
   oneline = sprintf('%s,''freq'',%s',oneline,freqname);
end


% ---------------- write code for parametric fit
function [blkf,showbounds,onplot] = ...
    writepfit(blkf,outname,ft,alpha,allprop,anycontinuous,anydiscrete,...
              exprlist,arglist,hlinename,hboundsname,xgridname,xcontname,xdiscretename,probdistspecname)

ds = ft.ds;
yname = expression2name(ds.yname,exprlist,arglist);
dist = ft.distspec;
ftype = ft.ftype;
onplot = true;

blkf{end+1} = sprintf('\n%% Fit this distribution to get parameter values');
[censname,freqname] = getcensfreqname(ds,exprlist,arglist);

% Exclude data if necessary
if ~isempty(ft.exclusionrule)
   [blkf,yname,censname,freqname] = applyexclusion(blkf,ft.exclusionrule,...
                                         yname,censname,freqname,arglist);
end

% Prepare data for fitting
[censname,freqname] = fixcensfreq(censname,freqname);

% Helpful note about using old results instead of fitting new data
if isequal(getfittype(ft),'param')
    blkf{end+1} = sprintf('%% To use parameter estimates from the original fit:');
    blkf{end+1} = sprintf('%%     %s = ProbDistUnivParam(''%s'',%s)', ...
                          outname,dist.code,cell2text(num2cell(ft.params)));
end

if isfield(dist,'fitfunc') && ~isempty(dist.fitfunc)
fname = func2str(dist.fitfunc);
onpath = exist(fname,'file');
else
    onpath = true;
end

prequired = ft.params(dist.prequired == 1); % binomial N, generalized Pareto threshold

% Construct argument list for fitdist
arglist = sprintf('%s, ''%s''', yname,dist.code);
arglist = addcensfreq(arglist,censname,freqname);

% Deal with distributions having fixed parameters
switch(dist.code)
    case 'binomial'
        arglist = sprintf('%s, ''n'', %d', arglist, prequired);
    case 'generalized pareto'
        arglist = sprintf('%s, ''theta'', %d', arglist, prequired);
end

blkf{end+1} = sprintf('%s = fitdist(%s); % fit distribution',outname,arglist);

% Get covariance matrix if we need confidence bounds
if ft.showbounds && ismember(ftype,{'cdf' 'survivor' 'cumhazard' 'icdf'})
   showbounds = true;
else
   showbounds = false;
end

% Sometimes we need the structure that describes the distribution
if ~onpath && (showbounds || isequal(ftype,'probplot'))
   blkf{end+1} = sprintf(...
          '\n%% Get a description of the %s distribution',...
          dist.name);
   blkf{end+1} = sprintf(...
          '%s = internal.stats.getdistributions(''%s'');\n',...
          probdistspecname,dist.code);
end

% Plot the fit and bounds if the original figure had them plotted
if isempty(ft.line) || ~ishghandle(ft.line)
   blkf{end+1} = '% This fit does not appear on the plot';
   onplot = false;
   return;
end

propvals = get(ft.line,allprop);
[c,m,l,w,s] = deal(propvals{:});

yplotname = genvarname('YPlot',arglist);
ylowername = genvarname('YLower',arglist);
yuppername = genvarname('YUpper',arglist);
tmpname = genvarname('YSwap',arglist);

switch(ftype)
 case {'pdf'}
   if anycontinuous && anydiscrete
      if ft.iscontinuous
         blkf{end+1} = sprintf('%s = %s;',xgridname,xcontname);
      else
         blkf{end+1} = sprintf('%s = %s;',xgridname,xdiscretename);
      end            
   end
   blkf{end+1} = sprintf('%s = pdf(%s,%s);',yplotname,outname,xgridname);
   
   blkf{end+1} = sprintf('%s = plot(%s,%s,''Color'',[%g %g %g],...',...
                         hlinename,xgridname,yplotname,c(1),c(2),c(3));
   blkf{end+1} = sprintf('          ''LineStyle'',''%s'', ''LineWidth'',%d,...',l,w);
   blkf{end+1} = sprintf('          ''Marker'',''%s'', ''MarkerSize'',%d);',m,s);
  
 case {'cdf' 'survivor' 'cumhazard' 'icdf'}
   if isequal(ftype,'icdf')
      fname = 'icdf';
   else
      fname = 'cdf';
   end

   if showbounds
       blkf{end+1} = sprintf('[%s,%s,%s] = %s(%s,%s,%g);',...
           yplotname,ylowername,yuppername,fname,outname,xgridname,alpha);
       
   else
       blkf{end+1} = sprintf('%s = %s(%s,%s);',...
           yplotname,fname,outname,xgridname);
   end

   if isequal(ftype,'survivor')
      blkf{end+1} = sprintf('%s = 1 - %s; % convert to survivor function',yplotname,yplotname);
      if showbounds
         blkf{end+1} = sprintf('%s = %s;',tmpname,ylowername);
         blkf{end+1} = sprintf('%s = 1 - %s;',ylowername,yuppername);
         blkf{end+1} = sprintf('%s = 1 - %s;',yuppername,tmpname);
      end
   elseif isequal(ftype,'cumhazard')
      blkf{end+1} = sprintf('%s = -log(1 - %s); % convert to cumulative hazard',yplotname,yplotname);
      if showbounds
         blkf{end+1} = sprintf('if ~isempty(%s)',ylowername);
         blkf{end+1} = sprintf('   %s = %s;',tmpname,ylowername);
         blkf{end+1} = sprintf('   %s = -log(1 - %s);',ylowername,yuppername);
         blkf{end+1} = sprintf('   %s = -log(1 - %s);',yuppername,tmpname);
         blkf{end+1} = 'end';
      end
   end
      
   blkf{end+1} = sprintf('%s = plot(%s,%s,''Color'',[%g %g %g],...',...
                         hlinename,xgridname,yplotname,c(1),c(2),c(3));
   blkf{end+1} = sprintf('          ''LineStyle'',''%s'', ''LineWidth'',%d,...',l,w);
   blkf{end+1} = sprintf('          ''Marker'',''%s'', ''MarkerSize'',%d);',m,s);

   if showbounds
      blkf{end+1} = sprintf('if ~isempty(%s)',ylowername);
      blkf{end+1} = sprintf('   %s = plot([%s(:); NaN; %s(:)], [%s(:); NaN; %s(:)],''Color'',[%g %g %g],...',...
                            hboundsname,xgridname,xgridname,ylowername,yuppername,c(1),c(2),c(3));
      blkf{end+1} = '             ''LineStyle'','':'', ''LineWidth'',1,...';
      blkf{end+1} = '             ''Marker'',''none'');';
      blkf{end+1} = 'end';
   end

 case 'probplot'
   stmt = sprintf('%s = probplot(gca,%s);',hlinename,outname);
   blkf{end+1} = stmt;
   blkf{end+1} = sprintf('set(%s,''Color'',[%g %g %g],''LineStyle'',''%s'', ''LineWidth'',%d);', ...
         hlinename,c(1),c(2),c(3),l,w);
end


% ---------------- write code for nonparametric fit
function [blkf,onplot] = ...
    writenpfit(blkf,outname,ft,allprop,anycontinuous,anydiscrete,exprlist,arglist,hlinename,xgridname,xcontname)

ds = ft.ds;
yname = expression2name(ds.yname,exprlist,arglist);
ftype = ft.ftype;

[censname,freqname] = getcensfreqname(ds,exprlist,arglist);
   
% Exclude data if necessary
if ~isempty(ft.exclusionrule)
   [blkf,yname,censname,freqname] = applyexclusion(blkf,ft.exclusionrule,...
                                          yname,censname,freqname,arglist);
end

% Prepare data for fitting
[censname,freqname] = fixcensfreq(censname,freqname);

kernel = sprintf('''%s''',ft.kernel);
if ft.bandwidthradio == 0
   width = '[]';
else
   width = ft.bandwidthtext;
end
if ischar(ft.support)
   spt = sprintf('''%s''',ft.support);
else
   spt = sprintf('[%g, %g]',ft.support);
end

% Create the fit object
argtxt = addcensfreq('',censname,freqname);
if ~isequal(ft,'unbounded')
    argtxt = sprintf('%s,''support'',%s',argtxt,spt);
end
if ~isequal(width,'[]')
    argtxt = sprintf('%s,''width'',%s',argtxt,width);
end
blkf{end+1} = sprintf('%s = fitdist(%s,''kernel'',''kernel'',%s%s);',...
    outname,yname,kernel,argtxt);

% Plot the fit and bounds if the original figure had them plotted
if isempty(ft.line) || ~ishghandle(ft.line)
   blkf{end+1} = '% This fit does not appear on the plot';
   onplot = false;
   return;
end
onplot = true;

propvals = get(ft.line,allprop);
[c,m,l,w,s] = deal(propvals{:});

switch(ftype)
 case {'pdf' 'icdf' 'cdf' 'survivor' 'cumhazard'}
   ydensname = genvarname('YPlot',arglist);
   if isequal(ftype,'pdf') && anycontinuous && anydiscrete
      blkf{end+1} = sprintf('%s = %s;',xgridname,xcontname);
   end
   switch(ftype)
       case 'pdf'
           rhs = sprintf('pdf(%s,%s)',outname,xgridname);
       case 'icdf'
           rhs = sprintf('icdf(%s,%s)',outname,xgridname);
       case 'cdf'
           rhs = sprintf('cdf(%s,%s)',outname,xgridname);
       case 'survivor'
           rhs = sprintf('1 - cdf(%s,%s)',outname,xgridname);
       case 'cumhazard'
           rhs = sprintf('-log(1 - cdf(%s,%s))',outname,xgridname);
   end
   blkf{end+1} = sprintf('%s = %s;',ydensname,rhs);
   blkf{end+1} = sprintf('%s = plot(%s,%s,''Color'',[%g %g %g],...',...
                         hlinename,xgridname,ydensname,c(1),c(2),c(3));
   blkf{end+1} = sprintf('          ''LineStyle'',''%s'', ''LineWidth'',%d,...',l,w);
   blkf{end+1} = sprintf('          ''Marker'',''%s'', ''MarkerSize'',%d);',m,s);
  

 case 'probplot'
   blkf{end+1} = sprintf('%s = probplot(gca,%s);',hlinename,outname);
   blkf{end+1} = sprintf('set(%s,''Color'',[%g %g %g],''LineStyle'',''%s'', ''LineWidth'',%d);', ...
         hlinename,c(1),c(2),c(3),l,w);
end

% --------------- get variable names for arg list
function [blkc,dsvarnames,exprlist,arglist] = getarglist(blkc,dsdb)

ds = down(dsdb);
dsnum = 0;

dsvarnames = {};  % each row contains the varnames for one dataset
exprlist = {};    % names and expressions of the data, censoring, frequency
arglist = {};     % variable names to use for each expression

descrtext = {'Y' 'Censoring' 'Frequency'};
while(~isempty(ds))
    % Get y,cens,freq names for each dataset
    dsnum = dsnum + 1;
    dsname = ds.name;
    yname = ds.yname;
    [censname,freqname] = getcensfreqname(ds);
    originalnames = {yname censname freqname};
    
    % Create comment text associating dataset with variable names
    blkc{end+1} = ' ';
    blkc{end+1} = sprintf('%% Data from dataset "%s":',dsname);
    
    % Each non-empty variable name becomes a function argument,
    % except expressions that are not valid variable names have
    % to be replaced by a variable name that we will select here
    for j=1:3
        exprj = originalnames{j};
        if isempty(exprj)
            continue;
        end
        
        % If we know the original expression, look up its index in our list
        exprnum = find(strcmp(exprj,exprlist));
        if isempty(exprnum)
            exprnum = length(exprlist) + 1;
            exprlist{exprnum,1} = exprj;
            if isvarname(exprj)
                namej = exprj;
            else
                namej = sprintf('arg_%d',exprnum);
            end
            arglist{exprnum,1} = namej;
        else
            namej = arglist{exprnum,1};
        end
        
        if isequal(namej,exprj)
            suffix = '';
        else
            exprj = strrep(exprj,'$GENERATED NAME$ ','');
            suffix = sprintf(' (originally %s)',exprj);
        end
        
        originalnames{j} = namej;
        blkc{end+1} = sprintf('%%    %s = %s%s',descrtext{j},namej,suffix);
    end
    
    dsvarnames(dsnum,1:3) = originalnames;
    ds = right(ds);
end


% --------------- write code for plotting data set
function [blkd,arglist,showbounds,onplot] = writedset(blkd,ds,arglist,...
    originalnames,allprop,alpha,hlinename,hboundsname,simplecase)

dsname = ds.name;
ftype = ds.ftype;
showbounds = false;
onplot = true;

yname = originalnames{1};
censname = originalnames{2};
freqname = originalnames{3};

havecens = ~isempty(censname);
havefreq = ~isempty(freqname);

% Create code to plot this dataset into the figure we have created
blkd{end+1} = '';
blkd{end+1} = sprintf('%% --- Plot data originally in dataset "%s"',dsname);
[censname,freqname] = fixcensfreq(censname,freqname);

dsline = ds.line;
if isempty(dsline) || ~ishghandle(dsline)
   blkd{end+1} = '% This dataset does not appear on the plot';
   onplot = false;
   return;
end

propvals = get(dsline,allprop);
[c,m,l,w,s] = deal(propvals{:});

switch(ftype)
 case 'pdf'
   % Generate code to compute the empirical cdf
   Xname = genvarname('CdfX',arglist);
   Fname = genvarname('CdfF',arglist);
   binheightname = genvarname('BinHeight',arglist);
   bincentername = genvarname('BinCenter',arglist);
   binedgename = genvarname('BinEdge',arglist);
   bininfoname = genvarname('BinInfo',arglist);
        
   oneline = sprintf('[%s,%s] = ecdf(%s,''Function'',''cdf''', Fname,Xname,yname);
   if havecens
      oneline = sprintf('%s,''cens'',%s',oneline,censname);
   end
   if havefreq
      oneline = sprintf('%s,''freq'',%s',oneline,freqname);
   end
   blkd{end+1} = sprintf('%s);  %% compute empirical cdf',oneline)';
   
   % Generate code to duplicate the current histogram bin width selection
   bininfo = ds.binDlgInfo;
   if isempty(bininfo)           % use default in case this is empty
      bininfo.rule = 1;
   end
   blkd{end+1} = sprintf('%s.rule = %d;', bininfoname,bininfo.rule);

   switch bininfo.rule
    case 3
      blkd{end+1} = sprintf('%s.nbins = %d;',bininfoname,bininfo.nbins);

    case 5
      blkd{end+1} = sprintf('%s.width = %g;',bininfoname,bininfo.width);
      blkd{end+1} = sprintf('%s.placementRule = %d;',bininfoname,bininfo.placementRule);
      if bininfo.placementRule ~= 1
         blkd{end+1} = sprintf('%s.anchor = %g;',bininfoname,bininfo.anchor);
      end
   end
   
   blkd{end+1} = sprintf('[~,%s] = internal.stats.histbins(%s,%s,%s,%s,%s,%s);',...
                         binedgename,yname,censname,freqname,bininfoname,Fname,Xname);

   % Generate code to produce the histogram
   blkd{end+1} = sprintf('[%s,%s] = ecdfhist(%s,%s,''edges'',%s); % empirical pdf from cdf',...
                         binheightname,bincentername,Fname,Xname,binedgename);
   blkd{end+1} = sprintf('%s = bar(%s,%s,''hist'');',hlinename,bincentername,binheightname);
   blkd{end+1} = sprintf('set(%s,''FaceColor'',''none'',''EdgeColor'',[%g %g %g],...', ...
                         hlinename,c(1),c(2),c(3));
   blkd{end+1} = sprintf('       ''LineStyle'',''%s'', ''LineWidth'',%d);', ...
         l,w);
   blkd{end+1} = 'xlabel(''Data'');';
   blkd{end+1} = 'ylabel(''Density'')';
  
 case {'cdf' 'survivor' 'cumhazard'}
   xcdfname = genvarname('CdfX',arglist);
   ycdfname = genvarname('CdfY',arglist);
   lowername = genvarname('CdfLower',arglist);
   uppername = genvarname('CdfUpper',arglist);

   showbounds = ds.showbounds;
   if showbounds
      oneline = sprintf('[%s,%s,%s,%s] = ecdf(%s,''Function'',''%s'',''alpha'',%g',...
                            ycdfname,xcdfname,lowername,uppername,yname, ftype,alpha);
   else
      oneline = sprintf('[%s,%s] = ecdf(%s,''Function'',''%s''',...
                            ycdfname,xcdfname,yname, ftype);
   end
   oneline = addcensfreq(oneline,censname,freqname);
   blkd{end+1} = sprintf('%s);  %% compute empirical function',oneline);
   blkd{end+1} = sprintf('%s = stairs(%s,%s,''Color'',[%g %g %g],''LineStyle'',''%s'', ''LineWidth'',%d);', ...
         hlinename,xcdfname,ycdfname,c(1),c(2),c(3),l,w);
   if showbounds
      xstairslowername = genvarname('StairsXlower',arglist);
      ystairslowername = genvarname('StairsYlower',arglist);
      blkd{end+1} = sprintf('[%s,%s] = stairs(%s,%s);',xstairslowername,ystairslowername,xcdfname,lowername);
      xstairsuppername = genvarname('StairsXupper',arglist);
      ystairsuppername = genvarname('StairsYupper',arglist);
      blkd{end+1} = sprintf('[%s,%s] = stairs(%s,%s);',xstairsuppername,ystairsuppername,xcdfname,uppername);
      blkd{end+1} = sprintf('%s = plot([%s(:); NaN; %s(:)], [%s(:); NaN; %s(:)],...',...
                            hboundsname,xstairslowername,xstairsuppername,ystairslowername,ystairsuppername);
      blkd{end+1} = sprintf('   ''Color'',[%g %g %g],''LineStyle'','':'', ''LineWidth'',1);', ...
         c(1),c(2),c(3));
   end
   blkd{end+1} = 'xlabel(''Data'');';
   switch(ftype)
      case 'cdf'
          blkd{end+1} = 'ylabel(''Cumulative probability'')';
      case 'survivor'
          blkd{end+1} = 'ylabel(''Survivor function'')';
      case 'cumhazard'
          blkd{end+1} = 'ylabel(''Cumulative hazard'')';
   end

 case 'icdf'
   xcdfname = genvarname('CdfX',arglist);
   ycdfname = genvarname('CdfY',arglist);
   oneline = sprintf('[%s,%s] = ecdf(%s,''Function'',''cdf''', ycdfname,xcdfname,yname);
   oneline = addcensfreq(oneline,censname,freqname);
   blkd{end+1} = sprintf('%s);  %% compute empirical cdf',oneline);
   blkd{end+1} = sprintf('%s = stairs(%s,[%s(2:end);%s(end)],''Color'',[%g %g %g],''LineStyle'',''%s'', ''LineWidth'',%d);', ...
         hlinename,ycdfname,xcdfname,xcdfname,c(1),c(2),c(3),l,w);
   blkd{end+1} = 'xlabel(''Probability'');';
   blkd{end+1} = 'ylabel(''Quantile'')';

 case 'probplot'
   [censname,freqname] = fixcensfreq(censname,freqname);
   if simplecase
       dtype = dfgetset('dtype');
       blkd{end+1} = sprintf('%s = probplot(''%s'',%s,%s,%s,''noref'');', ...
                         hlinename,dtype,yname,censname,freqname);
   else
       blkd{end+1} = sprintf('%s = probplot(gca,%s,%s,%s,''noref''); %% add data to existing plot', ...
                         hlinename,yname,censname,freqname);
   end
   blkd{end+1} = sprintf('set(%s,''Color'',[%g %g %g],''Marker'',''%s'', ''MarkerSize'',%d);', ...
         hlinename,c(1),c(2),c(3),m,s);
   blkd{end+1} = 'xlabel(''Data'');';
   blkd{end+1} = 'ylabel(''Probability'')';
end


% -----------------------------
function [blkf,yname,censname,freqname]=applyexclusion(blkf,exclrule,...
                                          yname,censname,freqname,arglist)
%APPLYEXCLUSION Change var names to use indexing to apply exclusion rule

% Create expressions for inclusion rules
if isempty(exclrule.ylow)
   e1 = '';
else
   ylow = str2double(exclrule.ylow);
   if exclrule.ylowlessequal==0
      e1 = sprintf('%s > %g', yname, ylow);
   else
      e1 = sprintf('%s >= %g', yname, ylow);
   end
end
if isempty(exclrule.yhigh)
   e2 = '';
else
   yhigh = str2double(exclrule.yhigh);
   if exclrule.yhighgreaterequal==0
      e2 = sprintf('%s < %g', yname, yhigh);
   else
      e2 = sprintf('%s <= %g', yname, yhigh);
   end
end

% Combine exclusion expressions
if isempty(e1)
   if isempty(e2)
      etxt = '';
   else
      etxt = e2;
   end
else
   if isempty(e2)
      etxt = e1;
   else
      etxt = sprintf('%s & %s',e1,e2);
   end
end

% Create code to generate index vector and reduce all variables
if ~isempty(etxt)
   exclname = genvarname('Excluded',arglist);
   blkf{end+1} = sprintf('\n%% Create vector for exclusion rule ''%s''',...
                         exclrule.name);
   blkf{end+1} =         '% Vector indexes the points that are included';
   blkf{end+1} = sprintf('%s = (%s);\n', exclname,etxt);

   newyname = genvarname('Data',arglist);
   blkf{end+1} = sprintf('%s = %s(%s);',newyname,yname,exclname);
   yname = newyname;
   if ~isempty(censname) && ~isequal(censname,'[]')
      newcensname = genvarname('Cens',arglist);
      blkf{end+1} = sprintf('%s = %s(%s);',newcensname,censname,exclname);
      censname = newcensname;
   end
   if ~isempty(freqname) && ~isequal(freqname,'[]')
      newfreqname = genvarname('Freq',arglist);
      blkf{end+1} = sprintf('%s = %s(%s);',newfreqname,freqname,exclname);
      freqname = newfreqname;
   end
end

% -----------------------------------------
function [censname,freqname] = getcensfreqname(ds,exprlist,arglist)
%GETCENSFREQNAME Get censoring and frequency names

censname = ds.censname;
freqname = ds.freqname;
if strcmp(censname,'(none)')
   censname = '';
end
if strcmp(freqname,'(none)')
   freqname = '';
end

if ~isvarname(censname) && ~isempty(ds.censored)
   % We have a censoring expression, so create a fake non-empty name
   censname = sprintf('$GENERATED NAME$ %s %s',ds.name,'censored');
end
if ~isvarname(freqname) && ~isempty(ds.frequency)
   % We have a frequency expression, so create a fake non-empty name
   freqname = sprintf('$GENERATED NAME$ %s %s',ds.name,'frequency');
end

if nargin>=3
   censname = expression2name(censname,exprlist,arglist);
   freqname = expression2name(freqname,exprlist,arglist);
end
   
  
% -------------------------------------------
function nm = expression2name(expr,exprlist,arglist)
%EXPRESSION2NAME Find out what name we're using in place of this expression

nm = expr;
if ~isempty(expr)
   j = find(strcmp(expr,exprlist));
   if isscalar(j)
      nm = arglist{j};
   end
end


% -------------------------------------------
function [censname,freqname] = fixcensfreq(censname,freqname)
%FIXCENSFREQ Fix empty cens/freq values for use in code

if isempty(censname)
   censname = '[]';
end
if isempty(freqname)
   freqname = '[]';
end

% ------------------------------------------
function pdprefix = getpdprefix(arglist)
% Get a prefix to use for the probdist var names
pdprefix = 'pd';
while(true)
   t = strncmp(pdprefix,arglist,length(pdprefix));
   if ~any(t) % no vars start with this prefix, okay
       return
   end
   
   ok = true;
   tfind = find(t);
   n = length(pdprefix);
   for j=1:length(tfind)
       argj = arglist{tfind(j)};
       if length(argj)>n && ismember(argj(n+1),'123456789')
           ok = false; % found a var with this prefix, then a number
           break
       end
   end
   if ok
       return
   end
   pdprefix = [pdprefix '_'];
end
