function dfexport2workspace()
% DFEXPORT2WORKSPACE Helper file for Distribution Fitter
%    DFEXPORT2WORKSPACE gets the saved evaluated results and passes them to
%    export2wsdlg

%   Copyright 2003-2016 The MathWorks, Inc.


results = dfgetset('evaluateResults');
export2wsdlg({getString(message('stats:dfstrings:cellstr_SaveToVariable'))}, ...
             {'evaluateresults'}, {results});
