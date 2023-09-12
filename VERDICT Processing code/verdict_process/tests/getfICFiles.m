function [Aref, B] = getfICFiles(testCase, outputFolder, resultsFolderName, fICSeriesNumber)
% GETFICFILES Get the files with fIC values. Used in testing
%
% [Aref, B] = getfICFiles(testCase, outputFolder, resultsFolderName, fICSeriesNumber)
%

dreffIC = dfparse(getAllFiles(testCase.reffICFolder)) ;
[vreffIC, mref] = d2mat(dreffIC,{'slice'},'op','dv') ;

doutput = dfparse(getAllFiles(fullfile(outputFolder,resultsFolderName,'DICOM'))) ;
[voutfIC, mop] = d2mat(doutput,{'slice','series'},'series', ...
    fICSeriesNumber, 'op','dv') ;

Aref = vreffIC(testCase.evalCoords{:}) ;
B    = voutfIC(testCase.evalCoords{:}) ;

end

