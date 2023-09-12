function fileNames = getAllFiles(fileAndFolderNames, opts)
% GETALLFILES Returns cell array of all files in input list of files and
% folders.
%
% Inputs:
%   - fileAndFolderNames: a cell array of file and/or folder names
%   - Name, Value       : includeSubFolders {true}, minBytes {132}
%
% Outputs:
%   - fileNames: a cell array of full file names
%


arguments
    fileAndFolderNames  % cell array of file and/or folder names
    opts.includeSubFolders = true
    opts.minBytes = 132
end

if ~iscell(fileAndFolderNames)
    fileAndFolderNames = {fileAndFolderNames} ;
end


fileNames = {} ;
ifn = 0 ; % running total of file numbers

% Loop over entries in fileAndFolderNames
for ient = 1:numel(fileAndFolderNames)
    thisEntry = fileAndFolderNames{ient} ;

    if isfolder(thisEntry)
        if opts.includeSubFolders
            dstruc = dir([thisEntry, filesep, '**', filesep, '*']) ;
        else
            dstruc = dir(thisEntry) ;
        end

        for idstruc = 1:numel(dstruc)
            if dstruc(idstruc).isdir
                % ignore (if sub-folders, they will be in dstruc)
            elseif dstruc(idstruc).bytes> opts.minBytes
                ifn=ifn+1;
                fileNames{ifn}  = fullfile(dstruc(idstruc).folder, dstruc(idstruc).name) ;
            end
        end
    else % Entry was a file
        ifn=ifn+1;
        fileNames{ifn} = thisEntry ;
    end
end

