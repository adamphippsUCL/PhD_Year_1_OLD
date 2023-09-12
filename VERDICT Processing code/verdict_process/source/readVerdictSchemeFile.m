function [version, data] = readVerdictSchemeFile(fileName)
% READVERDICTSCHEMEFILE Read a text file and return the version number and data as a MATLAB table
%   [version, data] = READVERDICTSCHEMEFILE(fileName) reads the text file specified by fileName and
%   returns the version number and data as a MATLAB table. The first line of the text file
%   should be in the format "VERSION: X" where X is the version number. The remaining lines
%   of the file should contain data in columns separated by whitespace. The fifth column
%   should be labeled "DELTA" and the sixth column should be labeled "delta".
%
%   Inputs:
%       fileName - the name of the text file to read (string)
%
%   Outputs:
%       version - the version number read from the file (integer)
%       data - a MATLAB table containing the data read from the file

% Open the file for reading
fileID = fopen(fileName, 'r');

% Read the first line of the file to get the version number
versionLine = fgetl(fileID);
versionTokens = split(versionLine, ':');
version = str2double(versionTokens{2});

% Read the rest of the file into a MATLAB table
data = readtable(fileName);

% Rename the fifth and sixth columns
data.Properties.VariableNames{5} = 'DELTA';
data.Properties.VariableNames{6} = 'delta';

% Close the file
fclose(fileID);
end