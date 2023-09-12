function compare_results
% COMPARE_RESULTS Compare VERDICT processing results visually
%
% Both HMU_008 and HMU_026 below show very samll differences
%
%
% David Atkinson
%
% See also verdict_fit

% HMU_008
refPath = '/Users/davidatkinson/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data/HMU/HMU_008/DENIFTIFIED' ;
currPath = '/Users/davidatkinson/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data/HMU/HMU_008/res-2023-03-19T172808' ;

% HMU_026
refPath = '/Users/davidatkinson/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data/HMU/HMU_026/assessors' ;
currPath = '/Users/davidatkinson/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data/HMU/HMU_026/res-2023-03-19T173547' ;


refFile = 'ficXNAT.mat' ;
refVar = 'ficXNAT' ;

currFile = 'outputs.mat' ;
currVar = 'fIC' ;

refFFN = fullfile(refPath, refFile) ;
currFFN = fullfile(currPath, currFile) ;

refdat  = load(refFFN) ;
currdat = load(currFFN) ;

refval = refdat.(refVar) ;
currval = currdat.(currVar) ;

sviewer(refval-currval, [], CLim=[-0.1 0.1], Name='difference')
