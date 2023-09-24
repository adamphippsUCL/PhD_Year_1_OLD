% MATLAB script to run VERDICT processing on specified patients

% Define list of patient numbers
PatNums = [...
%     "BAR_003",...
%     "BAR_005",...
%     "BAR_006",...
%     "BAR_009",...
%     "BAR_004",...
    "BAR_033",...
    "INN_019",...
%     "BAR_012"
    ];


% Model type
ModelType = "Original";


% Run VERDICT processing
for PatNum = PatNums
    
    VERDICT(PatNum, ModelType);
end