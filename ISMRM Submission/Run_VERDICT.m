% MATLAB script to run VERDICT processing on specified patients

% Define list of patient numbers
PatNums = [...
    "BAR_003",...
%     "BAR_005",...
%     "BAR_006",...
%     "BAR_009",...
%     "BAR_004",...
%     "BAR_033",...
%     "INN_019",...
%     "BAR_012"
%     "INN_087",...
%     "INN_099",...
%     "INN_109",...
%     "BAR_035",...
%     "BAR_017",...
%     "BAR_028"
    ];


%% Original model 

% Model type
ModelType = "Original";
solver = 'lsqnonnegTikonhov';

% Run VERDICT processing
for PatNum = PatNums
    disp(["----------->" PatNum])
    VERDICT(PatNum, ModelType, solver = solver);
end

%% No VASC model

% Model type
ModelType = "No VASC";
solver = 'lsqnonnegTikonhov';

% Run VERDICT processing
for PatNum = PatNums
    disp(["----------->" PatNum])
    VERDICT(PatNum, ModelType, solver = solver);
end

%% No VASC Reduced Rs 1

% Model type
ModelType = "No VASC Reduced Rs 1";
solver = 'lsqnonneg';

% Run VERDICT processing
for PatNum = PatNums
    disp(["----------->" PatNum])
    VERDICT(PatNum, ModelType, solver = solver);
end

%% No VASC Reduced Rs 2

% Model type
ModelType = "No VASC Reduced Rs 2";
solver = 'lsqnonneg';

% Run VERDICT processing
for PatNum = PatNums
    disp(["----------->" PatNum])
    VERDICT(PatNum, ModelType, solver = solver);
end

%% No VASC Reduced Rs 3

% Model type
ModelType = "No VASC Reduced Rs 3";
solver = 'lsqnonneg';

% Run VERDICT processing
for PatNum = PatNums
    disp(["----------->" PatNum])
    VERDICT(PatNum, ModelType, solver = solver);
end


%% No VASC Reduced Rs 4

% Model type
ModelType = "No VASC Reduced Rs 4";
solver = 'lsqnonneg';

% Run VERDICT processing
for PatNum = PatNums
    disp(["----------->" PatNum])
    VERDICT(PatNum, ModelType, solver = solver);
end