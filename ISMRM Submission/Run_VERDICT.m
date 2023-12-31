% MATLAB script to run VERDICT processing on specified patients

% Define list of patient numbers
PatNums = [...
%     "BAR_006",...
%     "BAR_010",...
%     "BAR_012",...
%     "BAR_019",...
%     "BAR_023",...
%     "BAR_025",...
%     "BAR_030",...
%     "BAR_033",...
%     "BAR_034",...
%     "BAR_035",...
    "BAR_038",...
%     "BAR_039",...
%     "BAR_040",...
%     "BAR_046",...
%     "BAR_047",...
%     "BAR_049",...
%     "BAR_050",...
%     "BAR_052"


    ];

% 
%% Model 1: Original model 

% Model type
ModelType = "Model 1";
solver = 'lsqnonnegTikonhov';

% Run VERDICT processing
for PatNum = PatNums
    disp(["----------->" PatNum])
    VERDICT(PatNum, ModelType, solver = solver);
end

% %% Model 2: No VASC model
% 
% % Model type
% ModelType = "Model 2";
% solver = 'lsqnonnegTikonhov';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% %% Model 3: No VASC Reduced Rs [0.1,5.1,10.1,15.1]
% 

% % Model type
% ModelType = "Model 3";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end

% % Model 4: No VASC Reduced Rs [3,6,9,12]

% % Model type
% ModelType = "Model 4";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% % 
% %% Model 5: No VASC Reduced Rs [6,9,12,15]
% 
% % Model type
% ModelType = "Model 5";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% 
% %% Model 6: No VASC Reduced Rs [0.1,3.1,6.1,9.1]
% 
% % Model type
% ModelType = "Model 6";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% 
% %% Model 7: No VASC, No 3000, Reduced Rs [0.1,7.6,15.1]
% 
% % Model type
% ModelType = "Model 7";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% 
% %% Model 8: No VASC, No 3000, Reduced Rs [3.75, 7.5, 11.25]
% 
% % Model type
% ModelType = "Model 8";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% %% Model 9: No VASC, No 3000, Reduced Rs [7.5, 11.25, 15]
% 
% % Model type
% ModelType = "Model 9";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% %% Model 10: No VASC, No 3000, Reduced Rs [0.1, 3.85, 7.6]
% 
% % Model type
% ModelType = "Model 10";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end


% % Model type
% ModelType = "Model 11";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% % Model type
% ModelType = "Model 12";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% % Model type
% ModelType = "Model 13";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end
% 
% % Model type
% ModelType = "Model 14";
% solver = 'lsqnonneg';
% 
% % Run VERDICT processing
% for PatNum = PatNums
%     disp(["----------->" PatNum])
%     VERDICT(PatNum, ModelType, solver = solver);
% end