function createPCode
%createPCode Creates p file in psource from source
% Currently only verdict_process

project = matlab.project.currentProject;
currd = pwd ;

cd(fullfile(project.RootFolder,'psource'))

pcode(fullfile(project.RootFolder,'source','verdict_process.m'))

cd(currd)

end