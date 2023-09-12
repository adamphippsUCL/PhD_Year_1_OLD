function explore_test(exec)
% Temporary function to explore set-up of test for verdict_process and
% Docker etc

dataFolderRoot = '/Users/davidatkinson/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data/VERDICT_test_DICOMs' ;

% Using pXNAT pre-registered data, check new fIC DICOMs match those saved from assessors

dataFolder   = fullfile(dataFolderRoot,'invivo','INN065','INN-065-pXNAT/data') ;
reffICFolder = fullfile(dataFolder,'DICOM/INN065 [INN065]/20160725 162419 [ - INNOVATE]/Series 802 [MR - fIC  with RescaleSlope]' );

evalCoords = {60:120, 60:120, 7};

resultsFolderName = 'res-pXNAT0065' ;
fICSeriesNumber = 56 ;

outputFolder = tempdir ;

args = { ...
    'register', 'false', ...
    'swapinvXNAT', 'true', ...
    'usedirecdiff', 'true', ...
    'forceXNATscheme', 'true' , ...
    'solver', 'SPAMS', ...
    'fICSeriesNumber', num2str(fICSeriesNumber) } ;

diff_rng = [-0.005 0.005] ;


switch exec
    case 'MATLAB'
        % MATLAB
        verdict(dataFolder, outputFolder, args{:}, 'resultsFolderName', resultsFolderName)
    case 'Docker'
        % DOCKER
        resultsFolderName = [resultsFolderName,'-Docker'] ;

        status = system('xhost +') ;
        if status ~= 0
            warning('No X detected, figures will not display')
        end

        docker_cmd = [ 'docker run -d --rm -e "DISPLAY=host.docker.internal:0" ', ...
            '-v ''/tmp/.X11-unix'':''/tmp/.X11-unix'' ', ...
            '-v ''',dataFolder, ''':/home/appuser/dfolder ', ...
            '-v ''',outputFolder, ''':/home/appuser/outputfolder ', ...
            'verdict /home/appuser/dfolder /home/appuser/outputfolder ' ] ;
        for iarg = 1:length(args)
            docker_cmd = [docker_cmd, args{iarg}, ' '] ;
        end
        docker_cmd = [docker_cmd, 'resultsFolderName ', resultsFolderName] 

        status = system(docker_cmd) ;
        if status ~= 0
            warning('Docker stats non-zero')
        end
    otherwise
        error('Unknown exec')
end
% read in pXNAT fIC
% read in testFIC from DICOM
% Compare

dreffIC = dfparse(getAllFiles(reffICFolder)) ;
[vreffIC, mref] = d2mat(dreffIC,{'slice'},'op','dv') ;

doutput = dfparse(getAllFiles(fullfile(outputFolder,resultsFolderName,'DICOM'))) ;
[voutfIC, mop] = d2mat(doutput,{'slice','series'},'series', ...
    fICSeriesNumber, 'op','dv') ;

% set box region and slice,
% evaluate largest absolute error 

compfig(vreffIC(evalCoords{:}), voutfIC(evalCoords{:}), [0 1], diff_rng)


disp('At end')

end




