classdef testverdict_process < matlab.unittest.TestCase
    %TESTVERDICT_PROCESS Test for verdict processing code
    %   
    % To test an individual case, e.g.
    %    runtests('testverdict_process/testGESigna')
    %
    %  Copyright 2023. David Atkinson. University College London
    %
    % See Also compfig getfICFiles verdict verdict_process
    
    methods(TestClassSetup)
        function prepareData(testCase)
            testCase.dataFolderRoot = '/Users/davidatkinson/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data/VERDICT_test_DICOMs' ;
            testCase.dataFolder   = fullfile(testCase.dataFolderRoot,'invivo','INN065','INN-065-pXNAT/data') ;
            testCase.reffICFolder = fullfile(testCase.dataFolder,'DICOM/INN065 [INN065]/20160725 162419 [ - INNOVATE]/Series 802 [MR - fIC  with RescaleSlope]' );
            testCase.evalCoords = {60:120, 60:120, 7};

            testCase.PhilipsClassicFolder = fullfile(testCase.dataFolderRoot,'phantoms/PhilipsClassic' );
            testCase.PhilipsEnhancedFolder = fullfile(testCase.dataFolderRoot,'phantoms/PhilipsEnhanced') ;
            testCase.GESignaFolder = fullfile(testCase.dataFolderRoot,'phantoms/GESigna') ;
        end
    end

    properties
        dataFolderRoot
        dataFolder
        reffICFolder
        evalCoords 

        PhilipsClassicFolder
        PhilipsEnhancedFolder
        GESignaFolder
    end

    % import matlab.unittest.constraints.IsEqualTo
    % import matlab.unittest.constraints.IsLessThan
    % import matlab.unittest.constraints.AbsoluteTolerance
    % import matlab.unittest.constraints.RelativeTolerance

    methods(Test)
        function testpXNAT(testCase)
            % Pre-registered data with most similar processing to
            % prostate-XNAT. MATLAB only (not Docker due to SPAMS)
            
            import matlab.unittest.constraints.IsLessThan

            resultsFolderName = 'res-pXNAT0065' ;

            fICReconNumber  = 5  ;  % DICOM fIC recon number for saved fIC 
            vBaseSeriesNumber = 8 ; % DICOM Series number for base Series

            fICSeriesNumber = 100*vBaseSeriesNumber + fICReconNumber ;

            outputFolder = tempdir ;

            args = { ...
                'quiet', 'true', ...
                'register', 'false', ...
                'swapinvXNAT', 'true', ...
                'usedirecdiff', 'true', ...
                'forceXNATscheme', 'true' , ...
                'solver', 'SPAMS', ...
                'fICReconNumber', num2str(fICReconNumber), ...
                'vBaseSeriesNumber', num2str(vBaseSeriesNumber) } ;

            diff_rng = [-0.005 0.005] ;
            maxdiff  = 0.01 ;

            verdict(testCase.dataFolder, outputFolder, args{:}, ...
                'resultsFolderName', resultsFolderName)

            [Aref, B] = getfICFiles(testCase, outputFolder, ...
                resultsFolderName, fICSeriesNumber) ;

            compfig(Aref, B, [0 1], diff_rng)

            resdiff = Aref - B ;
            
            testCase.verifyThat(max(abs(resdiff(:))) , IsLessThan(maxdiff) ) 

        end

        function testpXNATlsq( testCase )
            % As above but using lsqnonneg optimisation
            import matlab.unittest.constraints.IsLessThan

            resultsFolderName = 'res-pXNAT0065' ;

            fICReconNumber  = 6  ;  % DICOM fIC recon number for saved fIC 
            vBaseSeriesNumber = 8 ; % DICOM Series number for base Series

            fICSeriesNumber = 100*vBaseSeriesNumber + fICReconNumber ;

            outputFolder = tempdir ;

            args = { ...
                'quiet', 'true', ...
                'register', 'false', ...
                'swapinvXNAT', 'true', ...
                'usedirecdiff', 'true', ...  
                'forceXNATscheme', 'true' , ...
                'solver', 'lsqnonnegTikonhov', ...
                'fICReconNumber', num2str(fICReconNumber) , ...
                'vBaseSeriesNumber', num2str(vBaseSeriesNumber)} ;

            diff_rng = [-0.1 0.1] ;
            maxdiff  = 0.5 ;

            verdict(testCase.dataFolder, outputFolder, args{:}, ...
                'resultsFolderName', resultsFolderName)

            [Aref, B] = getfICFiles(testCase, outputFolder, ...
                resultsFolderName, fICSeriesNumber) ;

            compfig(Aref, B, [0 1], diff_rng)
            resdiff = Aref - B ;
            
            testCase.verifyThat(max(abs(resdiff(:))) , IsLessThan(maxdiff) ) 
        end

        

        function testpipeline1( testCase )
            % Full MATLAB pipeline (inc reg. and opt.)
            import matlab.unittest.constraints.IsLessThan

            resultsFolderName = 'res-pXNAT0065' ;
            fICReconNumber = 7 ;
            vBaseSeriesNumber = 8 ; % DICOM Series number for base Series

            fICSeriesNumber = 100*vBaseSeriesNumber + fICReconNumber ;

            outputFolder = tempdir ;

            args = { ...
                'quiet', 'true', ...
                'fICReconNumber', num2str(fICReconNumber) , ...
                'vBaseSeriesNumber', num2str(vBaseSeriesNumber) } ;
                

            diff_rng = [-0.5 0.5] ;
            maxdiff  = 1.1 ;

            verdict(testCase.dataFolder, outputFolder, args{:}, ...
                'resultsFolderName', resultsFolderName)

            [Aref, B] = getfICFiles(testCase, outputFolder, ...
                resultsFolderName, fICSeriesNumber) ;

            compfig(Aref, B, [0 1], diff_rng)
            resdiff = Aref - B ;

            testCase.verifyThat(max(abs(resdiff(:))) , IsLessThan(maxdiff) )
        end

        function testClassicEnhancedPhilips( testCase )
            % Full MATLAB pipeline on Classic and Enhanced and writing for
            % Philips scanner
            import matlab.unittest.constraints.IsLessThan

            resultsClassicFolderName = 'res-Classic' ;
            resultsEnhancedFolderName = 'res-Enhanced' ;

            fICReconNumber = 8 ;
            vBaseSeriesNumber = 9 ; % DICOM Series number for base Series

            fICSeriesNumber = 100*vBaseSeriesNumber + fICReconNumber ;

            outputFolder = tempdir ;

            args = { ...
                'quiet', 'true', ...
                'fICReconNumber', num2str(fICReconNumber) , ...
                'vBaseSeriesNumber', num2str(vBaseSeriesNumber) } ;
                

           verdict(testCase.PhilipsClassicFolder, outputFolder, args{:}, ...
                'resultsFolderName', resultsClassicFolderName)

            verdict(testCase.PhilipsEnhancedFolder, outputFolder, args{:}, ...
                'resultsFolderName', resultsEnhancedFolderName, ...
                'addPhilipsPrivate', true)

            AClassic  = load(fullfile(outputFolder,resultsClassicFolderName,'outputs.mat') ,'fIC') ;
            BEnhanced = load(fullfile(outputFolder,resultsEnhancedFolderName,'outputs.mat') ,'fIC') ;

            AfIC = AClassic.fIC(testCase.evalCoords{:}) ;
            BfIC = BEnhanced.fIC(testCase.evalCoords{:}) ;

            AfIC(~isfinite(AfIC))=0;
            BfIC(~isfinite(BfIC))=0;

            diff_rng = [-0.2 0.2] ;
            maxdiff  = 0.1 ;

            compfig(AfIC, BfIC, [0 1], diff_rng)
            resdiff =  AfIC- BfIC;

            testCase.verifyThat(max(abs(resdiff(:))) , IsLessThan(maxdiff) )

            % Also test that .mat is same as DICOM
            % Note the DICOM has had fIC clipped 
            doutput = dfparse(getAllFiles(fullfile(outputFolder,resultsEnhancedFolderName,'DICOM'))) ;
            [voutfIC, mop] = d2mat(doutput,{'slice','series'},'series', ...
              fICSeriesNumber, 'op','dv') ;

            BfIC(BfIC>1)=1;
            BfIC(BfIC<0)=0;

            compfig(BfIC, voutfIC(testCase.evalCoords{:}), [0 1], diff_rng)
            resdiff =  BfIC - voutfIC(testCase.evalCoords{:}) ;

            maxdiff  = 1.1 ;
            testCase.verifyThat(max(abs(resdiff(:))) , IsLessThan(maxdiff) )

        end

        function testGESigna( testCase )
            import matlab.unittest.constraints.IsLessThan

            resultsFolderName = 'res-GESigna' ;
            fICReconNumber = 2 ;
            vBaseSeriesNumber = 3 ; % DICOM Series number for base Series

            fICSeriesNumber = 100*vBaseSeriesNumber + fICReconNumber ;

            outputFolder = tempdir ;

            args = { ...
                'quiet', 'true', ...
                'fICReconNumber', num2str(fICReconNumber) , ...
                'vBaseSeriesNumber', num2str(vBaseSeriesNumber) } ;

            verdict(testCase.GESignaFolder, outputFolder, args{:}, ...
                'resultsFolderName', resultsFolderName, 'bvSortDirection', 'ascend')

            Dthis = load(fullfile(outputFolder,resultsFolderName,'outputs.mat'));

            Dexpected = load(fullfile(testCase.GESignaFolder,'GESigna_outputs.mat')) ;

            resdiff = Dthis.fIC - Dexpected.fIC ;
            maxdiff = 0.02 ; % can increase if methodology changes.
            testCase.verifyThat(max(abs(resdiff(:))) , IsLessThan(maxdiff) )


        end

        function testDocker( testCase )
            % Full MATLAB pipeline (inc reg. and opt.) vs Docker
            import matlab.unittest.constraints.IsLessThan

            resultsFolderName = 'res-065' ;
            resultsFolderNameDocker = 'res-065-Docker' ;
            fICReconNumber = 8 ;
            vBaseSeriesNumber = 8 ; % DICOM Series number for base Series

            fICSeriesNumber = 100*vBaseSeriesNumber + fICReconNumber ;

            outputFolder = tempdir ;

            args = { ...
                'quiet', 'true', ...
                'fICReconNumber', num2str(fICReconNumber), ...
                'vBaseSeriesNumber', num2str(vBaseSeriesNumber) } ;

            diff_rng = [-0.005 0.005] ;
            maxdiff  = 0.01 ;

            % Docker
            status = system('xhost +') ;
            if status ~= 0
                warning('No X detected, figures will not display')
            end

            docker_exec = '/usr/local/bin/docker' ;
            docker_cmd = [ docker_exec, ' run -d --rm -e "DISPLAY=host.docker.internal:0" ', ...
                '-v ''/tmp/.X11-unix'':''/tmp/.X11-unix'' ', ...
                '-v ''',testCase.dataFolder, ''':/home/appuser/dfolder ', ...
                '-v ''',outputFolder, ''':/home/appuser/outputfolder ', ...
                'verdict /home/appuser/dfolder /home/appuser/outputfolder ' ] ;
            for iarg = 1:length(args)
                docker_cmd = [docker_cmd, args{iarg}, ' '] ;
            end
            docker_cmd = [docker_cmd, 'resultsFolderName ', ...
                resultsFolderNameDocker] 

            status = system(docker_cmd) ;
            if status ~= 0
                error('Docker stats non-zero')
            end

            txt = input("Press Enter when Docker has finished.","s") ;

             % MATLAB
            verdict(testCase.dataFolder, outputFolder, args{:}, ...
                'resultsFolderName', resultsFolderName)


            [Aref, B] = getfICFiles(testCase, outputFolder, ...
                resultsFolderName, fICSeriesNumber) ;
            
            [Aref, BDocker] = getfICFiles(testCase, outputFolder, ...
                resultsFolderNameDocker, fICSeriesNumber) ;

            compfig(B, BDocker, [0 1], diff_rng)
            resdiff = BDocker - B ;

            testCase.verifyThat(max(abs(resdiff(:))) , IsLessThan(maxdiff) )
        end

end



end

