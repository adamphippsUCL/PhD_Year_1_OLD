classdef vresampleTest < matlab.unittest.TestCase
    % vresampleTest tests vresample
    %
    % Copyright, 2020, David Atkinson
    % D.Atkinson@ucl.ac.uk
    
    
    methods (Test)
        function test(testCase)
            
            geom1.IPP = [0 0 0];
            geom1.IOP = [1 0 0 0 1 0] ;
            geom1.PixelSpacing_HW = [1 1] ;
            geom1.Height = 128 ;
            geom1.Width = 128 ;
            geom1.SliceThickness = 2 ;
            
            [vol, geom] = analytic_phantom(geom1, 30) ;
            
            geom_resamp = geom_change_ps(geom, [2 2]) ;
            
            Vresamp = vresample(vol, geom, geom_resamp) ;
            
            [Vresz,geomt] = dinplanet(vol, geom,'imresize',[64 64]) ;
            
            testCase.verifyEqual(double(Vresz),double(Vresamp), 'RelTol',0.0001, 'AbsTol',1e-4);
            
            testCase.verifyEqual([geom_resamp.IPP], [geomt.IPP], 'RelTol', 0.001, 'AbsTol',1e-7) ;
            
            testCase.verifyEqual([geom_resamp.IOP], [geomt.IOP], 'RelTol', 0.001, 'AbsTol',1e-7) ;
           
            testCase.verifyEqual([geom_resamp.PixelSpacing_HW], ...
                [geomt.PixelSpacing_HW], 'RelTol', 0.001, 'AbsTol',1e-7) ;
        end
    end
end

      