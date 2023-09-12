classdef dinplanetTest < matlab.unittest.TestCase
    % dinplanetTest tests dinplanet code
    %
    % Copyright, 2020, David Atkinson
    % D.Atkinson@ucl.ac.uk
    
    
    methods (Test)
        function testCircshift(testCase)
           import matlab.unittest.fixtures.PathFixture
            DICOM_Folder = get_data_location('DICOM_20191211geometryIngenia','dir') ;
            f = testCase.applyFixture(PathFixture(DICOM_Folder));
            
            dinfo = datparse('IM_0008') ;
            [vold, md] = d2mat(dinfo, {'slice'}, 'op', 'dv') ;
                
            % round trip circshift
            
            [datat, geomt ] = dinplanet( vold, md.geom, 'circshift', [3 9] ) ;
            [datat, geomt ] = dinplanet( datat, geomt, 'circshift', [-3 -9] ) ;
            
            testCase.verifyEqual(double(datat),double(vold), 'RelTol',0.0001, 'AbsTol',1e-7);
            
            testCase.verifyEqual([geomt.IPP], [md.geom.IPP], 'RelTol', 0.001, 'AbsTol',1e-7) ;
                
            testCase.verifyEqual([geomt.IOP], [md.geom.IOP], 'RelTol', 0.001, 'AbsTol',1e-5) ;
                
            testCase.verifyEqual([geomt.PixelSpacing_HW], ...
                    [md.geom.PixelSpacing_HW], 'RelTol', 0.001, 'AbsTol',1e-7) ;
                   
        end
        
        function testRot(testCase)
           import matlab.unittest.fixtures.PathFixture
            DICOM_Folder = get_data_location('DICOM_20191211geometryIngenia','dir') ;
            f = testCase.applyFixture(PathFixture(DICOM_Folder));
            
            dinfo = datparse('IM_0008') ;
            [vold, md] = d2mat(dinfo, {'slice'}, 'op', 'dv') ;
                
            % round trip rotate
            
            [datat, geomt ] = dinplanet( vold, md.geom, 'rot', 1 ) ;
            [datat, geomt ] = dinplanet( datat, geomt, 'rot', -1 ) ;
            
            testCase.verifyEqual(double(datat),double(vold), 'RelTol',0.0001, 'AbsTol',1e-7);
            
            testCase.verifyEqual([geomt.IPP], [md.geom.IPP], 'RelTol', 0.001, 'AbsTol',1e-7) ;
                
            testCase.verifyEqual([geomt.IOP], [md.geom.IOP], 'RelTol', 0.001, 'AbsTol',1e-5) ;
                
            testCase.verifyEqual([geomt.PixelSpacing_HW], ...
                    [md.geom.PixelSpacing_HW], 'RelTol', 0.001, 'AbsTol',1e-7) ;
                   
        end
        
        function testResize(testCase)
           import matlab.unittest.fixtures.PathFixture
            DICOM_Folder = get_data_location('DICOM_20191211geometryIngenia','dir') ;
            f = testCase.applyFixture(PathFixture(DICOM_Folder));
            
            dinfo = datparse('IM_0008') ;
            [vold, md] = d2mat(dinfo, {'slice'}, 'op', 'dv') ;
                
            % round trip resize
            
            [datat, geomt ] = dinplanet( vold, md.geom, 'imresize', [64 94] ) ;
            [datat, geomt ] = dinplanet( datat, geomt, 'imresize', [size(vold,1) size(vold,2)] ) ;
            
            testCase.verifyEqual([geomt.IPP], [md.geom.IPP], 'RelTol', 0.001, 'AbsTol',1e-7) ;
                
            testCase.verifyEqual([geomt.IOP], [md.geom.IOP], 'RelTol', 0.001, 'AbsTol',1e-5) ;
                
            testCase.verifyEqual([geomt.PixelSpacing_HW], ...
                    [md.geom.PixelSpacing_HW], 'RelTol', 0.001, 'AbsTol',1e-7) ;
        end
        
        function testFlip(testCase)
            import matlab.unittest.fixtures.PathFixture
            DICOM_Folder = get_data_location('DICOM_20191211geometryIngenia','dir') ;
            f = testCase.applyFixture(PathFixture(DICOM_Folder));
            
            testem(testCase,  'IM_0011') % oblique
            testem(testCase,  'IM_0008') % sag
            testem(testCase,  'IM_0005') % tra 
            testem(testCase,  'IM_0002') % cor
            
            function testem(testCase, dfn)
                
                dinfo = datparse(dfn) ;
                
                [vold, md] = d2mat(dinfo, {'slice'}, 'op', 'dv') ;
                
                % round trip flipping!
                [datat, geomt ] = dinplanet( vold, md.geom, 'flip', 1 ) ;
                [datat, geomt ] = dinplanet( datat, geomt , 'flip', 2 ) ;
                [datat, geomt ] = dinplanet( datat, geomt , 'flip', 3 ) ;
                [datat, geomt ] = dinplanet( datat, geomt , 'flip', 3 ) ;
                [datat, geomt ] = dinplanet( datat, geomt , 'flip', 2 ) ;
                [datat, geomt ] = dinplanet( datat, geomt , 'flip', 1 ) ;
                
                
                testCase.verifyEqual(double(datat),double(vold), 'RelTol',0.0001, 'AbsTol',1e-7);
                
                testCase.verifyEqual([geomt.IPP], [md.geom.IPP], 'RelTol', 0.001, 'AbsTol',1e-7) ;
                
                testCase.verifyEqual([geomt.IOP], [md.geom.IOP], 'RelTol', 0.001, 'AbsTol',1e-5) ;
                
                testCase.verifyEqual([geomt.PixelSpacing_HW], ...
                    [md.geom.PixelSpacing_HW], 'RelTol', 0.001, 'AbsTol',1e-7) ;
            end
            
        end
        
        
    end
    
end