classdef NiftiDicomTest < matlab.unittest.TestCase
    % NiftiDicomTest tests Nifti code vs Dicom
    %
    % Copyright, 2020, David Atkinson
    % D.Atkinson@ucl.ac.uk
    
    
    methods (Test)
        function test(testCase)
            import matlab.unittest.fixtures.PathFixture
            DICOM_Folder = get_data_location('DICOM_20191211geometryIngenia','dir') ;
            f = testCase.applyFixture(PathFixture(DICOM_Folder));
            NIFTI_Folder = get_data_location('NIFTI_20191211geometryIngenia','dir') ;
            f = testCase.applyFixture(PathFixture(NIFTI_Folder));
            
           testem(testCase, 'OBJECT_phantom_T2W_TSE_obl_19_1.nii', 'IM_0011')
           testem(testCase, 'OBJECT_phantom_T2W_TSE_Sag_18_1.nii', 'IM_0008')
           testem(testCase, 'OBJECT_phantom_T2W_TSE_Tra_17_1.nii', 'IM_0005', 'reverse_slices', true)
           testem(testCase, 'OBJECT_phantom_T2W_TSE_Cor_14_1.nii', 'IM_0002')
            
            function testem(testCase, nfn, dfn, varargin)
                ninfo = niftiinfo(nfn);
                dinfo = datparse(dfn) ;
            
                [voln, mn] = n2mat(ninfo, 'op','dv',varargin{:}) ;
                [vold, md] = d2mat(dinfo, {'slice'}, 'op', 'dv') ;
                
                testCase.verifyEqual(double(voln),double(vold), 'RelTol',0.0001, 'AbsTol',1e-7);
                
                testCase.verifyEqual([mn.geom.IPP], [md.geom.IPP], 'RelTol', 0.001, 'AbsTol',1e-7) ;
                
                testCase.verifyEqual([mn.geom.IOP], [md.geom.IOP], 'RelTol', 0.001, 'AbsTol',1e-5) ;
                
                testCase.verifyEqual([mn.geom.PixelSpacing_HW], ...
                    [md.geom.PixelSpacing_HW], 'RelTol', 0.001, 'AbsTol',1e-7) ;
            end
            
        end
        
        
    end
    
end