classdef GeomTest < matlab.unittest.TestCase
    % GeomTest tests geometry related code
    %
    % Copyright, 2019, David Atkinson
    % D.Atkinson@ucl.ac.uk
    
    
    methods (Test)
        function testSetPlane(testCase)
            actSolution = set_plane('pvv',[0 0 0],[1 0 0],[0 1 0]);
            expNormal = [0 0 1];
            testCase.verifyEqual(actSolution.normal,expNormal);
        end
        
        function testOriInfo(testCase)
            
            import matlab.unittest.fixtures.PathFixture
            dataFolder = get_data_location('DICOM_20191211geometryIngenia','dir') ;
            f = testCase.applyFixture(PathFixture(dataFolder));
            
            
            % Confirmed against OsiriX
            lfn = 'IM_0008' ; ori_str = 'SAG' ;
            dinfo = datparse(lfn) ;
            [v,m] = d2mat(dinfo,{'slice'}) ;
            actSolution = ori_info(m.geom(1).IOP)
            testCase.verifyEqual(actSolution.plane_str,ori_str);
            
            lfn = 'IM_0002' ; ori_str = 'COR' ; east_str = 'LF' ;
            dinfo = datparse(lfn) ;
            [v,m] = d2mat(dinfo,{'slice'}) ;
            actSolution = ori_info(m.geom(1).IOP)
            testCase.verifyEqual(actSolution.plane_str,ori_str);
            testCase.verifyEqual(actSolution.east_str,east_str);
            
            lfn = 'IM_0005' ; ori_str = 'TRA' ; west_str = 'R' ;
            dinfo = datparse(lfn) ;
            [v,m] = d2mat(dinfo,{'slice'}) ;
            actSolution = ori_info(m.geom(1).IOP, 'label_tolerance', 0.03)
            testCase.verifyEqual(actSolution.plane_str,ori_str);
            testCase.verifyEqual(actSolution.west_str,west_str);
            
            lfn = 'IM_0011' ; east_str = 'PHR' ;  north_str = 'HLA' ;
            dinfo = datparse(lfn) ;
            [v,m] = d2mat(dinfo,{'slice'}) ;
            actSolution = ori_info(m.geom(1).IOP)
            testCase.verifyEqual(actSolution.east_str,east_str);
            testCase.verifyEqual(actSolution.north_str,north_str);
            
        end
    end
    
end