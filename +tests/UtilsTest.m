classdef UtilsTest < matlab.unittest.TestCase
    %UTILSTEST Test various small utilities
    
    % Test methods
    methods (Test)
        % Test function naninterp
        function test_naninterp(testCase)
            x = [1;2;NaN;4];
            y = naninterp(x);
            ver = isequal(y,(1:4)');
            verifyTrue(testCase,ver);
        end
        
        % Test function osnr_penalty
        function test_osnr_penalty(testCase)
            rosnr = osnr_penalty((1:10)',linspace(0.3,0.1,10)',0.2);
            verifyEqual(testCase,rosnr,5.5,'AbsTol',1e-10);
        end
    end
end