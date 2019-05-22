classdef FilterTest < matlab.unittest.TestCase
    %FILTERTEST Test filters
    %   Use this class to test the various filters
    
    % Test methods
    methods (Test)
        % Ideal low-pass filter
        function test_ideal_lpf(testCase)
            H = ideal_lpf(16,1/8,1/4,1);
            H_true = [zeros(2,1);ones(3,1);zeros(11,1)];
            ver = isequal(H,H_true);
            verifyTrue(testCase,ver);
        end
        
        % Raised cosine filter
        function test_raiscos_lpf(testCase)
            H_true = [1.0000;    0.9952;    0.9808;    0.9569;    0.9239;...
                0.8819;    0.8315;    0.7730;    0.7071;    0.7730;...
                0.8315;    0.8819;    0.9239;    0.9569;    0.9808;...
                0.9952];
            H = raiscos_lpf(16,1,1);
            er = var(H-H_true);
            
            verifyEqual(testCase,er,0,'AbsTol',1e-6);
        end
        
        % Two-pole RC filter
        function test_rc_2nd_lpf(testCase)            
            H = rc_2nd_lpf(8, 0.2, 1);
            H_true = [1.0000 + 0.0000i;  -0.3185 - 0.4074i;...
                -0.1376 + 0.0652i;  -0.0186 + 0.0454i;   0.0009 - 0.0190i;...
                -0.0186 - 0.0454i;  -0.1376 - 0.0652i;  -0.3185 + 0.4074i];
            er = var(H-H_true);
            verifyEqual(testCase,er,0,'AbsTol',1e-6);
        end
        
        % Single-pole RC filter
        function test_rc_lpf(testCase)
            H_true = [1.0000 + 0.0000i;   0.7191 - 0.4494i;   0.3902 - 0.4878i;...
                0.2215 - 0.4152i;   0.1379 + 0.3448i;   0.2215 + 0.4152i;....
                0.3902 + 0.4878i;   0.7191 + 0.4494i];
            H = rc_lpf(8, 0.2, 1);
            er = var(H-H_true);
            verifyEqual(testCase,er,0,'AbsTol',1e-6);
        end
        
        % Sinc-shaped LPF
        function test_sinc_lpf(testCase)
            H = sinc_lpf(8, 0.25, 1);
            H_true = [1.0000;    0.9212;    0.7071;    0.4166;    0.1261;...
                0.4166;    0.7071;    0.9212];
            er = var(H-H_true);
            verifyEqual(testCase,er,0,'AbsTol',1e-6);
        end        
        
        % Super-Gaussian LPF
        function test_superg_lpf(testCase)
            H = supergauss_lpf(8,0,4,0.2,1);
            H_true = [1.0000;    0.9920;    0.1267;    0.0000;    0.0000;...
                0.0000;    0.1267;    0.9920];
            er = var(H-H_true);
            verifyEqual(testCase,er,0,'AbsTol',1e-6);
        end
    end
end
