classdef MonitoringTest < matlab.unittest.TestCase
    %MONITORING Test monitoring routines
    %   Use this class to test the various performance monitoring routines
    
    % Test methods
    methods (Test)
        function test_equivalent_snr(testCase)
            N = 2^12;
            t = (-N/2:N/2-1)'/N*512;
            x = sinc(t);

            s = rng;
            rng(1);
            y = x + randn(N,1)*0.01 + 1j*eps;
            rng(s);            
            
            snr = estimate_equivalent_snr(y,256,...
                2,8,0.01);

            verifyEqual(testCase,snr,22.8629,'AbsTol',1e-3);
        end
    end
end