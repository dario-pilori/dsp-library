classdef QuantizeTest < matlab.unittest.TestCase
    %QUANTIZETEST Test simulation of DAC/ADC
    %   Use this class to test the simulations of DAC, ADC and quantization
    %   in general
    
    % Test methods
    methods (Test)
        % Verify quantization
        function test_quantize(testCase)
            N_max = 10; % Maximum number of bits
            x = (0:2^N_max-1)'+1j*(2^N_max-1:-1:0)';
            for n = N_max:-1:1
                y = quantize(x,n);
                verifyEqual(testCase,numel(unique(y)),2^n);
            end
        end
        
        % Verify ENOB-noise
        function test_enob(testCase)
            x = repmat([0;1],2^12,1);
            y = add_enob(x,[6;6],[0;1],1);
            Pnois_est = var(y-x);
            Pnois_true = 2^(-12)/12;
            verifyEqual(testCase,Pnois_est,Pnois_true,'RelTol',5e-2);
        end
        
        % Verify clipping
        function test_clip(testCase)
            x = (1:10)';
            y = clip(x,5/std(x));
            eq = isequal(y,[1:4 5*ones(1,6)]');
            verifyTrue(testCase,eq);
        end
    end
end