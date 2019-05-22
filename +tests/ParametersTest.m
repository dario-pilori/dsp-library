classdef ParametersTest < matlab.unittest.TestCase
    %PARAMETERSTEST Test parameters
    %   This class tests the parameters
        
    methods (Test)
        function test_default_params(testCase)
            Parameters = get_default_cohdsp_params(4);
            verifyNotEmpty(testCase,Parameters);
        end
        
        function test_default_params_fromTX(testCase)
            TX.Ntx = 2;
            TX.Mconst = 4;
            [TX.C,TX.B] = constellations(4);
            TX.optimized_const = false;
            TX.symRate = 32e9;
            
            Parameters = get_default_cohdsp_params_fromTX(TX);
            verifyNotEmpty(testCase,Parameters);
        end
    end
end

