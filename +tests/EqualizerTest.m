classdef EqualizerTest < matlab.unittest.TestCase
    %EQUALIZERTEST Test equalizer
    %   Use this class to test the adaptive equalizer routines

    % Test methods
    methods (Test)
        % Verify complex LMS equalizer
        function test_lms_complex(testCase)
            % Parameters
            Ns = 1e3;
            Parameters.mus = 5e-2;
            Parameters.Ntaps = 4;
            Parameters.C = constellations(2);
            Parameters.Mlms = 1;
            Parameters.lambda = 0;
            Parameters.Ks = Inf;
            Parameters.eqmode = 'complex';
            Parameters.methods = {'lms'};
            
            % Generate signal
            h = [1;-0.2j];
            a = Parameters.C(randi(4,Ns,1));
            
            % Generate noise
            s = rng;
            rng(1);
            n = randn(Ns*2,2)*[1;1j]*sqrt(1e-3);
            rng(s);
            
            % Equalize
            y = adaptive_equalizer(filter(h,1,upsample(a,2))+n,a,Parameters);
            y = y(1:2:end);
            
            % Decode and evaluate symbol errors
            d = sign(real(y))+1j*sign(imag(y));
            se = sum(d(100:end)~=a(101:997));
            verifyEqual(testCase,se,0);
        end
        
        % Verify complex CMA equalizer
        function test_cma_complex(testCase)
            % Parameters
            Ns = 1e3;
            Parameters.mus = 5e-2;
            Parameters.Ntaps = 4;
            Parameters.C = constellations(2);
            Parameters.Mlms = 1;
            Parameters.lambda = 0;
            Parameters.Ks = Inf;
            Parameters.eqmode = 'complex';
            Parameters.methods = {'cma'};
            
            % Generate signal
            h = [1;-0.2j];
            a = Parameters.C(randi(4,Ns,1));
            
            % Generate noise
            s = rng;
            rng(1);
            n = randn(Ns*2,2)*[1;1j]*sqrt(1e-3);
            rng(s);
            
            % Equalize
            y = adaptive_equalizer(filter(h,1,upsample(a,2))+n,a,Parameters);
            y = y(1:2:end);
            
            % Decode and evaluate symbol errors
            d = sign(real(y))+1j*sign(imag(y));
            se = sum(conj(d(100:end))~=a(101:997));
            verifyEqual(testCase,se,0);
        end
        
        % Verify real LMS equalizer
        function test_lms_real(testCase)
            % Parameters
            Ns = 1e3;
            Parameters.mus = 5e-2;
            Parameters.Ntaps = 4;
            Parameters.C = constellations(2);
            Parameters.Mlms = 1;
            Parameters.lambda = 0;
            Parameters.Ks = Inf;
            Parameters.eqmode = 'real';
            Parameters.methods = {'lms'};
            
            % Generate signal
            h = [1;-0.2j];
            a = Parameters.C(randi(4,Ns,1));
            
            % Generate noise
            s = rng;
            rng(1);
            n = randn(Ns*2,2)*[1;1j]*sqrt(1e-3);
            rng(s);
            
            % Equalize
            y = adaptive_equalizer(filter(h,1,upsample(a,2))+n,a,Parameters);
            y = y(1:2:end);
            
            % Decode and evaluate symbol errors
            d = sign(real(y))+1j*sign(imag(y));
            se = sum(d(100:end)~=a(101:997));
            verifyEqual(testCase,se,0);
        end
        
    end
end