classdef FrontEndCorrTest < matlab.unittest.TestCase
    %FRONTENDCORRTEST Test front-end corrections
    %   Use this class to test the front-end corrections, applied before
    %   adaptive equalization
    
    methods (Test)
        % Test spectrum (only sum)
        function test_calculate_spectrum(testCase)
            x = (1:1024)/1024;
            s = calculate_spectrum(x,1);
            verifyEqual(testCase,sum(s),24.2814,'RelTol',1e-4);
        end
        
        % Apply and compensate for CD
        function test_cd_compensation(testCase)
            s = rng;
            rng(1);
            x = randi(2,2^10,1)*2-3;
            rng(s);
            Parameters.CD = -100e-3;
            Parameters.Rs = 32e9;
            Parameters.f0 = 193.4e12;
            
            y = CD_compensation(x,Parameters);
            Parameters.CD = 100e-3;
            z = real(CD_compensation(y,Parameters));
            
            err = std(z-x(15:end-14));
            verifyLessThanOrEqual(testCase,err,1e-3);
        end
        
        % Digital Back Propagation
        function test_dbp(testCase)
            s = rng;
            rng(1);
            x = randi(2,2^11,1)*2-3;
            rng(s);
            Parameters.Rs = 32e9;
            Parameters.f0 = 193.4e12;
            Parameters.Lspan = 100e3;
            Parameters.adB = 0.2e-3;
            Parameters.gam = 1.3e-3;
            Parameters.b2 = -21.27e-27;
            Parameters.DBPsteps = 2;
            Parameters.csi = 1;
            Parameters.Pch = 0;
            Parameters.Nspan = 5;
            
            y = DBP_compensation(x,Parameters,2,false);

            Parameters.gam = -1.3e-3;
            Parameters.b2 = +21.27e-27;
            z = real(DBP_compensation(y,Parameters,2,false));
            
            err = std(z(100:end-100)-x(100:end-100));
            verifyLessThanOrEqual(testCase,err,5e-2);
        end
        
        % Test complex to real
        function test_complex_to_real(testCase)
            x = complex_to_real([1+1j,2-1j;-1-1j, 3+4j]);
            x_true = [1,1,2,-1;-1,-1,3,4];
            verifyEqual(testCase,x,x_true);
        end
        
        % Test real to complex
        function test_real_to_complex(testCase)
            x = real_to_complex([1,1,2,-1;-1,-1,3,4]);
            x_true = [1+1j,2-1j;-1-1j, 3+4j];
            verifyEqual(testCase,x,x_true);
        end
        
        % Delay a sine wave
        function test_delay_compensation(testCase)
            t = (0:0.01:5)';
            x = sin(pi*t);
            p.fADC = 1/0.01;
            p.skew = 0.5;
            y = delay_compensation(x,p);
            err = var(y(1:401)-cos(pi*t(1:401)));
            verifyEqual(testCase,err,0,'AbsTol',1e-10);
        end
        
        % Verify find_delay_pm_2sps
        function test_delay_2sps_even(testCase)
            a = randi(2,16,1)*2-3;
            x = circshift(a,-4);
            x2 = reshape(repmat(reshape(x,[],1),1,2).',[],size(x,2));
            t0 = find_delay_pm_2sps(x2,a);
            verifyEqual(testCase,t0,8);
        end
        
        function test_delay_2sps_odd(testCase)
            a = randi(2,19,1)*2-3;
            x = circshift(a,-7);
            x2 = reshape(repmat(reshape(x,[],1),1,2).',[],size(x,2));
            t0 = find_delay_pm_2sps(x2,a);
            verifyEqual(testCase,t0,14);
        end
        
        % Verify training_align
        function test_training_align(testCase)
            a = randi(2,16,1)*2-3;
            x = circshift(a,-7);
            x2 = reshape(repmat(reshape(x,[],1),1,2).',[],size(x,2));
            out = training_align(x2,a,false,false,false);
            verifyTrue(testCase,isequal(out,circshift(x2,14)));
        end
        
        % Verify training_align
        function test_full_training_align(testCase)
            a = randi(2,16,1)*2-3;
            x = circshift(a,-7);
            out = full_training_align(x,a);
            verifyTrue(testCase,isequal(out,circshift(a,-7)));
        end
        
        % Test coarse FOC
        function test_coarse_foc(testCase)
            N = 2^12;
            t = (-N/2:N/2-1)'/N*100;
            x = sinc(t);
            y = (x+0.05).*exp(1j*2*pi*(0:length(x)-1)/8).';
            p.kDC = 1e4;
            z = coarse_FOC(y,p);
            verifyEqual(testCase,var(z-x),0,'AbsTol',1e-10);
        end
                
    end
    
end