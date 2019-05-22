classdef PhaseRecoveryTest < matlab.unittest.TestCase
    %PHASERECOVERYTEST Test equalizer
    %   Use this class to test the adaptive equalizer routines

    % Test methods
    methods (Test)
        % Verify blind BPS+ML
        function test_blind_bps_ml(testCase)
            % Parameters
            Ns = 1e3;
            Parameters.cpe_memory = 32;
            Parameters.bps_anglestep = 5*(pi/180);
            Parameters.C = constellations(4);
            Parameters.cpe_pilotlen = 0;
            Parameters.cpe_blocksize = 2^10;
            Parameters.cpe_sps = 1;
            Parameters.cpe_ml_dec = false;
            Parameters.Ntx = 1;
            Parameters.cpe_sym = pi/2;            
                        
            % Generate signal
            a = Parameters.C(randi(16,Ns,1));
            
            % Generate noise
            s = rng;
            rng(1);
            n = randn(Ns,2)*[1;1j]*sqrt(1e-2);
            phn = cumsum(randn(Ns,1)*sqrt(2*pi*200e3/32e9));
            rng(s);
            
            % Apply noise
            x = a.*exp(1j*phn) + n;
            
            % Apply phase recovery
            [y,ph] = phase_recovery_pilot(x,a,Parameters);
            
            % Verify result
            verifyEqual(testCase,harddec(y,4),a);
            verifyEqual(testCase,ph,phn,'AbsTol',5*pi/180);
        end
        
        % Verify pilot/based BPS+ML
        function test_pilot_bps_ml(testCase)
            % Parameters
            Ns = 2^12;
            Parameters.cpe_memory = 32;
            Parameters.bps_anglestep = 5*(pi/180);
            Parameters.C = constellations(4);
            Parameters.cpe_pilotlen = 4;
            Parameters.cpe_blocksize = 2^8;
            Parameters.cpe_sps = 1;
            Parameters.cpe_ml_dec = false;
            Parameters.Ntx = 1;
            Parameters.cpe_sym = pi/2;
            
            % Generate signal
            a = Parameters.C(randi(16,Ns,1));
            
            % Generate noise
            s = rng;
            rng(1);
            n = randn(Ns,2)*[1;1j]*sqrt(1e-2);
            phn = cumsum(randn(Ns,1)*sqrt(2*pi*200e3/32e9));
            rng(s);
            
            % Apply noise
            x = a.*exp(1j*phn) + n;
            
            % Apply phase recovery
            [y,ph] = phase_recovery_pilot(x,a,Parameters);
            
            % Verify result
            verifyEqual(testCase,harddec(y,4),a);
            verifyEqual(testCase,ph,phn,'AbsTol',5*pi/180);
        end
        
        % Verify pilot/based V&V
        function test_pilot_vv(testCase)
            % Parameters
            Ns = 2^12;
            Parameters.cpe_memory = 32;
            Parameters.C = constellations(2);
            Parameters.cpe_pilotlen = 4;
            Parameters.cpe_blocksize = 2^8;
            Parameters.cpe_sps = 1;
            Parameters.Ntx = 1;
            Parameters.cpe_sym = pi/2;
            
            % Generate signal
            a = Parameters.C(randi(4,Ns,1));
            
            % Generate noise
            s = rng;
            rng(1);
            n = randn(Ns,2)*[1;1j]*sqrt(1e-2);
            phn = cumsum(randn(Ns,1)*sqrt(2*pi*500e3/32e9));
            rng(s);
            
            % Apply noise
            x = a.*exp(1j*phn) + n;
            
            % Apply phase recovery
            y = phase_recovery_vv_psk(x,a,Parameters);
            
            % Verify result
            verifyEqual(testCase,harddec(y,2),a);
        end
        
        % Verify blind V&V
        function test_blind_vv(testCase)
            % Parameters
            Ns = 2^12;
            Parameters.cpe_memory = 32;
            Parameters.C = constellations(2);
            Parameters.cpe_pilotlen = 0;
            Parameters.cpe_blocksize = 2^8;
            Parameters.cpe_sps = 1;
            Parameters.Ntx = 1;
            Parameters.cpe_sym = pi/2;
            
            % Generate signal
            a = Parameters.C(randi(4,Ns,1));
            
            % Generate noise
            s = rng;
            rng(1);
            n = randn(Ns,2)*[1;1j]*sqrt(1e-2);
            phn = cumsum(randn(Ns,1)*sqrt(2*pi*200e3/32e9));
            rng(s);
            
            % Apply noise
            x = a.*exp(1j*phn) + n;
            
            % Apply phase recovery
            y = phase_recovery_vv_psk(x,a,Parameters);
            
            % Rotate phase
            ph_rot = (angle(y(1))-angle(a(1)));
            y = y*exp(-1j*ph_rot);
            
            % Verify result
            verifyEqual(testCase,harddec(y,2),a);
        end
        
    end
end