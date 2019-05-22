classdef DecodingTest < matlab.unittest.TestCase
    %DECODINGTEST Test decoding routines
    %   Use this class to test the various decoding routines
    
    % Test methods
    methods (Test)
        % Test ber_theory_bits over a square constellation
        function test_ber_theory_square(testCase)
            snr = (1:10)';
            ber = log10(ber_theory_bits(snr,4));
            ber_true = [-0.5809;...
                -0.6248;...
                -0.6733;...
                -0.7264;...
                -0.7847;...
                -0.8494;...
                -0.9228;...
                -1.0080;...
                -1.1087;...
                -1.2292];
            verifyEqual(testCase,ber,ber_true,'RelTol',1e-3);
        end
        
        % Test ber_theory_bits over a rectangular constellation
        function test_ber_theory_rect(testCase)
            snr = (5:15)';
            ber = log10(ber_theory_bits(snr,5));
            ber_true = [-0.6741;...
                -0.7103;...
                -0.7526;...
                -0.8025;...
                -0.8615;...
                -0.9315;...
                -1.0149;...
                -1.1148;...
                -1.2349;...
                -1.3799;...
                -1.5558];
            verifyEqual(testCase,ber,ber_true,'RelTol',1e-3);
        end
        
        % Decode a noiseless QAM constellation
        function test_decode_qam_noiseless(testCase)
            addpath('../../capacity-scripts')
            Parameters.Mconst = 4;
            Parameters.ber_pp = 4;
            [Parameters.C,Parameters.B] = constellations(Parameters.Mconst);
            Parameters.differential = false;
            Parameters.calculate_MI = false;
            Parameters.Ntx = 1;
            
            a = Parameters.C(end:-1:1);
            
            [Results,be,x,L] = decode(repmat(a,Parameters.ber_pp,1),a,Parameters);
            
            verifyEqual(testCase,Results.ber,0);
            verifyEqual(testCase,Results.GMI,Parameters.Mconst);
            verifyEqual(testCase,sum(be(:)),0);
            verifyTrue(testCase,isequal(x,repmat(a,Parameters.ber_pp,1)));
            verifyEqual(testCase,sum(be(:)),0);
            verifyTrue(testCase,isequal(sum(sign(L)),zeros(1,Parameters.Mconst)));
        end
        
        % Decode a noisy QAM constellation
        function test_decode_qam_noisy(testCase)
            addpath('../../capacity-scripts')
            Parameters.Mconst = 4;
            Parameters.ber_pp = 100;
            [Parameters.C,Parameters.B] = constellations(Parameters.Mconst);
            Parameters.differential = false;
            Parameters.calculate_MI = true;
            Parameters.Ntx = 1;
            SNR = 14;
            
            Es = mean(abs(Parameters.C).^2);
            a = Parameters.C(end:-1:1);
            
            s = rng;
            rng(1);
            x = repmat(a,Parameters.ber_pp,1)+randn(Parameters.ber_pp*2^Parameters.Mconst,2)*...
                [1;1j]*sqrt(Es/2)*10^(-SNR/20);
            rng(s);
            
            Results = decode(x,a,Parameters);
            [gmith,mith] = qam_gmi_mex(Parameters.C,SNR,ones(16,1)/16);
            verifyEqual(testCase,Results.GMI,gmith,'RelTol',1e-2);
            verifyEqual(testCase,Results.MI,mith,'RelTol',1e-2);
        end
        
        % Test hard decision
        function test_hard_decision(testCase)
            C = constellations(4);
            x = C + (0.1-0.1j);
            ver = isequal(C,harddec(x,4));
            verifyTrue(testCase,ver);
        end
        
        % Test evaluation of residual phase noise
        function test_eval_nlpn(testCase)
            C = constellations(4);
            a = repmat(C,100,1);
            x = a.*exp(1j*2*pi*repmat([-1;1],800,1)*0.1);
            [NCI,autocorr_phase,autocorr_abs] = eval_nlpn(x,a);
            err_phase = var(autocorr_phase,ones(1601,1));
            err_abs = var(autocorr_abs,ones(1601,1));
            
            verifyGreaterThanOrEqual(testCase,NCI,100);
            verifyEqual(testCase,err_phase,0,'AbsTol',1e-10);
            verifyEqual(testCase,err_abs,0,'AbsTol',1e-10);
        end
        
        % Differentially decode a noiseless QAM constellation
        function test_differential_decode_qam_noiseless(testCase)
            Parameters.Mconst = 4;
            Parameters.ber_pp = 4;
            [Parameters.C,Parameters.B] = constellations(Parameters.Mconst);
            Parameters.Ntx = 1;
            Parameters.differential = true;
            
            a = Parameters.C(end:-1:1);
            
            [Results,be] = differential_decode(repmat(a,Parameters.ber_pp,1),a,Parameters);
            
            verifyEqual(testCase,Results.ber,0);
            verifyEqual(testCase,sum(be(:)),0);
        end
        
        % Differentially decode a noiseless QAM constellation
        function test_differential_decode_qam_noisy(testCase)
            Parameters.Mconst = 4;
            Parameters.ber_pp = 4;
            [Parameters.C,Parameters.B] = constellations(Parameters.Mconst);
            Parameters.Ntx = 1;
            Parameters.differential = true;            
            SNR = 12;
            Es = mean(abs(Parameters.C).^2);
            a = Parameters.C(end:-1:1);
            
            s = rng;
            rng(1);
            x = repmat(a,Parameters.ber_pp,1)+randn(Parameters.ber_pp*2^Parameters.Mconst,2)*...
                [1;1j]*sqrt(Es/2)*10^(-SNR/20);
            rng(s);
            
            Results = differential_decode(x,a,Parameters);
            
            verifyEqual(testCase,Results.ber,0.0391,'AbsTol',0.01);
        end

    end
end