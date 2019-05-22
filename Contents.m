% FUNCTIONS
%
% Files
%   adaptive_equalizer               - Adaptive equalization at 2SPS
%   calculate_spectrum               - Calculates the spectrum of the input signal
%   CD_compensation                  - Compensate for chromatic dispersion
%   coarse_FOC                       - Coarse Frequency offset compensation
%   complex_to_real                  - Transform to real-valued signal
%   constellations                   - Constellation and bit mapping

%   delay_compensation               - Compensate for the skew between channels


%   normalize_DC_AGC                 - Remove DC and normalize power
%   performance_monitoring           - Interpret complex-valued equalizer taps


%   real_to_complex                  - Transform to complex-valued signal
%   resample_signal                  - Resample the signal at 2SPS


%   PRBS_generator                   - This function generates pseudo-random binary sequences (PRBS)








%   ideal_lpf                        - Ideal low-pass filter


%   naninterp                        - 'linear' interpolation of points with NaN


%   osnr_penalty                     - Calculates OSNR penalty


%   quantize                         - Quantization

%   raiscos_lpf                      - Root Raised-Cosine low-pass filter
%   rc_lpf                           - One-pole (RC) low-pass filter


%   sinc_lpf                         - Sinc-shaped low-pass filter

%   supergauss_lpf                   - Supergaussian low-pass filter






%   get_default_cohdsp_params        - Coherent-DSP default parameters structure


%   ber_theory_bits                  - Theoretical BER
%   phase_recovery_data_aided        - Data-aided phase recovery



%   add_enob                         - Add DAC/ADC noise
%   clip                             - ADC/DAC clipping

%   pas_decode                       - Decode signal encoded with PAS encoder
%   pas_encode                       - Encode bits using Probabilistic Amplitude Shaping (PAS)
%   phase_recovery_pilot             - BPS and ML phase recovery with pilots


%   eval_nlpn                        - Evaluate statistics of non-linear phase noise
%   harddec                          - Square QAM hard decision

%   rc_2nd_lpf                       - RC_LPF      Two-pole (RC) low-pass filter

%   find_delay_pm_2sps               - Find delay between two sequences
%   training_align                   - Align training sequence to input data




%   FOC_FD                           - FINE_FOC     Fine Frequency offset compensation

%   find_delay_pm_1sps               - Find delay between two sequences of symbols

%   get_default_cohdsp_params_fromTX - GET_DEFAULT_COHDSP_PARAMS     Coherent-DSP default parameters structure
%   phase_recovery_vv_psk            - V&V PHASE RECOVERY FOR PSK

%   decode                           - Deconding and performance evaluation
%   full_training_align              - Align training sequence to input data

%   DBP_compensation                 - Non-Linear DBP compensation

%   psk_constellations               - PSK constellation and bit mapping


%   run_tests                        - 


%   differential_decode              - Minimum distance differential deconding and BER computation
%   estimate_equivalent_snr          - Estimate the equivalent SNR from electrical spectrum
%   index_to_binary                  - Convert decimal indeces to binary




