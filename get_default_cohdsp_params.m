function Parameters = get_default_cohdsp_params(M,varargin)
%GET_DEFAULT_COHDSP_PARAMS     Coherent-DSP default parameters structure
%   This function returns the default parameters for coherent-DSP in the
%   Parameters structure
%
%   Always run this function FIRST, then modify each parameter as needed in
%   the main script.
%
%   ALL numerical values are in SI units (e.g. seconds, meters,
%   Hertz, ...) and not their multiples
%
%   INPUTS:
%   M   :=  Number of bits/symbol in constellation
%
%   OUTPUTS:
%   p   :=  Parameters structure

%   Mar. 2016 - Dario Pilori <dario.pilori@polito.it>

%% Initialize
Parameters = struct;                                                       % Initialize as empty struct
validateattributes(M,{'numeric'},{'scalar','positive','integer'},'Parameters','Mconst',2);
if nargin>=2&&varargin{1}
    optimized_constellations = true;
else
    optimized_constellations = false;
end

%% Constellation and bit mapping
Parameters.Ntx = 2;                                                        % number of transmitters
Parameters.Mconst = M;                                                     % bits/symbol of constellation (M=6 for 64-QAM)
[Parameters.C,Parameters.B] = constellations(Parameters.Mconst,optimized_constellations);% get constellation and bit mapping

%% Acquisition parameters
Parameters.fADC = 64e9;                                                    % ADC sampling frequency (Hz)
Parameters.Rs   = 32e9;                                                    % Signal symbol rate (Baud)

%% DC, gain and skew adjustments
Parameters.kDC  = 1e4;                                                     % memory for DC removal
Parameters.kAGC = 1e4;                                                     % memory for Automatic Gain Control
Parameters.skew = [0;0;0;0]*1e-12;                                         % skew to compensate for (s)

%% Chromatic dispersion removal
Parameters.f0   = 193.4e12;                                                % signal central frequency (Hz)
Parameters.CD   = 0;                                                       % cumulated dispersion (s/m) (1 ps/nm = 1e-3 s/m)

%% Timing recovery
% not yet implemented

%% Adaptive equalization
Parameters.eq.Ntaps = 16;                                                  % adaptive equalizer taps (T/2-spaced)
Parameters.eq.eqmode = 'complex';                                          % complex-valued or real-valued equalization
Parameters.eq.methods = {'lms';'lms'};                                     % equalizer methods (LMS, MMA-DA, MMA-DD, CMA), lowercase cell array [N x 1]
Parameters.eq.mus   = [8e-4;5e-4];                                         % adaptation coefficents [N x 1]
Parameters.eq.Ks    = [1e4;Inf];                                           % number of symbols to run the first N-1 equalizers [N x 1]
Parameters.eq.Mlms  = 10;                                                  % memory of data-aided LMS phase recovery
Parameters.eq.lambda = 0;                                                  % CMA leakage factor (between 0 and 1, 0 is regular CMA)
Parameters.eq.C = Parameters.C;                                            % constellation for equalizer

%% Phase recovery
Parameters.cpe_memory = 64;                                                % memory in samples of phase recovery
Parameters.cpe_sps = 1;                                                    % samples per symbol of phase recovery (1 or 2)
Parameters.cpe_conj = true(Parameters.Ntx,1);                              % which TX sequences are conjugated
Parameters.bps_anglestep = 5*(pi/180);                                     % BPS angle step (radians)
Parameters.cpe_pilotlen = 10;                                              % length of pilot sequence (samples)
Parameters.cpe_blocksize = 2^10;                                           % length of block for pilot sequences (samples)
Parameters.cpe_ml_dec = optimized_constellations;                          % use ML decision in case of geometrical shaping
Parameters.cpe_sym = pi/2;                                                 % constellation rotation symmetry (pi/2 for QAM)

%% Adaptive post-equalizer
% to be implemented

%% Alignment and BER counting
Parameters.ber_pp = 6;                                                     % number of pattern periods to calculate the BER
Parameters.differential = false;                                           % differential (quadrant-based) decoding to avoid cycle slips
Parameters.calculate_MI = false;                                           % do not calculate by default Mutual Information
