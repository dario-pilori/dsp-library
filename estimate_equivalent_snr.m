function [SNR_mmse_le,SNR_mmse_dfe] = estimate_equivalent_snr(x,Nfft,...
    OS_filt,OS,ro,varargin)
%ESTIMATE_EQUIVALENT_SNR    Estimate the equivalent SNR from electrical spectrum
% This function estimates the equivalent SNR after ideal Minimum-Mean
% Square Error (MMSE) Linear Equalization (LE) and Decision-Feedback
% Equalization (DFE). 'Ideal' means with infinite number of taps a no error
% propagation (for DFE only). Note that ideal MMSE-DFE is
% capacity-achieving, i.e. the channel capacity is exactly
% log2(1+SNR_MMSE_DFE).
% For details, see https://doi.org/10.1002/0471439002.ch2
%
%   INPUTS:
%   x              :=  Received electrical signal, column vector(s) (at least 2SPS)
%   Nfft           :=  FFT size for the periodogram
%   OS_filt        :=  Extra oversampling factor (usually 4)
%   OS             :=  Input oversampling factor (usually 2)
%   ro             :=  RRC roll-off factor
%   start_filt     :=  Tap where the optical filter starts (optional)
%   end_filt       :=  Tap where the optical filter ends (optional)
%
%   OUTPUTS
%   SNR_mmse_le    :=  Equivalent SNR after ideal MMSE-LE (dB)
%   SNR_mmse_dfe   :=  Equivalent SNR after ideal MMSE-DFE (dB)

%   Mar. 2019 - Dario Pilori <dario.pilori@polito.it>

%% Verify and get parameters
validateattributes(x,{'numeric'},{'2d'},'','x',1);
validateattributes(Nfft,{'numeric'},{'scalar','integer','positive'},'','Nfft',2);
validateattributes(OS_filt,{'numeric'},{'scalar','integer','>',1},'','OS_filt',3);
validateattributes(OS,{'numeric'},{'scalar','integer','>=',2},'','OS',4);
validateattributes(ro,{'numeric'},{'scalar','nonnegative','<=',1},'','ro',5);

%% Estimate signal and noise PSDs
N_ch = size(x,2);
pxx = pwelch(x,Nfft);
start_nois = ceil(length(pxx)/2/OS*(1+ro))+round(Nfft/100);
end_nois = floor(-start_nois+length(pxx))-round(Nfft/100);

%% Remove optical filter
if nargin>6
    start_filt = varargin{1};
    validateattributes(start_filt,{'numeric'},...
        {'scalar','integer','nonnegative'},'','start_filt',6);
    end_filt = varargin{2};
    validateattributes(end_filt,{'numeric'},...
        {'scalar','integer','>=',start_filt},'','end_filt',7);
else
    start_filt = end_nois;
    end_filt = end_nois;
end

%% Calculate PSDs
Pnois = mean([pxx(start_nois:start_filt,:);...
    pxx(end_filt:end_nois,:)]);
Ps = max(eps,pxx - Pnois);

%% Upsample and sum repetitions
Ps = [Ps(1:Nfft/2,:);...
    zeros(Nfft*(OS_filt-1),N_ch);...
    Ps(Nfft/2+1:end,:)];
num = zeros(OS_filt*Nfft,N_ch);
mu = (-2:2)';
for i = 1:length(mu)
    num = num + circshift(Ps,[mu(i)*Nfft/OS,0]);
end

%% Calculate equivalent SNRs
snr_per_freq = num ./ Pnois;
snr_per_freq = [snr_per_freq(1:Nfft/4,:);...
    snr_per_freq(end-Nfft/4+1:end,:)];
    
SNR_mmse_le = -10*log10(2/Nfft*sum(1./(snr_per_freq+1))).'; 
SNR_mmse_dfe = 10*log10(exp(1))*2/Nfft*sum(log(snr_per_freq+1)).';
end