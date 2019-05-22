function [pxx,f] = calculate_spectrum(x,fs)
%CALCULATE_SPECTRUM     Calculates the spectrum of the input signal
%   This function calculates the spectrum of the signals in x (column vectors),
%   using a periodogram
%
%   It requires pwelch, introduced before MATLAB R2006a in the
%   Signal Processing Toolbox
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process [N x M]
%   fs          :=  Sampling frequency of the signal (Hz) [1 x 1]
%   OUTPUTS:
%   pxx         :=  Power spectral density, zero at center [Lwch x M]
%   fs          :=  Frequency axis (Hz) [Lwch x 1]
%
%   For info see: PWELCH

%   Mar. 2016 - Dario Pilori <dario.pilori@polito.it>

%% Check input parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(fs,{'numeric'},{'scalar','positive'},'','fs',2);

%% Calculate PSD
Lwch = 512;
pxx = fftshift(pwelch(x+1j*eps,Lwch));
f = (-Lwch/2:Lwch/2-1)'/Lwch*fs;