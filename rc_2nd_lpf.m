function H = rc_2nd_lpf(L, f3dB, os)
%RC_LPF      Two-pole (RC) low-pass filter
% This function returns the FFT of a one-pole RC low-pass filter.
% INPUTS
% L         :=  Length of the filter (even)
% f3dB      :=  3-dB frequency [0 ... +0.5]*os
% os        := 	Oversampling factor
% OUTPUTS
% H         :=  LPF in frequency domain

% Initially written by Andrea Arduino, modified by Dario Pilori

f_ax = (-L/2:L/2-1)'/L*os;
tau = 1/(2*pi*f3dB);
H = ifftshift((1/tau)*(tau./((1+2j*pi*f_ax*tau).^4)));