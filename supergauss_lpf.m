function H = supergauss_lpf(L,f0,n,f3dB,os)
% Supergaussian low-pass filter
% From Andrea Arduino
% INPUTS
% L         :=  Length of the filter (even, usually a power of 2)
% f0        :=  Central frequency [-0.5 ... +0.5]*os
% n         :=  Filter order
% f3dB      :=  Filter 3-dB frequency [0 ... +0.5]*os
% os        := 	Oversampling factor
% OUTPUTS
% H         :=  Filter in frequency domain (fft convention)
f_ax = (-L/2:L/2-1)'/L*os;                                                 % frequency axis
sig = f3dB/(log(2)^(1/2/n));                                               % variance of Gaussian
H = exp(-1/2*((f_ax-f0)/sig).^(2*n));                                      % Generate Gaussian
H = ifftshift(H);                                                          % shift to follow fft's convention