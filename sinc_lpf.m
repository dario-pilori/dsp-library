function H = sinc_lpf(L, f3dB, os)
% Sinc-shaped low-pass filter
% From Andrea Arduino
% INPUTS
% L         :=  Length of the filter (even)
% f3dB      :=  3-dB frequency [0 ... +0.5]*os
% os        := 	Oversampling factor
% OUTPUTS
% H         :=  LPF in frequency domain
if f3dB>=Inf
    H = ones(L,1);
    return;
end
f_ax = (-L/2:L/2-1)'/L*os;
tau = f3dB/0.442946;
H = sin(pi*f_ax/tau)./(pi*f_ax/tau);
H((f_ax/tau)==0) = 1;
H = ifftshift(H);