function H = ideal_lpf(L,f_min,f_max,os)
%IDEAL_LPF   Ideal low-pass filter
% This function generates an ideal low-pass filter.
% INPUTS
% L         :=  Length of the filter (even)
% f_min     :=  Minimum frequency [-0.5 ... +0.5]*os
% f_max     :=  Maximum frequency [-0.5 ... +0.5]*os
% os        := 	Oversampling factor
% OUTPUTS
% H         :=  LPF in frequency domain

% Initially written by Andrea Arduino, modified by Dario Pilori

f_ax = (-L/2:L/2-1)'/L*os;
H = ones(size(f_ax));
H((f_ax<f_min-eps)|(f_ax>f_max+eps)) = 0;
H = ifftshift(H);