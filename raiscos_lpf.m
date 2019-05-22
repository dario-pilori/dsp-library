function H = raiscos_lpf(L,ro,os)
% Root Raised-Cosine low-pass filter
% INPUTS  %
% L         :=  Length of the filter (even)
% ro        :=  Roll-Off factor [0...1]
% os        := 	Oversampling factor
% OUTPUTS %
% H         :=  LPF in frequency domain
f_ax = (-L/2:L/2-1)'/L*os;
H = zeros(L,1);
H(abs(f_ax)<=(1+ro)/2&abs(f_ax)>=(1-ro)/2) =...
    0.5*(1+cos(pi/ro*(abs( f_ax(abs(f_ax)<=(1+ro)/2&abs(f_ax)>=(1-ro)/2)) -(1-ro)/2)));
H(abs(f_ax)<=(1-ro)/2) = 1;
H = sqrt(ifftshift(H));