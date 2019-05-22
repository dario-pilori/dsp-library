function y = add_enob(x,Enob,Enob_f,fs)
%ADD_ENOB   Add DAC/ADC noise
% This function adds DAC/ADC noise to the input signal x based on the Enob
% (equivalent-number-of-bits) measurements;
% Inputs:
% x     := input signal
% Enob  := enob measurements at different frequencies (bit)
% Enob_f:= frequencies of enob measurements (Hz)
% fs    := signal sampling frequency (Hz)
% Outputs:
% y     := output signal

% Dec. 2016 Dario Pilori <dario.pilori@polito.it>

%% Init
% If complex, add another dimension
if ~isreal(x)
    z = complex_to_real(x);
    cmplx = true;
else
    z = x;
    cmplx = false;
end

%% Add enob
Ns = size(z,1);                                                            % signal length
Delta = (max(z)-min(z))./(2.^interp1(Enob_f,Enob,(0:ceil(Ns/2))'/Ns*fs));  % quantization bin size
if any(isnan(Delta))
    error('Enob measurements not sufficient to generate noise profile');
end
Pnois = Delta.^2/12;                                                       % power of uniformly-distributed noise
if mod(Ns,2)
    Pnois = [Pnois;flipud(Pnois(2:end-2,:))];                              % two-sided SNR (dB)
else
    Pnois = [Pnois;flipud(Pnois(2:end-1,:))];                              % two-sided SNR (dB)
end
n =  ifft(fft(randn(Ns,size(z,2))).*sqrt(Pnois));                          % generate colored noise
y = z + n;                                                                 % add noise

%% End
% if input was complex, return to complex
if cmplx
    y = real_to_complex(y);
end
end