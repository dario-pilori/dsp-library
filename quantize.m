function y = quantize(x,Nbit)
%QUANTIZE   Quantization
% This function performs the quantization of the input signal
% Inputs:
% x     := input signal
% Nbit  := number of bit of the quantizer
% Outputs:
% y     := quantized output signal

% Dec. 2015 Dario Pilori <dario.pilori@polito.it>

%% Init
% If Nbit is too high, give up
if Nbit>=64
    y = x;
    return
end
% If complex, add another dimension
if ~isreal(x)
    z = complex_to_real(x);
    cmplx = true;
else
    z = x;
    cmplx = false;
end

%% Quantize
% Calculate range
max_r = max(max(z));
min_r = min(min(z));
y = round((z+min_r)/(max_r-min_r)*(2^Nbit-1)); % quantize to 2^Nbit
y = y/(2^Nbit-1)*(max_r-min_r)-min_r;          % rescale

%% End
% if input was complex, return to complex
if cmplx
    y = real_to_complex(y);
end
end