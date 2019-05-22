function y = clip(x,CLev)
%CLIP  ADC/DAC clipping
% This function clips the signal x between +-CLev*std(x)
% Inputs:
% x     := input signal
% CLev  := clipping level (positive)
% Outputs:
% y     := clipped output signal

% Dec. 2016 Dario Pilori <dario.pilori@polito.it>

%% Init
% If CLev is too high, give up
if CLev>=1e3
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

%% Clip
mu = CLev*std(z);                                                          % clipping threshold
y = NaN(size(z));                                                          % initialize output
for n = 1:size(z,2)
    zn = z(:,n);                                                           % n-th channel
    zn(abs(zn)>mu(n)) = mu(n)*sign(zn(abs(zn)>mu(n)));                     % clip
    y(:,n) = zn;                                                           % and go back
end

%% End
% if input was complex, return to complex
if cmplx
    y = real_to_complex(y);
end
end