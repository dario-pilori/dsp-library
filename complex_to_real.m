function y = complex_to_real(x)
%COMPLEX_TO_REAL     Transform to real-valued signal
%   This function transforms a complex-valued set of signals x into a
%   real-valued set of signals. For example, if the input signal is 
%   [xI+1j*xQ,yI+1j*yQ] it becomes [xI,xQ,yI,yQ]
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process
%
%   OUTPUTS
%   y           :=  Array of processed signals (column vectors)

%   Mar. 2016 - Dario Pilori <dario.pilori@polito.it>

y = NaN(size(x,1),2*size(x,2));                                            % preallocate output
y(:,1:2:end) = real(x);                                                    % real part
y(:,2:2:end) = imag(x);                                                    % imaginary part