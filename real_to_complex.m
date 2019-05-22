function y = real_to_complex(x)
%REAL_TO_COMPLEX     Transform to complex-valued signal
%   This function transforms a real-valued set of signals x into a
%   complex-valued set of signals. For example, if the input signal is 
%   [xI,xQ,yI,yQ] it becomes [xI+1j*xQ,yI+1j*yQ] 
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process
%
%   OUTPUTS
%   y           :=  Array of processed signals (column vectors)

%   Mar. 2016 - Dario Pilori <dario.pilori@polito.it>

%% Check input parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(size(x,2),{'numeric'},{'scalar','even'},'','x',1);

%% Real to complex
y = x(:,1:2:end) + x(:,2:2:end)*1j;