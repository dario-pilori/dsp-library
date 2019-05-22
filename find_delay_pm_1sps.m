function [t0,xc] = find_delay_pm_1sps(x,a)
%FIND_DELAY_PM_1SPS     Find delay between two sequences of symbols
% Sequences are complex, align is done separtely to real and imag part
%
%   INPUTS:
%   x           :=  Array of symbols (column vectors) to process
%   a           :=  Transmit sequence
%
%   OUTPUTS
%   s0          :=  Delay in symbols
%   xc          :=  Averaged cross-correlation


%% Check input parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(a,{'numeric'},{'2d'},'','a');
validateattributes(size(x,2),{'numeric'},{'>=',size(a,2)},'size','Nrx');

%% Prepare for synchronization
Lf = size(a,1);                                                            % Length of training frame
Nrx = size(x,2);                                                           % Number of RX signals
Ntx = size(a,2);                                                           % Number of TX signals
xc = NaN(Lf,Nrx*Ntx);                                                      % allocate correlation matrix

%% Synchronize

for n = 1:Ntx                                                              % for each real/imag transmitter
    an = a(:,n);                                                           % n-th TX sequence
    for m = 1:size(x,2)                                                    % for each real/transmitter receiver
        rm = x(1:Lf,m);                                                  % m-th RX sequence
        xc(:,(n-1)*Ntx+m) = corrx(an-mean(an),rm-mean(rm));                % cross-correlate each transmitter with each receiver
    end                                                                    %
end                                                                        %
xc = abs(sum(xc,2));                                                       % sum all cross-correlations
[~,t0] = max(xc);                                                          % find peak
end

%% Helper function: cross correlation
% This function is similar to MATLAB's xcorr, but much simpler (and faster)
function y = corrx(x,h)
N = max(size(x,1),size(h,1));                                              % size of circular convolution
y = ifft(fft(x,N).*fft(flip(conj(h),1),N));                                % result is the IFFT of element-wise product
end