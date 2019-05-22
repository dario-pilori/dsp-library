function [t0,xc] = find_delay_pm_2sps(x,a)
%FIND_DELAY_PM_2SPS     Find delay between two sequences
%   This function finds the delay t0 between two sequences, x at 2SPS and a
%   at 1SPS. a is assumed to be repeated several times in x, and we assume
%   that a is "mixed" inside x by PMD.
%
%   Put it more simply, use this function to align the sequence received by
%   a DualPol coherent receiver to thre transmit sequence.
%
%   Due to phase noise, is recommended to align the power of the two
%   sequences whenever possible.
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process at 2SPS
%   a           :=  Transmit sequence
%
%   OUTPUTS
%   t0          :=  Delay at 2SPS
%   xc          :=  Averaged cross-correlation

%   Apr. 2017 - Dario Pilori <dario.pilori@polito.it>

%% Check input parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(a,{'numeric'},{'2d'},'','a');
validateattributes(size(x,2),{'numeric'},{'>=',size(a,2)},'size','Nrx');

%% Prepare for synchronization
Lf = size(a,1);                                                            % Length of training frame
Nrx = size(x,2);                                                           % Number of RX signals
Ntx = size(a,2);                                                           % Number of TX signals
a2  = reshape(repmat(reshape(a,[],1),1,2).',[],size(a,2));                 % 2x upsample TX data
xc = NaN(2*Lf,Nrx*Ntx);                                                    % allocate correlation matrix

%% Synchronize
for n = 1:Ntx                                                              % for each complex transmitter
    an = a2(:,n);                                                          % n-th TX sequence
    for m = 1:size(x,2)                                                    % for each complex receiver
        rm = x(1:2*Lf,m);                                                  % m-th RX sequence
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