function [y] = FOC_FD(x,Parameters,a)
%FINE_FOC     Fine Frequency offset compensation
%   This function applies a fine frequency recovery scheme to help adaptive
%   equalization
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process at 2SPS
%   Parameters  :=  Parameter structure (see GET_DEFAULT_PARAMS)
%   a           :=  Array of transmitted signals (column vectors)
%   OUTPUTS
%   y           :=  Array of processed signals (column vectors) at 2SPS

%   Jan. 2018 - Antonino Nespola <nespola@ismb.it>

%% Retrieve parameters from stucture
kDC  = Parameters.kDC;                                                     % memory of DC-block

%% Check input parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(kDC,{'numeric'},{'scalar','positive','integer'},'Parameters','kDC',2);
validateattributes(a,{'double'},{'2d'},'','a');
validateattributes(size(x,2),{'numeric'},{'>=',size(a,2)},'size','Nrx');

%% Prepare for freq offset correction
Lf = size(a,1);                                                            % Length of training frame
a2  = rectpulse(a,2);                                                      % 2x upsample TX data
xt = x(1:2*Lf,:);                                                          % truncate received signals
Df = 0.2;                                                                   % minimum allowed freq offset 10% of sampling rate (2xSymbolRate)
Nf = 10000;                                                                % Number of test freq offset 

%% First coarse f0 offset search
test_f0 = linspace(-Df,Df,Nf);
Peaks = NaN(Nf,1);
for n = 1: Nf
    Peaks(n) = FindPeakCorr(Lf,test_f0(n),xt,a2);
end
[~,idx] = max(Peaks);
f0_fine = test_f0(idx);

%% Apply freq correction
ph = 2*pi*f0_fine.*repmat((0:size(x,1)-1)',1,size(x,2));                   % generate phase                         % generate phase axis
y = x.*exp(-1j*ph);                                                        % shift frequency
y = y-movmean(y,kDC,1);                                                    % remove DC
fprintf('Normalized freq. offset %d\n',f0_fine);
end

%% Helper function: cross correlation
% cross-correlation between  transmitters  and receivers.
% Phase offset is applyed to the last
function peak_xc=FindPeakCorr(Lf,f0,x,a)
Nrx = size(x,2);                                                           % Number of RX signals
Ntx = size(a,2);                                                           % Number of TX signals
xc = NaN(2*Lf,Nrx*Ntx);                                                    % allocate correlation matrix

for n = 1:Ntx                                                              % for each complex transmitter
    an = a(:,n);                                                           % n-th TX sequence
    an = an-mean(an);                                                      % remove DC
    for m = 1:Nrx                                                          % for each complex receiver
        rm = x(:,m);                                                       % m-th RX sequence
        rm = rm.*exp(-2j*pi*f0.*(0:2*Lf-1)');                                % phase correction
        rm =rm-mean(rm);                                                   % remove DC
        xc(:,(n-1)*Ntx+m) = corrx(an,rm);                                  % cross-correlate each transmitter with each receiver
    end
end
xc = abs(sum(xc,2));                                                       % sum all cross-correlations
[peak_xc,~] = max(xc);                                                     % find peak
end

%% Helper function: cross correlation
% This function is similar to MATLAB's xcorr, but much simpler (and faster)
function y = corrx(x,h)
N = max(size(x,1),size(h,1));                                              % size of circular convolution
y = ifft(fft(x,N).*fft(flip(conj(h),1),N));                                % result is the IFFT of element-wise product
end