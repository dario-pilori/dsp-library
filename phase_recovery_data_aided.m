function [y,phn] = phase_recovery_data_aided(x,a,Parameters)
%PHASE_RECOVERY_DATA_AIDED     Data-aided phase recovery
%   This function recovers the phase of the input signal x, using a
%   data-aided scheme (cheating). PAY ATTENTION that this function is not
%   able to detect a conjugated input (due to Mach-Zehnder modulator, or
%   blind equalizer)!
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process at 2SPS
%   a           :=  Transmit sequence
%   Parameters  :=  Parameter structure (see GET_DEFAULT_PARAMS)
%
%   OUTPUTS
%   y           :=  Array of processed signals (column vectors) at 2SPS
%   phn         :=  Estimated phase error (radians)

%   June 2016 - Dario Pilori <dario.pilori@polito.it>

%% Retrieve parameters
Nmem = Parameters.cpe_memory;                                              % CPE memory (samples)
sps = Parameters.cpe_sps;                                                  % Input samples/symbol
is_conj = Parameters.cpe_conj;                                             % Conjugate TX pattern
M = Parameters.Mconst;                                                     % Constellation bits
Ntx = Parameters.Ntx;                                                      % number of TX sequences
ph = 1;                                                                    % sampling phase for 2SPS processing

%% Check input parameters
validateattributes(a,{'double'},{'2d'},'','a',2);
validateattributes(x,{'double'},{'2d','ncols',size(a,2)},'','x',1);
validateattributes(Ntx,{'double'},{'scalar','integer','positive'},'Parameters','Ntx',3);
validateattributes(M,{'double'},{'scalar','integer','>',2},'Parameters','Mconst',3);
validateattributes(Nmem,{'numeric'},{'scalar','>',1,'integer'},'Parameters','cpe_memory',3);
validateattributes(is_conj,{'logical'},{'column'},'Parameters','is_conj',3);
validateattributes(sps,{'numeric'},{'scalar','>=',1,'<=',2,'integer'},'Parameters','sps',3);

%% Downsample to 1SPS
x_ds = x(ph:sps:end,:);

%% Align input data
a = power_tx_sequence_align(x_ds(size(a,1)+1:2*size(a,1),:),a);            % power-based align
a = repmat(a,ceil(size(x_ds,1)/size(a,1)),1);                              % repeat sequence
a = a(1:size(x_ds,1),:);                                                   % and cut it
for n = 1:Ntx                                                              % conjugation of transmitters...
    if is_conj(n)
        a(:,n:Ntx:end) = conj(a(:,n:Ntx:end));
    end
end

%% Compute phase error
phn = exp(1j*(atan2(imag(x_ds),real(x_ds))-atan2(imag(a),real(a))));       % calculate phase error
phn = movmean(phn,Nmem);                                                   % moving average
if sps==2                                                                  % if 2SPS...
    phn = reshape(repmat(reshape(phn,[],1),1,2).'...
        ,[],size(phn,2));                                                  % 2x upsample phase error
end

%% Correct phase
y = x.*conj(sign(phn));                                                    % apply correct phase
end

%% Helper function
% This function aligns the transmit pattern a to the received sequence, in
% order to apply phase recovery
function a = power_tx_sequence_align(x,a)
%% Prepare for synchronization
Lf = size(a,1);                                                            % Length of training frame
Nrx = size(x,2);                                                           % Number of RX signals
Ntx = size(a,2);                                                           % Number of TX signals
xc = NaN(Lf,Nrx*Ntx);                                                      % allocate correlation matrix

%% Synchronize
for n = 1:Ntx                                                              % for each complex transmitter
    an = abs(a(:,n));                                                      % power of n-th TX sequence
    for m = 1:size(x,2)                                                    % for each complex receiver
        rm = abs(x(1:Lf,m));                                               % power of m-th RX sequence
        xc(:,(n-1)*Ntx+m) = corrx(an-mean(an),rm-mean(rm));                % cross-correlate each transmitter with each receiver
    end                                                                    %
end                                                                        %
xc = abs(sum(xc,2));                                                       % sum all cross-correlations
[~,t0] = max(xc);                                                          % find peak
a = circshift(a,Lf-t0);                                                    % shift data
end

%% Helper function: cross correlation
% This function is similar to MATLAB's xcorr, but much simpler (and faster)
function y = corrx(x,h)
N = max(size(x,1),size(h,1));                                              % size of circular convolution
y = ifft(fft(x,N).*fft(flip(conj(h),1),N));                                % result is the IFFT of element-wise product
end