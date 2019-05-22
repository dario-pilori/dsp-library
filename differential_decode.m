function [Results,be,x] = differential_decode(x,a,Parameters)
%DIFFERENTIAL_DECODE     Minimum distance differential deconding and BER computation
%   This function aligns the received waveform with the transmit pattern,
%   decodes with hard decision the transmitted symbols, and calculates the
%   BER according to the given bit mapping. To avoid cycle slips, this
%   function encodes differentially the two most significant bits (i.e. the
%   quadrant). This function is for HARD DECISION ONLY! Therefore, no MI
%   nor GMI will be calculated.
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process at 1SPS
%   a           :=  Transmit sequence
%   Parameters  :=  Parameter structure (see GET_DEFAULT_PARAMS)
%
%   OUTPUTS
%   Results     :=  Results structure [Ntx x Nrx]
%   be          :=  Raw bit errors [Ns x Mc x Ntx]
%   x           :=  Fully aligned received waveform

%   Apr. 2016 - Dario Pilori <dario.pilori@polito.it>

%% Check if differential
differential = Parameters.differential;                                    % apply differential decoding
validateattributes(differential,{'logical'},{'scalar'},'Parameters','differential',2);
if ~differential
    error('This function is for differential decoding ONLY! Please use decode.m for any other purpose');
end

%% Retrieve parameters
Mc = Parameters.Mconst;                                                    % bits/symbol
pp = Parameters.ber_pp;                                                    % number of pattern periods to calculate the BER
C = Parameters.C;                                                          % constellation
Ntx = Parameters.Ntx;                                                      % number of transmitters
min_dist = false;                                                           % minimum-distance decoding or boundary-based

%% Check input parameters
validateattributes(a,{'double'},{'2d'},'','a',2);
validateattributes(x,{'double'},{'2d','ncols',size(a,2)},'','x',1);
validateattributes(Mc,{'numeric'},{'scalar','positive','integer'},'Parameters','Mconst',2);
validateattributes(C,{'numeric'},{'column','nrows',2^Mc},'Parameters','C',2);
validateattributes(pp,{'numeric'},{'scalar','positive','integer'},'Parameters','ber_pp',2);
validateattributes(size(x,1),{'numeric'},{'>=',pp*size(a,1)},'','size(x)',1);
validateattributes(Ntx,{'numeric'},{'scalar','>=',1,'integer'},'Parameters','Ntx',1);

%% Initialize
Nrx = size(x,2);                                                           % number of sequences to process
Ns = pp*size(a,1);                                                         % length of BER-counting sequence
x = x(end-Ns+1:end,:);                                                     % select last Ns symbols
min_dist = min_dist||rem(Mc,2);                                            % if using a strange constellation, minimum-distance is mandatory

%% Differential decoding
x = differential_decoding(x);                                          % apply differential decode to TX
a = differential_decoding(a);                                          % apply differential decode to RX

%% Align and rotate sequences
x = full_align(x,a);                                                       % align sequences

%% Decoding and error counting
if ~min_dist                                                               % check if no minimum-distance
    d = harddec(x,Mc);                                                     % hard decision of the rotated received symbols
end
be = false(Ns,Mc,Nrx);                                                     % allocate bit error matrix
mi = NaN(Nrx,1);                                                           % allocate mutual information matrix
gmi = NaN(Nrx,1);                                                          % allocate generalized mutual information matrix
nci = NaN(Nrx,1);                                                          % non-circularity index (dB)
for n = 1:Nrx                                                              % for each transmitter
    [~,ax] = max(a(:,n)==C.',[],2);                                        % find indices of transmit symbols
    if min_dist
        [~,dx] = min(abs(x(:,n)-C.'),[],2);                                % find indices of receive symbols
    else
        [~,dx] = max(d(:,n)==C.',[],2);                                    % find indices of receive symbols
    end
    be(:,:,n) = index_to_binary(repmat(ax,1,pp))~=...
        index_to_binary(dx);                                               % bit errors per symbol
    nci(n) = eval_nlpn(x(:,n),a(:,n),0);                                   % calculate non-circularity index
end                                                                        %
ne  = squeeze(sum(sum(be,2),1));                                           % total number of bit errors
ber = ne/log2(size(C,1))/Ns;                                               % bit error rate
snr = 10*log10(mean(abs(C).^2))-10*log10(var(x-repmat(a,pp,1))).';

%% Build results structure
Results.ber = reshape(ber,Ntx,[]);                                         % average BER per received signal
Results.ne  = reshape(ne,Ntx,[]);                                          % number of bit errors
Results.snr = reshape(snr,Ntx,[]);                                         % SNR (ES/N0, dB)
Results.evm = reshape(10.^(-snr/20+2),Ntx,[]);                             % EVM in percent
Results.MI  = reshape(mi,Ntx,[]);                                          % mutual information (bits)
Results.GMI = reshape(gmi,Ntx,[]);                                         % generalized mutual information (bits)
Results.nci = reshape(nci,Ntx,[]);                                         % non-circularity index (dB)
end

%% Helper function: full alignment
% This function fully aligns the signal x with the signal a (amplitude and
% phase), taking into account transmitter conjugation
% x has to be fully equalized with carrier recovery, at 1SPS
function y = full_align(x,a)
Nrx = size(x,2);                                                           % number of signals to receive
ty = NaN(Nrx,1);                                                           % delays
ry = NaN(Nrx,1);                                                           % rotations
rc = NaN(Nrx,1);                                                           % modulator conjugation
y = NaN(size(x));                                                          % allocate aligned sequence
for n = 1:Nrx                                                              % for each transmitter
    if size(a,2)==1
        an = a;                                                            % TX sequence (all the same)
    else
        an = a(:,n);                                                       % n-th TX sequence
    end
    rn = x(1:size(a,1),n);                                                 % n-th RX sequence
    cx = corrx(rn,an);                                                     % crosscorrelate received and transmitted data
    cc = corrx(conj(rn),an);                                               % crosscorrelate conjugated-received and transmitted data    
    [mx,tx] = max(abs(cx));                                                % find peak
    [mc,tc] = max(abs(cc));                                                % find peak    
    if mx >= mc                                                            % if non-conjugated ...
        rc(n) = 0;                                                         % conjugate = 0
        ty(n) = tx;                                                        % pattern delay
        ry(n) = sign(round(exp(1j*angle(cx(tx)))));                        % phase rotation
        y(:,n) = circshift(x(:,n)*conj(ry(n)),size(a,1)-ty(n));            % circularly shift and rotate the received data
    else                                                                   % if conjugate ...
        rc(n) = 1;                                                         % conjugate = 1
        ty(n) = tc;                                                        % pattern delay
        ry(n) = sign(round(exp(1j*angle(cc(tc)))));                        % phase rotation
        y(:,n) = circshift(conj(x(:,n))*conj(ry(n)),-ty(n));               % circularly shift, rotate, and conjugate the received data
    end                                                                    %
end                                                                        %
end

%% Helper function: cross correlation
% This function is similar to MATLAB's xcorr, but much simpler (and faster)
function y = corrx(x,h)
N = max(size(x,1),size(h,1));                                              % size of circular convolution
y = ifft(fft(x,N).*fft(flip(conj(h),1),N));                                % result is the IFFT of element-wise product
end

%% Helper function: differential decoding
% This function differentially-decode the signal x
function xd = differential_decoding(x)
quad = sign(real(x))+1j*sign(imag(x));                                     % find quadrant of RX symbols
xd = x.*[quad(2:end,:);(1+1j)*ones(1,size(x,2))].*conj(quad).^2*(1+1j)/4;  % differentially-decode RX symbols
end
