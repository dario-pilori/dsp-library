function [Results,be,x] = decoding_and_ber(x,a,Parameters,varargin)
%DECODING_AND_BER     Minimum distance deconding and BER computation (legacy)
%   This function aligns the received waveform with the transmit pattern,
%   decodes with hard decision the transmitted symbols, and calculates the
%   BER according to the given bit mapping.
%
%   THIS FUNCTION IS DEPRECATED
%   As a replacement, use decode.m
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process at 1SPS
%   a           :=  Transmit sequence
%   Parameters  :=  Parameter structure (see GET_DEFAULT_PARAMS)
%   i_TX        :=  Indexes of transmit symbols (optional)
%
%   OUTPUTS
%   Results     :=  Results structure [Ntx x Nrx]
%   be          :=  Raw bit errors [Ns x Mc x Ntx]
%   x           :=  Fully aligned received waveform

%   Apr. 2016 - Dario Pilori <dario.pilori@polito.it>

warning('This function is deprecated, see decode.m for its replacement.');

%% Retrieve parameters
Mc = Parameters.Mconst;                                                    % bits/symbol
pp = Parameters.ber_pp;                                                    % number of pattern periods to calculate the BER
C = Parameters.C;                                                          % constellation
B = Parameters.B;                                                          % bit mapping
differential = Parameters.differential;                                    % apply differential decoding
Ntx = Parameters.Ntx;                                                      % number of transmitters
min_dist = true;                                                           % minimum-distance decoding or boundary-based
if nargin>=4
    i_TX = varargin{1};                                                    % get indeces of transmit symbols
    validateattributes(i_TX,{'numeric'},{'2d','size',size(a),'integer','>=',1,'<=',2^Mc},'','i_TX',3);
else
    i_TX = [];
end

%% Check input parameters
validateattributes(a,{'double'},{'2d'},'','a',2);
validateattributes(x,{'double'},{'2d','ncols',size(a,2)},'','x',1);
validateattributes(Mc,{'numeric'},{'scalar','positive','integer'},'Parameters','Mconst',2);
validateattributes(C,{'numeric'},{'column','nrows',2^Mc},'Parameters','C',2);
validateattributes(B,{'numeric'},{'2d','binary','nrows',2^Mc,'ncols',Mc},'Parameters','B',2);
validateattributes(pp,{'numeric'},{'scalar','positive','integer'},'Parameters','ber_pp',2);
validateattributes(size(x,1),{'numeric'},{'>=',pp*size(a,1)},'','size(x)',1);
validateattributes(differential,{'logical'},{'scalar'},'Parameters','differential',2);
validateattributes(Ntx,{'numeric'},{'scalar','>=',1,'integer'},'Parameters','Ntx',1);

%% Initialize
Nrx = size(x,2);                                                           % number of sequences to process
Ns = pp*size(a,1);                                                         % length of BER-counting sequence
x = x(end-Ns+1:end,:);                                                     % select last Ns symbols
min_dist = min_dist||rem(Mc,2);                                            % if using a strange constellation, minimum-distance is mandatory

%% Differential decoding
if differential
    x = differential_decoding(x);                                          % apply differential decode to TX
    a = differential_decoding(a);                                          % apply differential decode to RX
end

%% Align and rotate sequences
x = full_align(x,a);                                                       % align sequences

%% Decoding and error counting
if ~min_dist                                                               % check if no minimum-distance
    d = harddec(x,Mc);                                                     % hard decision of the rotated received symbols
end
be = false(Ns,Mc,Nrx);                                                     % allocate bit error matrix
mi = NaN(Nrx,1);                                                           % allocate mutual information matrix
gmi = NaN(Nrx,1);                                                          % allocate generalized mutual information matrix
snr = NaN(Nrx,1);                                                          % allocate Es/No (dB)
nci = NaN(Nrx,1);                                                          % non-circularity index (dB)
for n = 1:Nrx                                                              % for each transmitter
    if isempty(i_TX)
        [~,ax] = max((repmat(a(:,n),1,size(C,1))==repmat(C.',size(a,1),1)).'); % find indices of transmit symbols
    else
        ax = i_TX(:,n);                                                    % get indeces of transmit symbols
    end
    if min_dist
        [mi(n),snr(n),dx,gmi(n)] = MI_eval(x(:,n),ax,Mc,C,B,differential); % calculate mutual information and indeces of symbols
    else
        [~,dx] = max((repmat(int8(d(:,n)),1,size(C,1))==repmat(int8(C.'),size(d,1),1)).'); % find indices of receive symbols
        [mi(n),snr(n)] = MI_eval(x(:,n),ax,Mc,C,B,differential);           % calculate mutual information
    end
    be(:,:,n) = B(repmat(ax,1,pp),:)~=B(dx,:);                             % bit errors per symbol
    nci(n) = eval_nlpn(x(:,n),a(:,n),0);                                   % calculate non-circularity index
end                                                                        %
ne  = squeeze(sum(sum(be,2),1));                                           % total number of bit errors
ber = ne/size(B,2)/Ns;                                                     % bit error rate

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
