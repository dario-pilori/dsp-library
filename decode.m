function [Results,be,x,L] = decode(x,a,Parameters,varargin)
%DECODE     Deconding and performance evaluation
%   This function aligns the received waveform with the transmit pattern,
%   estimates the signal-to-noise ratio and calculate log-likelihood ratios
%   (LLRs). From these, BER and GMI are calculated as performance
%   indicators.
%
%   This is a new version that is going to replace the old
%   decoding_and_ber. This version is not compatible with the old version,
%   and it requires C MEX functions to evaluate LLR. These functions works
%   under 64-bit Linux. Don't know/don't care about other OS.
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
%   L           :=  Log-likelihood ratios [??]

%   May 2018 - Dario Pilori <dario.pilori@polito.it>

%% Retrieve parameters
Mc = Parameters.Mconst;                                                    % bits/symbol
pp = Parameters.ber_pp;                                                    % number of pattern periods to calculate the BER
C = Parameters.C;                                                          % constellation
B = Parameters.B;                                                          % bit mapping
differential = Parameters.differential;                                    % apply differential decoding
calculate_MI = Parameters.calculate_MI;                                    % set to true to calculate MI other than GMI
Ntx = Parameters.Ntx;                                                      % number of transmitters
if nargin>=4
    i_TX = varargin{1};                                                    % get indeces of transmit symbols
    validateattributes(i_TX,{'numeric'},{'2d','integer','>=',1,'<=',2^Mc},'','i_TX',3);
else
    i_TX = [];
end

%% Check input parameters
validateattributes(a,{'double'},{'2d'},'','a',2);
validateattributes(x,{'double'},{'2d','ncols',size(a,2)},'','x',1);
validateattributes(Mc,{'numeric'},{'scalar','positive','integer'},'Parameters','Mconst',2);
M = 2^Mc;                                                                  % constellation size
validateattributes(C,{'numeric'},{'column','nrows',M},'Parameters','C',2);
validateattributes(B,{'numeric'},{'2d','binary','nrows',M,'ncols',Mc},'Parameters','B',2);
if ~issorted(bi2de(B,'left-msb'))
    warning('Bits are not sorted; this is not good');
    [~,idx] = sort(bi2de(B,'left-msb'));
    B = B(idx,:);
    C = C(idx);
end
validateattributes(pp,{'numeric'},{'scalar','positive','integer'},'Parameters','ber_pp',2);
validateattributes(size(x,1),{'numeric'},{'>=',pp*size(a,1)},'','size(x)',1);
validateattributes(differential,{'logical'},{'scalar'},'Parameters','differential',2);
if differential
    error('Differential decoding not yet supported; please use old function');
end
validateattributes(calculate_MI,{'logical'},{'scalar'},'Parameters','calculate_MI',2);
validateattributes(Ntx,{'numeric'},{'scalar','>=',1,'integer'},'Parameters','Ntx',1);

%% Initialize sequences
Nrx = size(x,2);                                                           % number of sequences to process
Ns = pp*size(a,1);                                                         % length of BER-counting sequence
x = x(end-Ns+1:end,:);                                                     % select last Ns symbols

%% Align and rotate sequences
x = full_align(x,a);                                                       % align sequences

%% Decoding and error counting
be = false(Ns,Mc,Nrx);                                                     % allocate bit error matrix
mi = NaN(Nrx,1);                                                           % allocate mutual information matrix
gmi = NaN(Nrx,1);                                                          % allocate generalized mutual information matrix
gmi_pn = NaN(Nrx,1);                                                       % allocate phase-noise generalized mutual information matrix
snr = NaN(Nrx,1);                                                          % allocate Es/No (dB)
nci = NaN(Nrx,1);                                                          % non-circularity index (dB)
for n = 1:Nrx                                                              % for each transmitter
    y = x(:,n);                                                            % get n-th received waveform
    
    %% Generate indeces of transmit symbols
    if isempty(i_TX)
        [~,ax] = max((repmat(a(:,n),1,size(C,1))==repmat(C.',size(a,1),1)).'); % find indices of transmit symbols
    else
        ax = i_TX(:,n).';                                                  % get indeces of transmit symbols
    end
    ax = repmat(ax,1,pp);                                                  % repeat transmit indeces pp times
    b_TX = de2bi(ax-1,'right-msb');                                        % get transmit bit sequence
    
    %% Get symbol probability and normalize power
    [P_symb,pos] = evaluate_probability(ax,M,Ns);
    y = y/std(y)*sqrt(sum(P_symb.*abs(C).^2));

    %% Calculate SNR
    [centroids,sigma2] = centroid_and_sigma(y,pos,M);                      % calculate centroids and noise variance
    entr = -P_symb.'*log2(P_symb);                                         % calculate constellation entropy
    snr(n) = 10*log10(sum(P_symb.*abs(centroids).^2)/median(sigma2));      % calculate average SNR (dB)
        
    %% Calculate MI
    if calculate_MI
        if isreal(C)
            mi(n) = pam_mi_montecarlo_mex(centroids, median(sigma2),...
                centroids(ax)-y, P_symb);                                  % real AWGN Mutual Information
        else
            mi(n) = qam_mi_montecarlo_mex(centroids, median(sigma2),...
                centroids(ax)-y, P_symb);                                  % complex AWGN Mutual Information
        end
    end
    
    %% Calculate LLRs
    if isreal(C)
        %% PAM LLLR
        L = reshape(pam_llr_mex(centroids, sigma2, y, P_symb),...
            log2(M),Ns).';                                                 % calculate real AWGN log-likelihood ratios
        snr(n) = snr(n) - 10*log10(2);                                     % use SNR definition for real-valued signals
    else
        %% QAM AWGN LLR
        L = reshape(qam_llr_mex(centroids, sigma2, y, P_symb),...
            log2(M),Ns).';                                                 % calculate complex AWGN log-likelihood ratios
        
        %% Estimate Von Mises parameters
        [sigma_pn,s2_est] = estimate_tikhonov(y,ax,C,sigma2);              % estimate Tikhonov distribution parameters
        
        %% QAM PN LLR
        L_pn = reshape(qam_llr_pn_mex(centroids,1/s2_est,1/sigma_pn^2,...
            y,P_symb),log2(M),Ns).';                                       % calculate LLR with full formula
        if any(isnan(L_pn(:)))                                             % in case of numerical errors...
            L_pn = reshape(qam_llr_pn_maxlog_mex(centroids,1/s2_est,...
                1/sigma_pn^2,y,P_symb),log2(M),Ns).';                      % use MAX-LOG approximation!
        end
    end
    
    %% Calculate GMI
    [~,tmp_gmi] = fminbnd(@(s) ...
        sum(sum(log2(1+exp((-1).^(~b_TX).*L*s))))/Ns,0.01,10);             % calculate GMI using MonteCarlo from LLR
    gmi(n) = max(0,entr - tmp_gmi);                                        % force non-negative GMI
    if ~isreal(C)
        [~,tmp_gmi] = fminbnd(@(s) ...
            sum(sum(log2(1+exp((-1).^(~b_TX).*L_pn*s))))/Ns,0.01,10);      % calculate GMI using MonteCarlo from LLR
        gmi_pn(n) = max(0,entr - tmp_gmi);                                 % force non-negative GMI
    end
    
    %% Calculate bit errors using the LLRs
    be(:,:,n) = b_TX ~= (-sign(L)+1)/2;                                    % perform hard-decision on LLR to find pre-FEC BER
    
    %% Calculate phase noise
    single_p = y./centroids(ax);                                           % condense constellation to a single point
    nci(n) = 10*log10(var(imag(single_p))./var(real(single_p)));           % evaluate non-circularity index
end
ne  = squeeze(sum(sum(be,2),1));                                           % total number of bit errors
ber = ne/size(B,2)/Ns;                                                     % bit error rate

%% Build results structure
Results.ber = reshape(ber,Ntx,[]);                                         % average BER per received signal
Results.ne  = reshape(ne,Ntx,[]);                                          % number of bit errors
Results.snr = reshape(snr,Ntx,[]);                                         % SNR (ES/N0, dB)
Results.evm = reshape(10.^(-snr/20+2),Ntx,[]);                             % EVM in percent
Results.MI  = reshape(mi,Ntx,[]);                                          % mutual information (bits)
Results.GMI = reshape(gmi,Ntx,[]);                                         % generalized mutual information (bits)
Results.GMIpn = reshape(gmi_pn,Ntx,[]);                                    % generalized mutual information (bits)
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

%% Helper function: probability calculation
function [P_symb,pos] = evaluate_probability(iTX,M,Ns)
pos = int16(iTX)==(1:M).';                                             % indexes of each point
P_symb = sum(pos,2)/Ns;                                                % probability of point
if any(P_symb<=eps)
    error('The transmit sequence does not contain all symbols. Please modify the constellation.');
end
end

%% Helper function: centroids and SNR calculation
function [centroids,sigma2] = centroid_and_sigma(x,pos,M)
centroids = NaN(M,1);                                                  % preallocate centroids
sigma2 = NaN(M,1);                                                     % preallocate variances
for i = 1:M                                                            % for all transmitted symbols
    rx_s = x(pos(i,:));                                                % get received points relative to i-th TX symbol
    centroids(i) = mean(rx_s);                                         % calculate centroid
    sigma2(i) = max(eps,var(rx_s));                                    % calculate variance (never to zero)
end
end

%% Helper function: Tikhonov distribution parameters
function [sigma_pn,sigma_n] = estimate_tikhonov(x,iTX,C,sigma2)
z = real(mean(x./C(iTX)));
if z>=1
    sigma_pn = 1e-4;
else
    sigma_pn = sqrt(8*(1-z)/(3+z));
end
sigma_n = mean(sigma2(abs(C)==min(unique(abs(C)))));
end

