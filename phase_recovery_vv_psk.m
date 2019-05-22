function [pr,errphi_filt,ph_pilot] = phase_recovery_vv_psk(x,a,Parameters)
%PHASE_RECOVERY_VV_PSK     V&V PHASE RECOVERY FOR PSK
%   This function recovers the phase of the input signal x, using a
%   Viterbi&Viterbi phase estimation. In order to avoid cycle
%   slips, which happen often with high symbol error rates, periodic pilot
%   symbols are inserted.
%
%   The training sequence is made of blocks of cpe_pilotlen symbols, spaced
%   by cpe_blocksize.
%
%   References: 10.1109/LPT.2010.2049644, 10.1364/OFC.2014.Th4D.1
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process (1SPS)
%   a           :=  Transmit sequence
%   Parameters  :=  Parameter structure (see GET_DEFAULT_PARAMS)
%
%   OUTPUTS
%   y           :=  Array of processed signals (column vectors) at (1SPS)
%   errphi_filt :=  Recovered phase (rad)

%   Dec 2017 - Dario Pilori <dario.pilori@polito.it>

%% Retrieve parameters
N = Parameters.cpe_memory;                                                 % CPE memory (samples)
pL = Parameters.cpe_pilotlen;                                              % pilot length
bS = Parameters.cpe_blocksize;                                             % pilot block size
sps = Parameters.cpe_sps;                                                  % Input samples/symbol
Ntx = Parameters.Ntx;                                                      % number of TX sequences
C = Parameters.C;                                                          % Constellation

%% Derived params
Mapsk = 2*pi/Parameters.cpe_sym;                                           % M-th power Viterbi&Viterbi
ph = 1;                                                                    % sampling phase for 2SPS processing
phi0 = min(abs(angle(C)));                                                 % get minimum phase

%% Check parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(a,{'double'},{'2d'},'','a',2);
validateattributes(N,{'numeric'},{'scalar','>=',1,'integer'},'Parameters','cpe_memory',3);
validateattributes(Mapsk,{'numeric'},{'scalar','>=',1,'integer'},'Parameters','sym',3);
validateattributes(C,{'numeric'},{'column'},'Parameters','C',3);
validateattributes(pL,{'numeric'},{'scalar','nonnegative','integer'},'Parameters','cpe_pilotlen',3);
validateattributes(bS,{'numeric'},{'scalar','>=',pL,'integer'},'Parameters','cpe_blocksize',3);
validateattributes(sps,{'numeric'},{'scalar','>=',1,'<=',2,'integer'},'Parameters','sps',3);
validateattributes(Ntx,{'double'},{'scalar','integer','positive'},'Parameters','Ntx',3);

%% Check if pilots are needed
Ntx = size(a,2);
if pL==0||bS==0
    pilot = false;
else
    pilot = true;
end

%% Downsample to 1SPS
x = x(ph:sps:end,:);

%% Prepare pilot sequence
a = repmat(a,ceil(size(x,1)/size(a,1)),1);                                 % repeat sequence
a = a(1:size(x,1),:);                                                      % and cut it
a_P = add_pilot(a,pL,bS);                                                  % prepare pilot sequence

%% M-th power
yM = sign(x.^Mapsk);                                                       % remove M-PSK and amplitude modulation

%% Estimate phase
yM = movmean(yM,N);                                                        % moving average of 4-th power signal
errphi_filt = angle(yM);                                                   % extract and unwrap phase

%% Pilot-tone-assisted cycle slip cancellation
if pilot
    errphi_filt = errphi_filt - phi0*Mapsk;                                % remove initial phase
    for n = 1:Ntx
        xn = x(:,n);
        a_idx = mod(n-1,Ntx)+1;                                            % get correct index of training sequence
        ph_pilot = angle(movmean(reshape(xn,size(a_P(:,a_idx),1),[]).*...
            conj(a_P(:,a_idx)),pL,'Endpoints','fill'));                    % calculate correct phase during pilot blocks
        pilot_pos = find(~isnan(ph_pilot(:,1)));                           % find pilot blocks
        ph_pilot(pilot_pos,:) = reshape(unwrap(reshape(ph_pilot(pilot_pos,:),[],1))...
            ,size(pilot_pos,1),[]);                                        % unwrap pilot phase
        ph_pilot = interp_pilot_phase(ph_pilot(:));                        % interpolate pilot-tones phase
        
        %% Apply correction
        errphi_filt(:,n) = errphi_filt(:,n) - ...
            2*pi*round((errphi_filt(:,n)-Mapsk*ph_pilot)/2/pi);            % detect and correct cycle slips using pilot tones
    end
else
    errphi_filt = unwrap(errphi_filt)-phi0*Mapsk;                          % blind unwrapping
    ph_pilot = NaN(size(errphi_filt));                                     % set to NaN pilot phase
end
errphi_filt = errphi_filt/Mapsk;                                           % divide phase by M

%% Correct phase
pr = x .* exp(-1j*errphi_filt);                                            % apply phase correction
end

%% Helper function: add pilot symbols
function a_P = add_pilot(a,cpe_pilotlen,blocksize)
a_P = NaN(size(a));
pilot_idx = reshape((0:blocksize:size(a,1)-blocksize)+(1:cpe_pilotlen)',[],1);
for n = 1:size(a,2)
    a_P(pilot_idx,n) = a(pilot_idx,n);
end
end

%% Helper function: Interpolate pilot phase
function x = interp_pilot_phase(x)
% Taken from: http://it.mathworks.com/matlabcentral/fileexchange/8225-naninterp
x(1) = x(find(~isnan(x),1)); % set first and last point
x(end) = x(find(~isnan(x),1,'last'));
x(isnan(x)) = interp1(find(~isnan(x)),x(~isnan(x)),find(isnan(x)),'linear');% linear interpolation
end