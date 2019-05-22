function [y,ph,ph_pilot] = phase_recovery_pilot(x,a,Parameters)
%PHASE_RECOVERY_PILOT     BPS and ML phase recovery with pilots
%   This function recovers the phase of the input signal x, using a
%   chain of BPS, followed by ML estimation. In order to avoid cycle
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
%   ph          :=  Recovered phase (rad)

%   Jan 2017 - Dario Pilori <dario.pilori@polito.it>

%% Retrieve parameters
N = Parameters.cpe_memory;                                                 % CPE memory (samples)
stp = Parameters.bps_anglestep;                                            % BPS angle step (radians)
C = Parameters.C;                                                          % Constellation
pL = Parameters.cpe_pilotlen;                                              % pilot length
bS = Parameters.cpe_blocksize;                                             % pilot block size
sps = Parameters.cpe_sps;                                                  % Input samples/symbol
ml_dec = Parameters.cpe_ml_dec;                                            % use ML decision instead of hard decision (slower)
Ntx = Parameters.Ntx;                                                      % number of TX sequences
sym = Parameters.cpe_sym;                                                  % constellation symmetry (2*pi, pi, pi/2)
ph = 1;                                                                    % sampling phase for 2SPS processing

%% Check parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(a,{'double'},{'2d'},'','a',2);
validateattributes(C,{'numeric'},{'column'},'Parameters','C',3);
validateattributes(N,{'numeric'},{'scalar','>=',1,'integer'},'Parameters','cpe_memory',3);
validateattributes(2*pi/sym,{'numeric'},{'scalar','>=',1,'integer'},'Parameters','sym',3);
validateattributes(stp,{'numeric'},{'scalar','positive','<',sym},'Parameters','bps_anglestep',3);
validateattributes(pL,{'numeric'},{'scalar','nonnegative','integer'},'Parameters','cpe_pilotlen',3);
validateattributes(bS,{'numeric'},{'scalar','>=',pL,'integer'},'Parameters','cpe_blocksize',3);
validateattributes(sps,{'numeric'},{'scalar','>=',1,'<=',2,'integer'},'Parameters','sps',3);
validateattributes(ml_dec,{'logical'},{'scalar'},'Parameters','ml_dec',3);
validateattributes(Ntx,{'double'},{'scalar','integer','positive'},'Parameters','Ntx',3);
if size(C,1)<=4
    error('Pilot-aided phase recovery does not work with QPSK/BPSK. Please use Viterbi-Viterbi.');
end

%% Check if pilot is needed
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

%% Initialize 
M = round(sym/stp);                                                        % calculate number of test angles
pt = (0:M-1).'/M*sym-sym/2;                                                % test phase vector in radians
ph = NaN(size(x));                                                         % allocate recovered phase

%% For all received waveforms
for n = 1:size(x,2)                                                        % for all received waveforms
    xn = x(:,n);                                                           % get n-th waveform
    
    %% Blind Phase Search
    x_rot = xn.*exp(1j*pt).';                                              % apply test phases to received waveform
    if ml_dec
        er = movmean(abs(x_rot-ML_dec(x_rot,C)).^2,N);                     % averaged squared distance from closest point
    else
        er = movmean(abs(x_rot-hard_dec(x_rot,C)).^2,N);                   % averaged squared distance from closest point
    end
    [~,idx] = min(er,[],2);                                                % get minimum error
    
    %% Pilot-tone-aided cycle slip cancellation
    if pilot
        % Use pilot tones to correctly unwrap BPS phase
        a_idx = mod(n-1,Ntx)+1;                                            % get correct index of training sequence
        ph_pilot = angle(movmean(reshape(xn,size(a_P(:,a_idx),1),[]).*...
            conj(a_P(:,a_idx)),pL,'Endpoints','fill'));                    % calculate correct phase during pilot blocks
        pilot_pos = find(~isnan(ph_pilot(:,1)));                           % find pilot blocks
        ph_pilot(pilot_pos,:) = reshape(unwrap(reshape(ph_pilot(pilot_pos,:),[],1))...
            ,size(pilot_pos,1),[]);                                        % unwrap pilot phase
        ph_pilot = interp_pilot_phase(ph_pilot(:));                        % interpolate pilot-tones phase
        
        %% Apply correction
        ph_bps = -pt(idx);                                                 % recovered phase after BPS
        ph_bps = ph_bps - sym*round((ph_bps-ph_pilot)/sym);                % detect and correct cycle slips using pilot tones
    else
        % If pilots are not used, unwrap BPS recovered phase
        ph_bps = -unwrap(pt(idx)*2*pi/sym)/(2*pi/sym);                     % recovered phase after BPS
        ph_pilot = NaN(size(x,1),1);                                       % return NaN for pilot phase
    end
    
    %% Maximum Likelihood
    y1 = xn.*exp(-1j*ph_bps);                                              % apply first phase correction
    ph(:,n) = unwrap(angle(movmean(conj(y1).*xn,N)));                      % ML phase recovery
end

%% Apply phase correction
y = x.*exp(-1j*ph);
end

%% Helper function: hard decision
% This function applies hard-decision to the received signal x
function c = hard_dec(x,C)
switch log2(size(C,1))
    case 2
        c = sign(real(x))+1j*sign(imag(x));
    case 4
        c = 2*(sign(real(x))+1j*sign(imag(x)));
        c = c+sign(real(x-c))+1j*sign(imag(x-c));
    case 6
        c = 4*(sign(real(x))+1j*sign(imag(x)));
        c = c+2*(sign(real(x-c))+1j*sign(imag(x-c)));
        c = c+sign(real(x-c))+1j*sign(imag(x-c));
    case 8
        c = 8*(sign(real(x))+1j*sign(imag(x)));
        c = c+4*(sign(real(x-c))+1j*sign(imag(x-c)));
        c = c+2*(sign(real(x-c))+1j*sign(imag(x-c)));
        c = c+sign(real(x-c))+1j*sign(imag(x-c));
    otherwise
        c = ML_dec(x,C);
end
end

%% Helper function: Maximum likelihood detection
function c = ML_dec(x,C)
N = size(x,2);
x = x(:);
[~,idx] = min(abs(x.'-C));
c = C(idx);
c = reshape(c,[],N);
end

%% Helper function: Interpolate pilot phase
function x = interp_pilot_phase(x)
% Taken from: http://it.mathworks.com/matlabcentral/fileexchange/8225-naninterp
x(1) = x(find(~isnan(x),1)); % set first and last point
x(end) = x(find(~isnan(x),1,'last'));
x(isnan(x)) = interp1(find(~isnan(x)),x(~isnan(x)),find(isnan(x)),'linear');% linear interpolation
end

%% Helper function: add pilot symbols
function a_P = add_pilot(a,cpe_pilotlen,blocksize)
a_P = NaN(size(a));
pilot_idx = reshape((0:blocksize:size(a,1)-blocksize)+(1:cpe_pilotlen)',[],1);
for n = 1:size(a,2)
    a_P(pilot_idx,n) = a(pilot_idx,n);
end
end