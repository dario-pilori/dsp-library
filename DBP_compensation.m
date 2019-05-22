function x = DBP_compensation(x,Parameters,varargin)
%DBP_COMPENSATION     Non-Linear DBP compensation
%   This function tries to compensate for self-channel non-linear
%   interference (and chromatic dispersion) using a simple Digital
%   Back-Propagation (DBP) algorithm.
%
%   INPUTS:
%   x           :=  Array of complex-valued signals (column vectors) to process at 2SPS
%   Parameters  :=  Parameter structure (see GET_DEFAULT_PARAMS)
%   SPS         :=  Number of samples per symbol of x (optional)
%   cut         :=  Cut output (optional)
%
%   OUTPUTS
%   y           :=  Array of processed signals (column vectors) at 2SPS

%   May 2018 - Dario Pilori <dario.pilori@polito.it>

%% Get optional input
if nargin>2
    SPS = varargin{1};
else
    SPS = 2;
end
if nargin>3
    cut = varargin{2};
else
    cut = true;
end

%% Retrieve parameters from stucture
fs   = SPS*Parameters.Rs;                                                  % sampling rate
b2 = Parameters.b2;                                                        % chromatic dispersion in (s/(Hz m))
Lspan = Parameters.Lspan;                                                  % span length (m)
Pch = Parameters.Pch;                                                      % per-channel power (dBm)
Nspan = Parameters.Nspan;                                                  % number of spans
adB = Parameters.adB;                                                      % fiber attenuation at the reference frequency (dB/m)
csi = Parameters.csi;                                                      % Gamma scaling parameter
gam = Parameters.gam;                                                      % non-linear coefficent (1/(W m))
Nstep  = Parameters.DBPsteps;                                              % number of steps per span
if Nspan<=0
    return;
end

%% Check input parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(fs,{'numeric'},{'scalar','positive'},'Parameters','Rs',2);
validateattributes(SPS,{'numeric'},{'scalar','>=',2},'Parameters','SPS',2);
validateattributes(b2,{'numeric'},{'scalar'},'Parameters','b2',2);
validateattributes(Lspan,{'numeric'},{'scalar','positive'},'Parameters','Lspan',2);
validateattributes(adB,{'numeric'},{'scalar','nonnegative'},'Parameters','adB',2);
validateattributes(gam,{'numeric'},{'scalar'},'Parameters','gam',2);
validateattributes(b2,{'numeric'},{'scalar'},'Parameters','b2',2);
validateattributes(Pch,{'numeric'},{'scalar'},'Parameters','Pch',2);
validateattributes(Nspan,{'numeric'},{'scalar','integer','positive'},'Parameters','Nspan',2);
validateattributes(csi,{'numeric'},{'scalar','nonnegative'},'Parameters','csi',2);
validateattributes(Nstep,{'numeric'},{'scalar','integer','>=',1},'Parameters','DBPsteps',2);

%% Evaluate CD delay
L_orig = size(x,1);                                                        % size of input signal
max_d = 4*Nspan*abs(b2)*pi*(fs^2)*Lspan;                                   % total delay per span (samples)
L_pad = 2^nextpow2(L_orig+1+max_d);                                        % length of filter
if cut 
    L_out = floor(L_orig-max_d);                                           % output length
else
    L_out = L_orig;
end

%% Pad signal to next power of two
padding = repmat(x,ceil((L_pad-L_orig)/2/L_orig),1);                       % padding (repeated signal)
x_pad = [padding(1+size(padding,1)-floor((L_pad-L_orig)/2):end,:);...
    x;...
    padding(1:ceil((L_pad-L_orig)/2),:)];

%% Frequency axis
f_step = fs/L_pad;                                                         % frequency resolution of DBP
f = fftshift((-fs/2:f_step:fs/2-f_step)).';                                % frequency axis

%% Calculate step lengths
att = adB/20/log10(exp(1));                                                % attenuation in Np/m
delta = (1-exp(-2*att*Lspan))/Nstep;                                       % effective span length
L1 = log((1-((1:Nstep)'-1)*delta)./(1-(1:Nstep)'*delta))/2/att;            % optimal step length

%% Back-propagate
for ii = 1:Nspan                                                           % for each span
    P = 10^(Pch/10)*1e-3*exp(-2*att*Lspan);                                % calculate power at end of span
    for k = 1:Nstep                                                        % for each step inside span
        %% Calculate steps
        L_step = L1(Nstep-k+1);                                            % step length
        Leff=(1-exp(-2*att*L_step))/2/att;                                 % effective length of the step
        
        %% Linear step
        x_pad = ifft(fft(x_pad).*exp(+1j*2*pi^2*b2*L_step*f.^2));          % chromatic dispersion
        
        %% Non-linear step
        P = P*exp(2*att*L_step);                                           % scale power at beginning of step
        x_pad = x_pad.*exp(1j*csi*gam*Leff*P/SPS*sum(abs(x_pad).^2,2));        % non-linear step
    end
end

%% Remove padding
x = double(x_pad(floor((L_pad-L_out)/2)+1:L_pad-ceil((L_pad-L_out)/2),:)); % remove padding and some samples
