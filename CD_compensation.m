function x = CD_compensation(x,Parameters,varargin)
%CD_COMPENSATION     Compensate for chromatic dispersion
%   This function applies an inverse chromatic dispersion filter to the
%   complex-valued array of column-vector signals x.
%   The signals are assumed to be at 2 samples per symbol.
%
%   INPUTS:
%   x           :=  Array of complex-valued signals (column vectors) to process at 2SPS
%   Parameters  :=  Parameter structure (see GET_DEFAULT_PARAMS)
%   SPS         :=  Number of samples per symbol of x (optional)
%   cut         :=  Cut output (optional)
%
%   OUTPUTS
%   y           :=  Array of processed signals (column vectors) at 2SPS

%   Mar. 2016 - Dario Pilori <dario.pilori@polito.it>

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
f0   = Parameters.f0;                                                      % central frequency (Hz)
CD   = Parameters.CD;                                                      % cumulated dispersion (s/m) (1 ps/nm = 1e-3 s/m)

%% Check input parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(f0,{'numeric'},{'scalar','positive'},'Parameters','f0',2);
validateattributes(fs,{'numeric'},{'scalar','positive'},'Parameters','Rs',2);
validateattributes(SPS,{'numeric'},{'scalar','>=',2},'Parameters','SPS',2);
validateattributes(CD,{'numeric'},{'scalar'},'Parameters','CD',2);
if all(CD==0)||any(isnan(CD))
    return;
end

%% Build CD filter
c0 = 299792458;                                                            % speed of light (m/s)
N = size(x,1);                                                             % size of input signal
max_d = 2*abs(CD)*c0*(fs/f0)^2*(2/SPS);                                    % CD filter delay within one WDM channel
[~,E] = log2(N+1);                                                         % get logarithm of length of filter
L_filt = 2^E;                                                              % length of filter
if cut 
    L_out = floor(N-2*max_d);                                              % output length
else
    L_out = N;
end
H = exp(-1j*pi*c0*CD*((-L_filt/2:L_filt/2-1)*fs/L_filt/f0).^2).';          % CD frequency response

%% Apply CD filter
x_pad = [zeros(floor((L_filt-N)/2),size(x,2));...                          % pad with zeros
    x;...
    zeros(ceil((L_filt-N)/2),size(x,2))];
x_pad = ifft(ifftshift(bsxfun(@times,fftshift(fft(x_pad)),H)));            % apply filter in frequency domain
x = x_pad(floor((L_filt-L_out)/2)+1:L_filt-ceil((L_filt-L_out)/2),:);      % remove zeros and samples