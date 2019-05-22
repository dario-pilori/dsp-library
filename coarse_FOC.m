function y = coarse_FOC(x,Parameters)
%COARSE_FOC     Coarse Frequency offset compensation
%   This function applies a coarse frequency recovery scheme to help adaptive
%   equalization
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process at 2SPS
%   Parameters  :=  Parameter structure (see GET_DEFAULT_PARAMS)
%
%   OUTPUTS
%   y           :=  Array of processed signals (column vectors) at 2SPS

%   Mar. 2016 - Dario Pilori <dario.pilori@polito.it>

%% Retrieve parameters from stucture
kDC  = Parameters.kDC;                                                     % memory of DC-block

%% Check input parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(kDC,{'numeric'},{'scalar','positive','integer'},'Parameters','kDC',2);

%% Frequency recovery
pxx = pwelch(x);                                                           % calculate spectrum
[~,idx] = max(sum(pxx,2),[],1);                                            % find peak
f0 = round((idx-1)*size(x,1)/size(pxx,1));                                 % calculate peak normalized frequency
if any(f0/size(x,1)>0.25&f0/size(x,1)<0.75)
    warning('Frequency offset too big. Setting to zero.');
    y = x;
    return;
end
ph = 2*pi*repmat((0:size(x,1)-1)',1,size(x,2))/size(x,1)...
    .*repmat(f0,size(x,1),1);                                              % generate phase axis
y = x.*exp(-1j*ph);                                                        % shift frequency
y = y-movmean(y,kDC,1);                                                    % remove DC