function y = resample_signal(x,Parameters)
%RESAMPLE_SIGNAL     Resample the signal at 2SPS
%   This function resamples the input signal at 2
%   samples per symbol
%
%   It requires resample, introduced before MATLAB R2006a in the
%   Signal Processing Toolbox
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process
%   Parameters  :=  Parameter structure (see GET_DEFAULT_PARAMS)
%
%   OUTPUTS
%   y           :=  Array of processed signals (column vectors)
%
%   For info see: RESAMPLE

%   Mar. 2016 - Dario Pilori <dario.pilori@polito.it>

%% Retrieve parameters from stucture
fADC = Parameters.fADC;                                                    % ADC sampling rate (Hz)
Rs   = Parameters.Rs;                                                      % Symbol rate (Baud)

%% Check input parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(fADC,{'numeric'},{'scalar','positive'},'Parameters','fADC',2);
validateattributes(Rs,{'numeric'},{'scalar','positive'},'Parameters','Rs',2);

%% Resample
if fADC < 1.25*Rs
    warning('ADC sampling rate is lower than 1.25 SPS. Expect problems.');
end
[n1,n2] = rat(2*Rs/fADC);                                                  % calculate fractional resampling factor
if n1/n2~=1
    y = resample(x,n1,n2);                                                 % low-pass filter and resample (see RESAMPLE)
else
    y = x;
end