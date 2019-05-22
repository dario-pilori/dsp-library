function y = normalize_DC_AGC(x,Parameters)
%NORMALIZE_DC_AGC     Remove DC and normalize power
%   Ths function removes the residual DC and gain from the input signals,
%   treated separately.
%
%   It requires movmean and movstd, therefore MATLAB R2016a+ is needed
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process
%   Parameters  :=  Parameter structure (see GET_DEFAULT_PARAMS)
%
%   OUTPUTS
%   y           :=  Array of processed signals (column vectors)

%   Mar. 2016 - Dario Pilori <dario.pilori@polito.it>

%% Retrieve parameters from stucture
kDC = Parameters.kDC;                                                      % memory of DC-block (samples)
kAGC = Parameters.kAGC;                                                    % memory of AGC (samples)

%% Check input parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(kDC,{'numeric'},{'scalar','positive','integer'},'Parameters','kDC',2);
validateattributes(kAGC,{'numeric'},{'scalar','positive','integer'},'Parameters','kAGC',2);

%% Normalize DC and AGC
x = x-movmean(x,kDC,1);                                                    % remove DC
y = x./movstd(x,kAGC,1);                                                   % AGC