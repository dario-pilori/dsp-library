function y = delay_compensation(x,Parameters)
%DELAY_COMPENSATION     Compensate for the skew between channels
%   Ths function removes the skew between the input signals using a cubic
%   interpolation
%
%   It requires interp1, introduced before MATLAB R2006a
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process
%   Parameters  :=  Parameter structure (see GET_DEFAULT_PARAMS)
%
%   OUTPUTS
%   y           :=  Array of processed signals (column vectors)

%   Mar. 2016 - Dario Pilori <dario.pilori@polito.it>

%% Retrieve parameters from stucture
dt = Parameters.skew;
fs = Parameters.fADC;

%% Check input parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(dt,{'numeric'},{'size',[size(x,2),1],'nonnan'},'Parameters','skew',2);
validateattributes(fs,{'numeric'},{'scalar','positive'},'Parameters','fADC',2);

%% Compensate for skew
y = NaN(size(x));
for ii = 1:size(x,2)
    if dt(ii)==0
        y(:,ii) = x(:,ii);
    else
        y(:,ii) = interp1((0:size(x,1)-1).'/fs, x(:,ii), (0:size(x,1)-1).'/fs+dt(ii), 'pchip');
    end
end