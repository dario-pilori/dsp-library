function performance_monitoring(w,a,Parameters)
%PERFORMANCE_MONITORING     Interpret complex-valued equalizer taps
%   This function interprets the adaptive equalizer taps w to detect PMD,
%   residual chromatic dispersion, and DGD.
%
%   It works only for complex-valued taps, otherwise it will erroneously
%   interpret transmitter/receiver imbalances as PMD.
%
%   INPUTS:
%   w           :=  Complex-valued equalizer taps
%   a           :=  Transmit sequence
%   Parameters  :=  Equalizer parameter structure (see GET_DEFAULT_PARAMS)
%
%   OUTPUTS

%   Apr. 2016 - Dario Pilori <dario.pilori@polito.it>

%% Retrieve parameters  
eqmode = Parameters.eqmode;                                                % real-valued or complex-valued equalizer
Ntaps = Parameters.Ntaps;                                                  % number of equalizer taps
C = Parameters.C;                                                          % get constellation

%% Check input parameters
validateattributes(w,{'double'},{'2d'},'','w',1);
validateattributes(a,{'numeric'},{'2d'},'','a');
validateattributes(C,{'numeric'},{'column'},'','C');
validateattributes(Ntaps,{'numeric'},{'integer'},'','Ntaps');
if ~strcmp(eqmode,'complex')
    error('performance_monitoring works only on complex-valued equalizers')
end

%% Process equalizer taps
Rs = 32e9;
f_ax = (-Ntaps/2:Ntaps/2-1)'/Ntaps*(2*Rs)*1e-9;
Ntx = size(a,2);                                                           % number of transmitted signals
w = reshape(w,Ntaps,[]);                                                   % reshape taps
Nrx = size(w,2)/Ntx;                                                       % number of received signals
W = fftshift(fft(w));                                                      % Fourier-transform them

%% Calculate PDL
PDL = NaN(Ntaps,1);
phi = NaN(Ntaps,1);
for i = 1:Ntaps                                                            % for each frequency
    s = svd(reshape(W(i,:),Ntx,Nrx));                                      % singular values of Jones matrix
    PDL(i) = 20*log10(abs(max(s)/min(s)));                                 % calculate PDL/MDL
    phi(i)=angle(sqrt(det(reshape(W(i,:),Ntx,Nrx))));                      % residual phase
end
phi = unwrap(phi);
CD = (phi+flip(phi))/2;
beta2 = polyfit(f_ax(1+Ntaps/4:end-Ntaps/4),CD(1+Ntaps/4:end-Ntaps/4)-max(CD),2)/2/pi^2*1e6;
disp(['Residual D=',num2str(-beta2(1)/1.274,'%3.2f'),' ps/nm'])

%% Show results
figure(1)
clf
plot(f_ax,PDL)
grid on
xlabel('f (GHz)')
ylabel('PDL (dB)')

figure(2)
clf
plot(f_ax,CD)
grid on
xlabel('f (GHz)')
ylabel('\phi (radians)')
