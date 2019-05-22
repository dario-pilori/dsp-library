function [NCI,autocorr_phase,autocorr_abs] = eval_nlpn(cst,a,varargin)
%EVAL_NLPN  Evaluate statistics of non-linear phase noise
%  This function evaluates the statistics (e.g. variance and
%  autocorrelation) of the non-linear phase noise which affects the
%  constellation cst. In order to have a correct evaluation, cst should not
%  have any added ASE noise, nor Lorentzian laser phase noise.
%  For details see: 10.1109/JLT.2016.2613893
%
%   INPUTS:
%   cst             :=  Constellations over which NLPN should be evaluated
%   a               :=  Transmitted constellation
%   cut             :=  Initial samples to cut (optional)
%
%   OUTPUTS
%   NCI             :=  Non-circularity index (dB)
%   autocorr_phase  :=  Autocorrelation of phase noise (normalized)
%   autocorr_abs    :=  Autocorrelation of amplitude noise (normalized)

%   Mar. 2017 - Dario Pilori <dario.pilori@polito.it>

%% Check input parameters
validateattributes(cst,{'numeric'},{'2d'},'','cst',1);
validateattributes(a,{'numeric'},{'ncols',size(cst,2)},'','a',2);
if nargin>2
    cut = varargin{1};
    validateattributes(cut,{'numeric'},{'scalar','integer','nonnegative'},'','cut',3);
else
    cut = 0;
end

%% Align sequences
cst = cst(1+cut:end,:);
del = finddelay(a,cst(1:size(a,1),:));
a = circshift(a,mean(del));
clear del

%% Repeat transmit pattern and normalize
tx = repmat(a,ceil(size(cst,1)/size(a,1)),1);
tx = tx(1:size(cst,1),:);
cst = cst./tx;
clear tx

%% Calculate NCI
NCI = 10*log10(var(imag(cst))./var(real(cst)));

%% Calculate autocorrelations
if nargout>1
    autocorr_abs = NaN(size(a,1)+1,size(cst,2));
    autocorr_phase = NaN(size(a,1)+1,size(cst,2));
    for pol = 1:size(cst,2)
        autocorr_abs(:,pol) = abs(xcorr(real(cst(:,pol))-mean(real(cst(:,pol))),size(a,1)/2,'unbiased'));
        autocorr_phase(:,pol) = abs(xcorr(imag(cst(:,pol))-mean(imag(cst(:,pol))),size(a,1)/2,'unbiased'));
    end
    autocorr_abs = autocorr_abs./max(autocorr_abs);
    autocorr_phase = autocorr_phase./max(autocorr_phase);
end