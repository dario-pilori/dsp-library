function [out,xc] = training_align(x,a,cut,power,align_tx)
%TRAINING_ALIGN     Align training sequence to input data
%   This function aligns the complex-valued received sequence
%   x to the transmit data a, at 2SPS, before the equalizer.
%   This function averages over all the received waveforms (e.g.
%   polarization).
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process at 2SPS
%   a           :=  Transmit sequence
%   cut         :=  Cut instead of circshift the sequence (logical)
%   power       :=  Align power of signals instead of full field (logical)
%   align_tx    :=  Align transmit sequence instead of received sequence (logical)
%
%   OUTPUTS
%   out         :=  Aligned signal
%   xc          :=  Averaged cross-correlation (useful for debug)

%   Apr. 2017 - Dario Pilori <dario.pilori@polito.it>

%% Check input parameters
validateattributes(cut,{'logical'},{'scalar'},'','cut',3);
validateattributes(power,{'logical'},{'scalar'},'','power',4);
validateattributes(align_tx,{'logical'},{'scalar'},'','align_tx',5);
if cut&&align_tx
    error('Cutting the transmit sequence does not seem a good idea');
end

%% Measure delay
if power
    [t0,xc] = find_delay_pm_2sps(abs(x),abs(a));
else
    [t0,xc] = find_delay_pm_2sps(x,a);
end

%% Compensate for delay
if cut
    out = x(2*size(a,1)-t0:end,:);                                         % cut head
else
    if align_tx
        out = circshift(a,-round(t0/2));                                   % shift TX data
    else
        out = circshift(x,t0-2*size(a,1));                                 % shift RX data
    end
end