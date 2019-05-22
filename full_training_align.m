function aTR = full_training_align(x,a)
%FULL_TRAINING_ALIGN     Align training sequence to input data
%   This function aligns the complex-valued received sequence
%   x to the transmit data a, at 1SPS, AFTER the equalizer.
%   This function operates independently on the two polarizations, and it
%   can handle phase conjugations due to modulator.
%
%   This function is useful after a fully blind equalizer (e.g. CMA).
%
%   INPUTS:
%   x           :=  Received sequences (1SPS, after equalizer)
%   a           :=  Transmit sequence
%
%   OUTPUTS
%   aTR         :=  Aligned transmit sequence

%   Apr. 2018 - Dario Pilori <dario.pilori@polito.it>

%% Check input parameters
validateattributes(x,{'numeric'},{'2d'},'','x',1);
validateattributes(a,{'numeric'},{'2d','ncols',size(x,2)},'','a',2);

%% Align
aTR = a;
for ll = 1:size(a,2)
    M = NaN(2*size(x,2),1);
    idx = NaN(2*size(x,2),1);
    for k = 1:size(x,2)
        [M(2*k),idx(2*k)] = max(abs(xcorr(x(1:length(a),k),a(:,ll))));     % regular correlation
        [M(2*k-1),idx(2*k-1)] = max(abs(xcorr(x(1:length(a),k),...
            conj(a(:,ll)))));                                              % correlation with conjugation
    end
    [~,best] = max(M);                                                     % get best correlation
    aTR(:,ll) = circshift(a(:,ceil(best/2)), idx(best)-length(a));         % align TX sequence
    if mod(best,2)
        aTR(:,ll) = conj(aTR(:,ll));                                       % conjugate if needed
    end
end