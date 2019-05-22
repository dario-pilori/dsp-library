function reqOSNR = osnr_penalty(OSNRs,BERs,tBER)
%OSNR_PENALTY   Calculates OSNR penalty
% This function calculates the required OSNR from a BER/OSNR curve using
% a linear interpolation in dB scale
% Nov. 2015 Dario Pilori <dario.pilori@polito.it>
% INPUTS: %
% OSNRs:    [Nosnr x 1] OSNR values
% BERs:     [Nosnr x M] BER values. Multiple curves are supported on the second dimension.
% tBER:     [1 x 1]     Target BER.
% OUTPUTS: %
% reqOSNR   [M x 1]     Vector of required OSNRs. NaN if error.

if size(BERs,1)~=size(OSNRs,1)
    error('Wrong dimension');
end
reqOSNR = NaN(size(BERs,2),1);
if tBER>0
    BERs(BERs>0.3) = NaN;
end
BERs(abs(BERs)==Inf) = NaN;
for n = 1:size(BERs,2)
    test_ber = BERs(:,n);
    if mean(diff(test_ber))<=0
        i_before = find(test_ber>tBER,1,'last');
        i_after = find(test_ber<tBER,1);
    else
        i_before = find(test_ber<tBER,1,'last');
        i_after = find(test_ber>tBER,1);
    end
    if isempty(i_after)||isempty(i_before)||...
            i_before>=size(test_ber,1)||i_after>size(test_ber,1)
        continue;
    end
    reqOSNR(n) = OSNRs(i_before)+(OSNRs(i_after)-OSNRs(i_before))*...
        (tBER-test_ber(i_before))/(test_ber(i_after)-test_ber(i_before));
end