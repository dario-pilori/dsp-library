function ber = ber_theory_bits(snr,b)
%BER_THEORY_BITS Theoretical BER
% this function finds the theoretical BER for a given SNR (Es/N0) in dB for the
% modulation format with b bits/symbol
% For square-QAM constellation, computation is exact
% For cross-QAM constellations, see https://dx.doi.org/10.1109/TWC.2005.857997
%
% INPUTS
% snr   :=  Es/No in dB (column vector)
% b     :=  Bits per constellation symbol (integer)
%
% OUTPUTS
% ber   :=  bit error rate

% June 2016 - Dario Pilori <dario.pilori@polito.it>

%% Check input parameters
validateattributes(snr,{'numeric'},{'column'},'','snr',1);
validateattributes(b,{'numeric'},{'scalar','integer','positive'},'','b',2);

%% Compute
C = constellations(b);
M = size(C,1);
if mod(log2(M),2)==0                                                       % square constellation
    ber = berawgn(10*log10(10.^(snr/10)/log2(M)),'qam',M);
else                                                                       % cross constellation
    if M>=8
        ber = (1+1/sqrt(2*M)+1/3/M)*(4-6/sqrt(2*M))/log2(M)/2*erfc(sqrt(48/(31*M-32)*10.^(snr/10)));% approx formula from Vitthaladevuni et al., Exact BER computation for cross QAM constellations
    else
        error('Constellation not supported');
    end
end