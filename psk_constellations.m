function [C,B] = psk_constellations(M)
%PSK_CONSTELLATIONS     PSK constellation and bit mapping
%   This function returns the PSK constellation and bit mapping for a given
%   number of bits per symbol M.
%
%   INPUTS:
%   M       :=  Number of bits per symbol
%   OUTPUTS:
%   C       :=  Complex-valued constellation (2^M x 1)
%   B       :=  Bit-mapping (2^M x M)

%   Dec 2017 - Dario Pilori <dario.pilori@polito.it>
validateattributes(M,{'numeric'},{'scalar','>=',1,'integer'},'Parameters','M',1);

Mpsk = 2^M;
C = exp(1j*2*pi*(0:Mpsk-1)'/Mpsk);
[~,idx] = sort(bin2gray(0:Mpsk-1,'psk',Mpsk));
C = C(idx);
B = de2bi(0:Mpsk-1,'left-msb');