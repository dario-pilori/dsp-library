function d = harddec(x,Mc)
%HARDDEC     Square QAM hard decision
%   This function applies hard-decision to the received signal x.
%
%   INPUTS:
%   x           :=  Array of symbols (column vectors)
%   Mc          :=  Constellation bits/symb
%
%   OUTPUTS
%   d           :=  Hard-decision symbols

%   Apr. 2016 - Dario Pilori <dario.pilori@polito.it>

switch Mc
    case 2
        d = sign(real(x))+1j*sign(imag(x));
    case 4
        d = 2*(sign(real(x))+1j*sign(imag(x)));
        d = d+sign(real(x-d))+1j*sign(imag(x-d));
    case 6
        d = 4*(sign(real(x))+1j*sign(imag(x)));
        d = d+2*(sign(real(x-d))+1j*sign(imag(x-d)));
        d = d+sign(real(x-d))+1j*sign(imag(x-d));
    case 8
        d = 8*(sign(real(x))+1j*sign(imag(x)));
        d = d+4*(sign(real(x-d))+1j*sign(imag(x-d)));
        d = d+2*(sign(real(x-d))+1j*sign(imag(x-d)));
        d = d+sign(real(x-d))+1j*sign(imag(x-d));
    otherwise
        error('Constellation not supported');
end
end