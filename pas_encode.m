function [sig,k,Beta] = pas_encode(bit_TX,k,m,n,gam,Ndiv,ps,hEnc)
%PAS_ENCODE     Encode bits using Probabilistic Amplitude Shaping (PAS)
% This function encodes a bit sequence bit_TX using a PAS encoder, that
% consists in probabilistic-shaping using the CCDM algorithm, followed by
% LDPC encoding. The encoding is done in a way to preserve the
% probabilistic shaping done by CCDM, by adding the parity-check bits as
% the sign of the constellation.

% References:
% - F. Buchali et al., "Rate Adaptation and Reach Increase by Probabilistically Shaped 64-QAM: An Experimental Demonstration", IEEE JLT 34.7, 1599-1609 (2016).
% - P. Schulte et al., "Constant Composition Distribution Matching", IEEE IT 62.1, 430-434 (2016)

%% Check input parameters
validateattributes(k,{'numeric'},{'scalar','integer','positive'})
validateattributes(m,{'numeric'},{'scalar','integer','positive'})
validateattributes(n,{'numeric'},{'scalar','integer','positive'})
validateattributes(gam,{'numeric'},{'scalar','>=',0,'<=',1})
validateattributes(Ndiv,{'numeric'},{'scalar','integer','positive'})
validateattributes(ps,{'struct'},{})
validateattributes(hEnc,{'comm.LDPCEncoder'},{})
validateattributes(bit_TX,{'logical'},{'size',[k+gam*n,2]});

%% For the two quadratures
sig = NaN(n,2);
for i = 1:2 
    %% Encode amplitudes
    a = NaN(n,1);                                                          % preallocate amplitudes
    for jj = 1:Ndiv
        a((jj-1)*n/Ndiv+(1:n/Ndiv)) = ccdm.encode(...
            bit_TX((jj-1)*k/Ndiv+(1:k/Ndiv),i),ps.n_ccdm).';               % Apply CCDM to generate amplitude levels
    end
    
    %% Encode sign
    Beta = circshift(de2bi(gray2bin(0:2^(m-1)-1,'pam',2^(m-1))),-1);       % beta lookup-table to convert amplitude to bits
    bit = Beta(a+1,:).';                                                   % encode amplitudes into bits
    bit_ENC = step(hEnc,[reshape(bit,n*(m-1),1);...
        bit_TX(k+1:end,i)]);                                               % LDPC encode
    sgn = sign([bit_TX(k+1:end,i);bit_ENC(end-(1-gam)*n+1:end)]-0.5);      % generate sign vector

    %% Generate m-ASK constellation
    sig(:,i) = (a*2+1).*sgn;                                               % modulate amplitude and sign
end

%% Generate m^2-QAM constellation
sig = sig(:,1)+1j*sig(:,2);