function bit_RX = pas_decode(x,m,n,k,gam,Ndiv,Beta,ps,C,Psymb,sigma2,hDec)
%PAS_DECODE     Decode signal encoded with PAS encoder
% This function dencodes a bit sequence encoded using a PAS encoder. It
% decodes the sequence using a LDPC decoder, then it recovers the bits by
% inverting the CCDM probabilistic-shaping algorithm.

% References:
% - F. Buchali et al., "Rate Adaptation and Reach Increase by Probabilistically Shaped 64-QAM: An Experimental Demonstration", IEEE JLT 34.7, 1599-1609 (2016).
% - P. Schulte et al., "Constant Composition Distribution Matching", IEEE IT 62.1, 430-434 (2016)

%% Check input parameters
validateattributes(k,{'numeric'},{'scalar','integer','positive'})
validateattributes(m,{'numeric'},{'scalar','integer','positive'})
M = 2^(2*m);
validateattributes(n,{'numeric'},{'scalar','integer','positive'})
validateattributes(gam,{'numeric'},{'scalar','>=',0,'<=',1})
validateattributes(Ndiv,{'numeric'},{'scalar','integer','positive'})
validateattributes(Beta,{'numeric'},{'size',[2^(m-1) m-1]})
validateattributes(ps,{'struct'},{})
validateattributes(C,{'numeric'},{'column','nrows',M})
validateattributes(Psymb,{'numeric'},{'column','nrows',M})
validateattributes(x,{'numeric'},{'column','nrows',n});

%% Calculate log-likelihood ratios
LLR = flipud(reshape(qam_llr_mex(C, sigma2*ones(size(C)), x, Psymb),...
    2*m,n));                                                               % calculate log-likelihood ratios

%% For all the quadratures
bit_RX = false(round(k+gam*n),2);                                          % preallocate receiver bits
for i = 1:2                                                                % for all the quadratures                                               
    q = (2-i)*m;                                                           % select i-th quadrature
    
    %% LDPC decode
    bit_DEC = step(hDec,[reshape(LLR(2+q:q+m,:),n*(m-1),1);reshape(LLR(1+q,:),n,1)]);% LDPC decoding
    
    %% Decode using CCDM
    [~,Locb] = ismember(reshape(bit_DEC(1:end-round(gam*n)),m-1,n).',Beta,'rows');% invert Beta
    ccdm_out_bits = NaN(k,1);
    for jj = 1:Ndiv
        ccdm_out_bits((jj-1)*k/Ndiv+(1:k/Ndiv)) = ...
            ccdm.decode(Locb((jj-1)*n/Ndiv+(1:n/Ndiv))-1,ps.n_ccdm,k/Ndiv).';   % invert CCDM
    end
    
    %% Build output vector
    bit_RX(:,i) = logical([ccdm_out_bits;bit_DEC(end-round(gam*n)+1:end)]); % construct output vector
end