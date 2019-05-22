function [MI,GMI] = MI_eval_GH(C,B,snr,varargin)
%MI_EVAL_GH     Evaluation of mutual information using Gauss-Hermite quadrature
%   This function evaluates the mutual information (MI) and the generalized
%   mutual information (GMI) of the 2D (complex) constellation C over an
%   AWGN channel with standard deviation of the complex-valued noise equal
%   to sigma_n. Evaluation is performed using the Gauss-Hermite quadrature
%   which allows fast numerical integration of the mutual information. For
%   GMI calculation the bit mapping B is used.
%
%   Computations taken from A. Alvarado et al., "Achievable Information Rates
%   for Fiber Optics: Applications and Computations", J. Lightw. Technol. 36(2) pp. 424-439
%   https://dx.doi.org/10.1109/JLT.2017.2786351
%   Parameters of Gauss-Hermite taken from: http://keisan.casio.com/exec/system/1281195844
%
%   INPUTS:
%   C   :=      Constellation [M x 1]
%   B   :=      Bit mapping [M x log2(M)]
%   snr :=      Signal-to-noise ratios (Es/No, dB unit) [N x 1]
%   Pk  :=      Probability of each constellation point (optional) [M x 1]
%
%   OUTPUTS
%   MI      :=      Mutual information in bits [N x 1]
%   GMI     :=      Generalized mutual information in bits [N x 1]

%   February 2018 - Dario Pilori

%% Check input parameters
M = size(C,1);
validateattributes(C,{'numeric'},{'column'},'Parameters','C',1);
validateattributes(B,{'numeric','logical'},{'binary','nrows',M,'ncols',log2(M)},'','B',2);
validateattributes(snr,{'numeric'},{'column','nonnegative'},'','sigma_n',3);
if nargin>3
    Pk = varargin{1};
    validateattributes(Pk,{'numeric'},{'column','nonnegative','nrows',M},'','Pk',4);
else
    Pk = ones(M,1)/M; % equiprobable
end

%% Calculate sigma_n
sigma_n = sqrt(sum(Pk.*abs(C).^2))*10.^(-snr/20);

%% Params for Gauss-Hermite quadrature
gh = [-3.436159118837737603327	7.64043285523262062916E-6;
    -2.532731674232789796409	0.001343645746781232692202;
    -1.756683649299881773451	0.0338743944554810631362;
    -1.036610829789513654178	0.2401386110823146864165;
    -0.3429013272237046087892	0.6108626337353257987836;
    0.3429013272237046087892	0.6108626337353257987836;
    1.036610829789513654178	    0.2401386110823146864165;
    1.756683649299881773451	    0.03387439445548106313617;
    2.532731674232789796409	    0.001343645746781232692202;
    3.436159118837737603327	    7.64043285523262062916E-6];
x = gh(:,1);
w = gh(:,2);
clear gh

%% Evaluate Mutual Information
MI = zeros(size(sigma_n,1),1);
for i = 1:size(sigma_n,1)
    for l = 1:M
        for m = 1:size(x,1)
            MI(i) = MI(i)-Pk(l)/pi*w(m)*sum(...
                w.*log2(sum(Pk'.*exp(...
                -(abs(C(l)-C.').^2-2*sigma_n(i)*real((x+1j*x(m)).*(C(l)-C.')))...
                /sigma_n(i)^2),2)));
        end
    end
end

%% Evaluate Generalized Mutual Information
% To be optimized...
if nargout>1
    GMI = zeros(size(sigma_n,1),1);
    for z = 1:size(sigma_n,1)
        for k = 1:log2(M)
            for b = [false,true]
                I = C(B(:,k)==b);
                PI = Pk(B(:,k)==b);
                PIs = sum(PI);
                for i = 1:size(I,1)
                    for l = 1:size(x,1)
                        GMI(z) = GMI(z)-w(l)*PI(i)/pi*sum(w.*log2(...
                            sum(Pk'.*exp(-(abs(I(i)-C.').^2 - 2*sigma_n(z)*...
                            real((x+1j*x(l)).*(I(i)-C.')))/sigma_n(z)^2),2)./...
                            sum(PI'.*exp(-(abs(I(i)-I.').^2 - 2*sigma_n(z)*...
                            real((x+1j*x(l)).*(I(i)-I.')))/sigma_n(z)^2),2)*PIs));
                    end
                end
            end
        end
    end
    Pb = reshape(sum(Pk.*[~B B],1),[],2);
    GMI = GMI - Pk'*log2(Pk) + sum(sum(Pb.*log2(Pb)));
end