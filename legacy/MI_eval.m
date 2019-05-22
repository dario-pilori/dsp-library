function [MI,SNR,i_rx,GMI] = MI_eval(x,ax,Mc,C,B,dif)
%MI_EVAL     Evaluation of mutual information
%   This function evaluates the mutual information between  x and ax.
%   This function is meant to be run only by DECODING_AND_BER
%
%   INPUTS:
%   x   :=      Received symbols [Ns x Ntx]
%   ax  :=      Indexes of transmitted symbols [Ntx x Ns/pp] or [1 x Ns/pp]
%   Mc  :=      Bits/symbol of the constellation [scalar]
%   C   :=      Constellation [2^Mc x 1]
%   B   :=      Bit mapping [2^Mc x Mc]
%   dif :=      Differentially encoded symbols [scalar logical]
%
%   OUTPUTS
%   MI      :=      Mutual information in bits [Ntx x 1]
%   SNR     :=      Signal-to-noise ratio Es/No in dB [Ntx x 1]
%   i_rx    :=      Indices of received symbols [Ns x Ntx]
%   GMI     :=      Generalized mutual information in bits [Ntx x 1]

%   June 2016 - Gabriella Bosco and Dario Pilori

%% Check input parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(size(x,1)/size(ax,2),{'numeric'},{'integer'},'','length(x)/length(ax)',1);
validateattributes(ax,{'double'},{'2d'},'','ax',2);
validateattributes(Mc,{'numeric'},{'scalar','positive','integer'},'Parameters','Mconst',3);
validateattributes(C,{'numeric'},{'column','nrows',2^Mc},'Parameters','C',4);
validateattributes(B,{'numeric','logical'},{'binary','nrows',2^Mc,'ncols',Mc},'','B',5);
validateattributes(dif,{'logical'},{'scalar'},'','differential',6);
debug = false; % set to true to enable some debugging features

%% Initialize
Ns = size(x,1);                                                            % number of received symbols
Ntx = size(x,2);                                                           % number of received signals
M = 2^Mc;                                                                  % number of points in the constellation
centroids = NaN(M,1);                                                      % centroids of constellation
sigma2 = NaN(M,1);                                                         % variance of each centroid
MI = NaN(Ntx,1);                                                           % mutual information
GMI = NaN(Ntx,1);                                                          % generalized mutual information
SNR = NaN(Ntx,1);                                                          % SNR Es/No (dB)
i_rx = NaN(Ns,Ntx);                                                        % indices of received symbols

%% For all received signals
for nn = 1:Ntx
    %% Calculate centroids, variance and MI
    y = x(:,nn);                                                           % get nn-th transmit signal
    i_tx = repmat(ax(mod(nn-1,size(ax,1))+1,:),1,Ns/size(ax,2));           % nn-th transmit indexes
    pos = bsxfun(@eq,int16(i_tx),(1:M).');                                 % indexes of this point
    P_symb = sum(pos,2)/Ns;                                                % probability of point
    if any(P_symb<=eps)
        warning('The transmit sequence does not contain all symbols. Please use a longer sequence.');
        return;
    end
    for i = 1:M                                                            % for all transmitted symbols
        rx_s = y(pos(i,:));                                                % get received points relative to i-th TX symbol
        centroids(i) = mean(rx_s);                                         % calculate centroid
        sigma2(i) = var(rx_s)+eps;                                         % calculate variance (never to zero)
    end
    SNR(nn) = 10*log10(sum(P_symb.*abs(centroids).^2)./sum(P_symb.*sigma2));             % calculate SNR
    
    %% Calculate MI
    if ~dif
        MI_den = sum(bsxfun(@times,exp(bsxfun(@rdivide,-abs(bsxfun(@minus,y,centroids.')).^2,sigma2.')),P_symb.'),2);% cumulate MI
        MI(nn) = -mean(log2(MI_den))-1/log(2);                             % convert to bits
    end
    
    %% Calculate generalized MI (if needed)
    if (nargout>=4)&&(~dif)
        %% Bit probabilities
        P_symb_bit = bsxfun(@times,P_symb,[~B B]);
        P_bit = reshape(sum(P_symb_bit,1),[],2).';
        P_symb_bit = bsxfun(@rdivide,P_symb_bit,sum(P_symb_bit,1));
        
        %% GMI
        GMI_num = zeros(Mc,Ns);
        for nsymb=1:M
            P_bit_nsymb = P_symb_bit(nsymb,:);
            GMI_num = GMI_num + bsxfun(@times, exp(-abs(y-centroids(nsymb)).^2/sigma2(nsymb)), P_bit_nsymb(bsxfun(@plus, B(i_tx,:)*Mc, (1:Mc)))).';
        end
        GMI(nn) = sum(mean(log2(GMI_num),2)-mean(log2(MI_den)));
        
        %% Correction factor for statistical shaping
        H_symb = sum(P_symb.*log2(1./P_symb));                             % entropy of each symbol
        H_bit = sum(sum(P_bit.*log2(1./P_bit)));                           % sum of entropies of each bit
        corr = H_bit-H_symb;                                               % correction factor
        GMI(nn) = GMI(nn)-corr;                                            % apply correction factor
    end
    
    %% Minimum distance decoding if needed
    if nargout>=3
        [~,i_rx(:,nn)] = min(abs(repmat(y,1,M)-repmat(centroids.',Ns,1)),[],2);
    end
    
    %% Print useful debugging stuff
    if debug
        figure
        plot(y,'.')
        axis equal
        grid on
        hold on
        x_R = bsxfun(@plus,real(centroids),sigma2*[-0.5,0.5]);
        x_I = bsxfun(@plus,imag(centroids),sigma2*[-0.5,0.5]);
        for i = 1:M
            rectangle('Position',[x_R(i,1) x_I(i,1) x_R(i,2)-x_R(i,1) x_I(i,2)-x_I(i,1)],'Curvature',[1 1])
        end
        plot(centroids,'.','MarkerSize',18)
        hold off
    end
end

%% Correct SNR
if isreal(C)
    SNR = SNR-10*log10(2);
end