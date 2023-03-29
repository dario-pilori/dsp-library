function [y,w,e] = adaptive_equalizer(x,a,Parameters)
%ADAPTIVE_EQUALIZER     Adaptive equalization at 2SPS
%   This function applies a time-domain adaptive equalizer to the received
%   signal at 2SPS. The equalizer can either be CMA-based or LMS-based,
%   real-valued or complex-valued. In case of data-aided transmission, 
%   an aligned training sequence is required (see POWER_TRAINING_ALIGN).
%
%   INPUTS:
%   x           :=  Array of signals (column vectors) to process at 2SPS
%   a           :=  Transmit sequence
%   Parameters  :=  Equalizer parameters sub-structure (see GET_DEFAULT_PARAMS)
%
%   OUTPUTS
%   y           :=  Array of processed signals (column vectors) at 2SPS
%   w           :=  Equalizer taps [Nrx*Nd*Ntaps x Ntx*Nd]
%   e           :=  Equalizer error

%   Mar. 2016 - Dario Pilori <dario.pilori@polito.it>

%% Retrieve parameters  
Ntaps = Parameters.Ntaps;                                                  % taps of the adaptive equalizer
eqmode = Parameters.eqmode;                                                % real-valued or complex-valued equalizer
mus = Parameters.mus;                                                      % adaptation coefficents
methods = Parameters.methods;                                              % equalizer methodds
Ks = Parameters.Ks;                                                        % number of symbols to run the equalizer
M = Parameters.Mlms;                                                       % memory of LMS phase recovery
lambda = Parameters.lambda;                                                % leakage factor for CMA
C = Parameters.C;                                                          % constellation
if isfield(Parameters,'Mpam')
    Mpam = Parameters.Mpam;                                                % if PAM, then get number of levels
    validateattributes(log2(Mpam),{'numeric'},{'scalar','positive','integer'},'Parameters','Mpam',3);
end

%% Check input parameters
validateattributes(x,{'double'},{'2d'},'','x',1);
validateattributes(Ntaps,{'numeric'},{'scalar','positive','integer','even'},'Parameters','Ntaps',2);
validateattributes(M,{'numeric'},{'scalar','>=',1,'integer'},'Parameters','Mlms',2);
validateattributes(methods,{'cell'},{'column'},'Parameters','methods',2);
validateattributes(mus,{'numeric'},{'nonnegative','size',[size(methods,1) 1]},'Parameters','mus',2);
validateattributes(Ks,{'numeric'},{'nonnegative','size',[size(methods,1) 1]},'Parameters','Ks',2);
validateattributes(lambda,{'numeric'},{'scalar','<',1,'>=',0},'Parameters','lambda',2);
validateattributes(a,{'numeric'},{'2d'},'','a');
validateattributes(C,{'numeric'},{'column'},'','C');

switch eqmode
    case 'complex'
        Nd = 1;
    case 'real'
        Nd = 2;
        x = complex_to_real(x);
    otherwise
        error('eqmode must be ether complex or real.');
end

%% Precalculate parameters
Lx = size(x,1);                                                            % number of samples to process
Lpk = floor(Ntaps/2);                                                      % peak position in equalizer taps
Nrx = size(x,2)/Nd;                                                        % number of received signals
Ntx = size(a,2);                                                           % number of transmitted signals

%% Initialize
y = zeros(Lx,Ntx*Nd);                                                      % output signal
e = zeros(Lx,Ntx*Nd);                                                      % error signal
fin = 0;                                                                   % initialize end index
if any(strcmp('mma-dd',methods))
    Rdd = repmat(unique(round(abs(C).^2,8)),1,Ntx);                        % MMA-DD radii
end
if any(strcmp('cma',methods))
    Rcma = repmat(mean(abs(C).^2),1,Ntx);                                  % CMA single radius
end

%% Initialize equalizer taps
w = zeros(Nrx*Nd*Ntaps,Ntx*Nd);                                            % equalizer taps
if ~(strcmp('lms',methods{1})||strcmp('r-lms-dd',methods{1}))              % LMS has to be initialized to zero
    for nn = 1:Ntx*Nd
        w(Lpk+(nn-1)*Ntaps,nn) = 1;                                        % initialize as an impulse (CMA)
    end
end

%% Run equalizers
for nn = 1:size(methods,1)                                                 % for all equalizers
    %% Inizialize parameters
    iniz = max(Ntaps,fin+1);                                               % starting point of equalizer
    fin = min(Lx-Ntaps-1,Ks(nn));                                          % ending point of equalizer
    
    %% Find correct equalizer and run it
    if Nd==1                                                               % Complex-valued equalization
        switch methods{nn}
            case 'mma-da'
                [y,w,e] = eq_mmada_complex(x,y,a,w,e,mus(nn),iniz,fin,Lpk,Ntaps,Nrx,lambda);
            case 'mma-dd'
                [y,w,e] = eq_cma_complex(x,y,w,e,mus(nn),iniz,fin,Rdd,Lpk,Ntaps,Nrx,lambda);
            case 'lms'
                [y,w,e] = eq_lms_complex(x,y,a,w,e,mus(nn),iniz,fin,M,Lpk,Ntaps,Ntx,Nrx);
            case 'r-lms-dd'
                [y,w,e] = eq_lms_dd_complex(x,y,w,e,mus(nn),iniz,fin,Lpk,Ntaps,Nrx,Mpam);
            case 'cma'
                [y,w,e] = eq_cma_complex(x,y,w,e,mus(nn),iniz,fin,Rcma,Lpk,Ntaps,Nrx,lambda);
            otherwise
                error(['Complex-valued equalizer method ',methods{nn},' not supported.']);
        end
    elseif Nd==2                                                           % Real-valued equalization
        switch methods{nn}
            case 'mma-da'
                [y,w,e] = eq_mmda_real(x,y,a,w,e,mus(nn),iniz,fin,Lpk,Ntaps,Ntx,Nrx,lambda);
            case 'mma-dd'
                [y,w,e] = eq_cma_real(x,y,w,e,mus(nn),iniz,fin,Rdd,Lpk,Ntaps,Ntx,Nrx,lambda);                
            case 'lms'
                [y,w,e] = eq_lms_real(x,y,a,w,e,mus(nn),iniz,fin,M,Lpk,Ntaps,Ntx,Nrx);
            case 'cma'
                [y,w,e] = eq_cma_real(x,y,w,e,mus(nn),iniz,fin,Rcma,Lpk,Ntaps,Ntx,Nrx,lambda);
            otherwise
                error(['Real-valued equalizer method ',methods{nn},' not supported.']);
        end
    end
end

%% End
y = y(Ntaps:end-Ntaps-1,:);                                                % cut head and tail
if Nd==2
    y = real_to_complex(y);                                                % switch to complex numbers
end
end
% End of adaptive equalizer

% Start of helper functions
%% Complex-valued (e.g. 2x2) data-aided MMA
function [y,w,e] = eq_mmada_complex(x,y,a,w,e,mu,iniz,fin,Lpk,Ntaps,Nrx,lambda)%#codegen
%% Initialize
vec = -Lpk:1:Lpk-1;                                                        % preallocate data vector

%% For all samples
for nn = iniz:fin
    %% Apply equalizer
    u = reshape(x(nn+vec,:),Nrx*Ntaps,1);                                  % get signal
    y(nn,:)= u'*w;                                                         % apply equalizer
    
    %% Calculate error on even samples
    if rem(nn,2)==0
        i_a = mod(nn/2-1,size(a,1))+1;                                     % index in training sequence
        e(nn,:) = abs(y(nn,:)).^2-abs(a(i_a,:)).^2;                        % calculate error
        w = (1-lambda)*w-mu*u*(e(nn,:).*y(nn,:));                          % apply equalizer
    end
end
end

%% Real-valued (e.g. 4x4) data-aided MMA
function [y,w,e] = eq_mmda_real(x,y,a,w,e,mu,iniz,fin,Lpk,Ntaps,Ntx,Nrx,lambda)%#codegen
%% Initialize
vec = -Lpk:1:Lpk-1;                                                        % preallocate data vector

%% Run
for nn = iniz:fin
    %% Apply equalizer
    u = reshape(x(nn+vec,:),Nrx*2*Ntaps,1);                                % get signal
    y(nn,:)= u.'*w;                                                        % apply equalizer
    
    %% Calculate error on even samples
    if rem(nn,2)==0
        i_a = mod(nn/2-1,size(a,1))+1;                                     % index in training sequence
        e(nn,:) = reshape(repmat(y(nn,1:2:Ntx*2).^2+y(nn,2:2:Ntx*2).^2-...
            abs(a(i_a,:)).^2,2,1),1,Ntx*2);                                % calculate error
        w = (1-lambda)*w-mu*u*(e(nn,:).*y(nn,:));                          % update equalizer
    end
end
end

%% Complex-valued (e.g. 2x2) decision-directed CMA
function [y,w,e] = eq_cma_complex(x,y,w,e,mu,iniz,fin,Rdd,Lpk,Ntaps,Nrx,lambda)%#codegen
%% Initialize
vec = -Lpk:1:Lpk-1;                                                        % preallocate data vector

%% Run
for nn = iniz:fin
    %% Apply equalizer
    u = reshape(x(nn+vec,:),Nrx*Ntaps,1);                                  % get signal
    y(nn,:)= u'*w;                                                         % apply equalizer
    
    %% Calculate error on even samples
    if rem(nn,2)==0
        [~,idx]=min(abs(Rdd...
            -repmat(abs(y(nn,:)).^2,size(Rdd,1),1)),[],1);                 % get best radius
        e(nn,:) = abs(y(nn,:)).^2-Rdd(idx);                                % calculate error
    end
    
    %% Update equalizer
    w = (1-lambda)*w-mu*u*(e(nn,:).*y(nn,:));
end
end

%% Real-valued (e.g. 4x4) decision-directed CMA
function [y,w,e] = eq_cma_real(x,y,w,e,mu,iniz,fin,Rdd,Lpk,Ntaps,Ntx,Nrx,lambda)%#codegen
%% Initialize
vec = -Lpk:1:Lpk-1;                                                        % preallocate data vector

%% Run
for nn = iniz:fin
    %% Apply equalizer
    u = reshape(x(nn+vec,:),Nrx*2*Ntaps,1);                                % get signal
    y(nn,:)= u.'*w;                                                        % apply equalizer
    
    %% Calculate error on even samples
    if rem(nn,2)==0
        [~,idx]=min(abs(Rdd-repmat(y(nn,1:2:Ntx*2).^2+y(nn,2:2:Ntx*2).^2,size(Rdd,1),1)),[],1);% get best radius
        e(nn,:) = reshape(repmat(y(nn,1:2:Ntx*2).^2+y(nn,2:2:Ntx*2).^2-Rdd(idx)...
            ,2,1),1,Ntx*2);                                                % calculate error
        w = (1-lambda)*w-mu*u*(e(nn,:).*y(nn,:));                          % update equalizer
    end
end
end

%% Complex-valued (e.g. 2x2) data-aided LMS
function [y,w,e] = eq_lms_complex(x,y,a,w,e,mu,iniz,fin,M,Lpk,Ntaps,Ntx,Nrx)%#codegen
%% Initialize
vec = -Lpk:1:Lpk-1;                                                        % preallocate data vector
ph = ones(M,Ntx);                                                          % vector of phase errors

%% Run
for nn = iniz:fin
    %% Apply equalizer
    u = reshape(x(nn+vec,:),Nrx*Ntaps,1);                                  % get signal
    y(nn,:)= w'*u;                                                         % apply equalizer
    
    %% Calculate error on even samples
    if rem(nn,2)==0
        i_a = mod(nn/2-1,size(a,1))+1;
        if M>1
            ph = circshift(ph,[1,0]);
            ph(1,:) = exp(-1j*(atan2(imag(y(nn,:)),real(y(nn,:)))-atan2(imag(a(i_a,:)),real(a(i_a,:)))));
        end
        e(nn,:) = y(nn,:) - a(i_a,:).*conj(sign(sum(ph,1)));
        w = w-mu*u*conj(e(nn,:));                                          % apply equalizer
    end
end
end

%% Real-valued (e.g. 4x4) data-aided LMS
function [y,w,e] = eq_lms_real(x,y,a,w,e,mu,iniz,fin,M,Lpk,Ntaps,Ntx,Nrx)%#codegen
%% Initialize
vec = -Lpk:1:Lpk-1;                                                        % preallocate data vector
ph = ones(M,Ntx);                                                          % vector of phase errors

%% Run
for nn = iniz:fin
    %% Apply equalizer
    u = reshape(x(nn+vec,:),2*Nrx*Ntaps,1);                                % get signal
    y(nn,:)= w'*u;                                                         % apply equalizer
    
    %% Calculate error on even samples
    if rem(nn,2)==0
        i_a = mod(nn/2-1,size(a,1))+1;                                     % index of training symbol
        if M>1                                                             % if phase recovery
            ph = circshift(ph,[1,0]);                                           % shift phase moving average
            ph(1,:) = exp(-1j*(atan2(y(nn,2:2:end),y(nn,1:2:end))-atan2(imag(a(i_a,:)),real(a(i_a,:)))));% calculate phase error
        end
        e(nn,:) = y(nn,:) - complex_to_real(a(i_a,:).*conj(sign(sum(ph,1))));% calculate equalizer error
        w = w-mu*u*e(nn,:);                                                % apply equalizer
    end
end
end

%% Complex-valued (e.g. 2x2) decision-directed LMS
% This function works ONLY with real-valued modulation formats (e.g. PAM)
function [y,w,e] = eq_lms_dd_complex(x,y,w,e,mu,iniz,fin,Lpk,Ntaps,Nrx,Mpam)%#codegen
%% Initialize
vec = -Lpk:1:Lpk-1;                                                        % preallocate data vector

%% Run
for nn = iniz:fin
    %% Apply equalizer
    u = reshape(x(nn+vec,:),Nrx*Ntaps,1);                                  % get signal
    y(nn,:)= w'*u;                                                         % apply equalizer
    
    %% Calculate error on even samples
    if rem(nn,2)==0
        e(nn,:) = y(nn,:) - pam_hard_decision(y(nn,:),Mpam);
        w = w-mu*u*conj(e(nn,:));                                          % apply equalizer
    end
end
end

%% PAM hard decision
function d = pam_hard_decision(x,Mpam)
switch Mpam    
    case 2
        d = sign(x);
    case 4
        d = 2*sign(x);
        d = d+sign(x-d);
    case 8
        d = 4*sign(x);
        d = d+2*sign(x-d);
        d = d+sign(x-d);
    case 16
        d = 8*sign(x);
        d = d+4*sign(x-d);
        d = d+2*sign(x-d);
        d = d+sign(x-d);
    otherwise
        error('Number of PAM levels not supported.');
end
end