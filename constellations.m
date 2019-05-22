function [C,B] = constellations(M,varargin)
%CONSTELLATIONS     Constellation and bit mapping
%   This function returns the constellation and bit mapping for a given
%   number of bits per symbol M.
%
%   INPUTS:
%   M       :=  Number of bits per symbol
%   OUTPUTS:
%   C       :=  Complex-valued constellation (2^M x 1)
%   B       :=  Bit-mapping (2^M x M)

%   Mar. 2016 - Dario Pilori <dario.pilori@polito.it>

%% Load probabilistic shaping
if nargin>=2&&varargin{1}
    % Load geometrically-shaped constellation
    fname = ['~/Documents/MATLAB/dsp-library/optimized_constellations/',num2str(2^M),'-QAM_optimized.mat'];
    if ~exist(fname,'file')
        warning('Optimized constellation NOT found, using standard constellation...');
    else
        load(fname,'C','B');
        return
    end
end

%% Load text constellation
cur_dir = fileparts(mfilename('fullpath'));
if M==0
    C = 0;
    B = 0;
else
    switch M
        case -1
            filestr = 'text_constellations/nrz.txt';
        case 1
            filestr = 'text_constellations/bpsk.txt';
        case 2
            filestr = 'text_constellations/qpsk.txt';
        otherwise
            filestr = sprintf('text_constellations/%dqam.txt', 2^M);
    end
    
    C = load(fullfile(cur_dir, filestr));
    C = C(:,1) + 1j*C(:,2);
    B = de2bi(0:(2^abs(M)-1),'left-msb');
end