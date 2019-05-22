function b = index_to_binary(iTX)
%INDEX_TO_BINARY Convert decimal indeces to binary
% This function converts a column vector of decimal indeces to a binary
% vector. This function is similar to MATLAB's de2bi.
%
% INPUTS
% iTX   :=  Vector to convert to binary [N x 1]
%
% OUTPUTS
% b     :=  Binary output [N x N_bit]

% March 2019 - Dario Pilori <dario.pilori@polito.it>

%% Check input
validateattributes(iTX,{'numeric'},{'column','integer','nonnegative'},'','iTX',1);

%% Convert
N_bit = ceil(log2(max(iTX)));
b = NaN(size(iTX,1),N_bit);
for i = 1:N_bit
    tmp = fix(iTX/2);
    b(:,i) = iTX - 2*tmp;
    iTX = tmp;
end
end