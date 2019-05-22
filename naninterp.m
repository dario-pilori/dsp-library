function x = naninterp(x)
%NANINTERP  'linear' interpolation of points with NaN
% Taken from: http://it.mathworks.com/matlabcentral/fileexchange/8225-naninterp
x(isnan(x)) = interp1(find(~isnan(x)),x(~isnan(x)),find(isnan(x)),'linear'); 