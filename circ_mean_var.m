function [ dir_mean, dir_var ] = circ_mean_var( theta_array, mag_array )
% mean & std in circular space, 180deg range

theta_array = reshape(theta_array,1,length(theta_array));
if nargin < 2
    mag_array = ones(size(theta_array));
end

Mdir = zeros(1,3);
Mdir(1) = mag_array*cosd(2*theta_array');
Mdir(2) = mag_array*sind(2*theta_array');
Mdir(3) = sqrt(Mdir(2)^2 + Mdir(1)^2);
dir_var = 1 - Mdir(3)/sum(mag_array);
dir_mean = acosd(Mdir(1)/Mdir(3));
if Mdir(2) < 0
    dir_mean = 360 - dir_mean;
end
dir_mean = dir_mean/2;
if dir_mean >= 180
    dir_mean = dir_mean - 180;
elseif dir_mean < 0
    dir_mean = dir_mean + 180;
end

end

