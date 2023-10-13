function [ dir_mean, dir_var ] = color_mean_var( theta_array, mag_array )
% mean & std in circular space, 360deg range

theta_array = reshape(theta_array,1,length(theta_array));
if nargin < 2
    mag_array = ones(size(theta_array));
end

Mdir = zeros(1,3);
Mdir(1) = mag_array*cosd(theta_array');
Mdir(2) = mag_array*sind(theta_array');
Mdir(3) = sqrt(Mdir(2)^2 + Mdir(1)^2);
dir_var = 1 - Mdir(3)/sum(mag_array);
dir_mean = acosd(Mdir(1)/Mdir(3));
if Mdir(2) < 0
    dir_mean = 360 - dir_mean;
end
if dir_mean >= 360
    dir_mean = dir_mean - 360;
elseif dir_mean < 0
    dir_mean = dir_mean + 360;
end

end

