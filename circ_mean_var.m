function [ theta_mean, theta_var ] = circ_mean_var( theta_array, mag_array )
% mean & std in circular space
% theta_array: orientation, 180 deg range
% mag_array: weight of each instance in theta_array

if nargin < 2
    mag_array = ones(size(theta_array));
end

Mdir = zeros(1,3);
Mdir(1) = sum(mag_array.*cosd(2*theta_array));
Mdir(2) = sum(mag_array.*sind(2*theta_array));
Mdir(3) = sqrt(Mdir(2)^2 + Mdir(1)^2);
theta_var = 1 - Mdir(3)/sum(mag_array);
theta_mean = acos(Mdir(1)/Mdir(3))/pi*180;
if Mdir(2) < 0
    theta_mean = 360 - theta_mean;
end
theta_mean = theta_mean/2;
if theta_mean >= 180
    theta_mean = theta_mean - 180;
elseif theta_mean < 0
    theta_mean = theta_mean + 180;
end

end

