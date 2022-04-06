function [ theta ] = circ90( theta )
%wrap theta to [-90, 90)
% theta: scaler, vector, or 2D matrix

large = theta >= 90;
small = theta < -90;
while (sum(sum(large))) > 0
    theta(large) = theta(large)-180;
    large = theta >= 90;
end
while (sum(sum(small))) > 0
    theta(small) = theta(small)+180;
    small = theta < -90;
end

end

