function [ theta ] = circ180( theta )
%wrap theta to [0, 180)
% theta: scaler, vector, or 2D matrix

large = theta >= 180;
small = theta < 0;
while (sum(sum(large))) > 0
    theta(large) = theta(large)-180;
    large = theta >= 180;
end
while (sum(sum(small))) > 0
    theta(small) = theta(small)+180;
    small = theta < 0;
end

end

