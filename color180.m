function [ out ] = color180( in )
%wrap color to [-180,180)
large = in >= 180;
small = in < -180;
while sum(large,'all') >0
    in(large) = in(large)-360;
    large = in >= 180;
end
while sum(small,'all') >0
    in(small) = in(small)+360;
    small = in < -180;
end
out = in;

end

