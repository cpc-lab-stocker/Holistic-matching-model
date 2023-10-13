function [ out ] = circ90( in )
%wrap orientation to [-90,90)
large = in >= 90;
small = in < -90;
while sum(large,'all') >0
    in(large) = in(large)-180;
    large = in >= 90;
end
while sum(small,'all') >0
    in(small) = in(small)+180;
    small = in < -90;
end
out = in;


end

