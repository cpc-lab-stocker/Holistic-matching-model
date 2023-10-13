function [ out ] = circ180( in )
%wrap orientation to [0,180)
large = in >= 180;
small = in < 0;
while sum(large,'all') >0
    in(large) = in(large)-180;
    large = in >= 180;
end
while sum(small,'all') >0
    in(small) = in(small)+180;
    small = in < 0;
end
out = in;

end

