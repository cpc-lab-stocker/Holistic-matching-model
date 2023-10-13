function [ out ] = color360( in )
%wrap color to [0,360)
large = in >= 360;
small = in < 0;
while sum(large,'all') >0
    in(large) = in(large)-360;
    large = in >= 360;
end
while sum(small,'all') >0
    in(small) = in(small)+360;
    small = in < 0;
end
out = in;


end

