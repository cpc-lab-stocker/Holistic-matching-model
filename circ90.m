function [ out ] = circ90( in )
%UNTITLED3 Summary of this function goes here
%   in must be a row vector, not a column vector
% size_in = size(in);
% for i = 1:size_in(1)
%     for j = 1:size_in(2)
%         while in(i,j) < -90
%             in(i,j) = in(i,j)+180;
%         end
% 
%         while in(i,j) >= 90
%             in(i,j) = in(i,j)-180;
%         end
%     end
% end
% out = in;

large = in >= 90;
small = in < -90;
while (sum(sum(large))) >0
    in(large) = in(large)-180;
    large = in >= 90;
end
while (sum(sum(small))) >0
    in(small) = in(small)+180;
    small = in < -90;
end
out = in;


end

