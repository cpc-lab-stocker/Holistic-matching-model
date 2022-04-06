function [ out ] = circ180( in )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% size_in = size(in);
% for i = 1:size_in(1)
%     for j = 1:size_in(2)
%         while in(i,j) <0
%             in(i,j) = in(i,j)+180;
%         end
% 
%         while in(i,j) >= 180
%             in(i,j) = in(i,j)-180;
%         end
% %         out(i,j) = in(i,j);
%     end
% end
% out = in;

large = in >= 180;
small = in < 0;
while (sum(sum(large))) >0
    in(large) = in(large)-180;
    large = in >= 180;
end
while (sum(sum(small))) >0
    in(small) = in(small)+180;
    small = in < 0;
end
out = in;

end

