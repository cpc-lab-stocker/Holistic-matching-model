function [err_grid, sti_grid, prob_err_sti_T_smooth] = prob_smooth(err_grid0, sti_grid0, prob_err_sti_T, sigma, err_vec, sti_vec)
%smooth error distribution by convolving with a gaussian with sigma

[err_grid, sti_grid] = meshgrid(err_vec, sti_vec); %(-89.5:89.5, 0.5:1:179.5);
prob_err_sti_T_smooth = zeros(size(err_grid));
for i = 1:size(prob_err_sti_T,1)
    for j = 1:size(prob_err_sti_T,2)
        prob_err_sti_T_smooth = prob_err_sti_T_smooth + prob_err_sti_T(i,j) * normpdf(sqrt(circ90(sti_grid-sti_grid0(i,1)).^2+circ90(err_grid-err_grid0(1,j)).^2),0,sigma);
    end
end
prob_err_sti_T_smooth = prob_err_sti_T_smooth./sum(prob_err_sti_T_smooth,2); % sti * err
end

