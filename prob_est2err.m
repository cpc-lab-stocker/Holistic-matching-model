function [err_grid, sti_grid, prob_err_sti_T] = prob_est2err(sti_vec, est_vec, prob_est_sti_T, err_vec)
%convert estimate distribution to error distribution. error = estimate - stimulus.
[est_grid_0, sti_grid_0] = meshgrid([est_vec, est_vec(1)+180], [sti_vec, sti_vec(1)+180]);
prob_est_sti_T = [prob_est_sti_T, prob_est_sti_T(:,1)];
prob_est_sti_T = [prob_est_sti_T; prob_est_sti_T(1,:)];
[err_grid, sti_grid] = meshgrid(err_vec, sti_vec); %(-89.5:89.5, 0.5:1:179.5);
est_grid = NaN(length(sti_vec), length(err_vec));
for i = 1:size(sti_grid,1)
    est_grid(i,:) = circ180(err_grid(i,:) + sti_grid(i,1));
end
prob_err_sti_T = interp2(est_grid_0, sti_grid_0, prob_est_sti_T, est_grid, sti_grid);
end

