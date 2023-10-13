function [bias_smooth, var_smooth, std_smooth] = bias_var_std_from_prob(err_vec, prob_err_sti_T_smooth)
%compute bias, variance & SD from error distribution
bias_smooth = NaN(1,size(prob_err_sti_T_smooth,1));
var_smooth = NaN(1,size(prob_err_sti_T_smooth,1));
for i = 1:size(prob_err_sti_T_smooth,1)
    [ bias_smooth(i), var_smooth(i) ] = circ_mean_var( err_vec, prob_err_sti_T_smooth(i,:) );
    bias_smooth(i) = circ90(bias_smooth(i));
end
std_smooth = sqrt(-2*log(1-var_smooth))*180/pi/2;
end

