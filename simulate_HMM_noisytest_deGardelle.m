% simulate model for noisy test 
% w/ best fitting parameters from deGardelle dataset
clear all;

load('Coppola_1998_analysis/Coppola_outdoor_indoor_fit_flip_90.mat')

Dur = [3,4,5,6];
load(['fit_deGardelle_HMM_combodur' num2str(Dur) '.mat'])

stimulus = 0:179;
res = 0.5;

probErrSti_pred_smooth_group = NaN(180, 180, length(Dur));
bias_pred_smooth_group = NaN(length(Dur), 180);
var_pred_smooth_group = NaN(length(Dur), 180);
std_pred_smooth_group = NaN(length(Dur), 180);

figure(1);
set(gcf,'Position',[0, 0, 1000, 750]);

figure(2);
set(gcf,'Position',[0, 0, 1200, 450]);

for kk = 1:length(Dur)
dur = Dur(kk);
kappa_i = kappa_i_group(kk);

%% simulate
[ adjust1, p1_adj_theta_T, ~, ~ ] = HMM_noisytest( fit_spline_90, nrml_90, kappa_i, kappa_e, kappa_b, p_c, kappa_m, w_c, res, stimulus );

%% estimate distribution -> error distribution
err_vec0 = -90:res:90-res;
[err_grid0, sti_grid0, probErrSti_T] = prob_est2err(stimulus, adjust1, p1_adj_theta_T, err_vec0);

% smooth model predicted distr to make comparable with data distr
sigma = 5;
err_vec = -89.5:89.5;
sti_vec = 0.5:1:179.5;
[err_grid, sti_grid, probErrSti_T_smooth] = prob_smooth(err_grid0, sti_grid0, probErrSti_T, sigma, err_vec, sti_vec);
probErrSti_pred_smooth_group(:,:,kk) = probErrSti_T_smooth';

figure(1);
subplot(2,2,length(Dur)+1-kk)
imagesc(sti_vec([1,end]), err_vec([1,end]), probErrSti_T_smooth')
colormap(flipud(gray))
xlabel('Test orientation (deg)')
ylabel('Error (deg)')
set(gca, 'XTick',0:45:180, 'YTick',-90:45:90)
set(gca,'YDir','normal')
set(gca, 'TickLength',[0;0])
set(gca, 'FontSize', 18)
colorbar

%% bias, SD
[bias_pred_smooth, var_pred_smooth, std_pred_smooth] = bias_var_std_from_prob(err_vec, probErrSti_T_smooth);

bias_pred_smooth_group(kk,:) = bias_pred_smooth;
var_pred_smooth_group(kk,:) = var_pred_smooth;
std_pred_smooth_group(kk,:) = std_pred_smooth;

figure(2);
subplot(1,2,1)
hold on
p2 = plot(sti_vec, bias_pred_smooth, 'LineWidth', 2, 'Color', [216, 82, 24]/255);
p2.Color(4) = (length(Dur)-kk+1)/length(Dur);

subplot(1,2,2)
hold on
p2 = plot(sti_vec, std_pred_smooth, 'LineWidth', 2, 'Color', [216, 82, 24]/255);
p2.Color(4) = (length(Dur)-kk+1)/length(Dur);

end

figure(2)
subplot(1,2,1)
refline = plot([0 180], [0 0], 'k', 'LineWidth', 2);
refline.Color(4) = 0.2;
xlim([0 180])
set(gca,'XTick',0:45:180)
ylim([-15 15])
set(gca,'YTick',-15:5:15)
xlabel('Test orientation (deg)')
ylabel('Bias (deg)')
set(gca, 'FontSize', 24)

subplot(1,2,2)
xlim([0 180])
set(gca,'XTick',0:45:180)
ylim([0 40])
set(gca,'YTick',0:10:40)
xlabel('Test orientation (deg)')
ylabel('SD (deg)')
set(gca, 'FontSize', 24)
