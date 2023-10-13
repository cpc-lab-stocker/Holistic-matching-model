% simulate model for undelayed & delayed color estimation.
% w/ best fitting parameters from Bae dataset
clear all;

load('color_prior.mat', 'color_angle', 'prior_CR')
load('color_category.mat', 'colors_name', 'n_colors', 'theta_name', 'kappa_name')
load('fit_color_HMM.mat', 'kappa_i', 'kappa_e', 'kappa_c', 'kappa_m', 'w_c')

%% simulate
stimulus = 0:1:360;
res = 1;

[ bias_pred_un, var_pred_un, adjust1, p_adj_theta_T_un ] = HMM_color_noisytest( color_angle, prior_CR, kappa_i(1), kappa_e, theta_name, kappa_name, kappa_c, kappa_m, w_c(1), res, stimulus );
[ bias_pred_de, var_pred_de, adjust1, p_adj_theta_T_de ] = HMM_color_noisytest( color_angle, prior_CR, kappa_i(2), kappa_e, theta_name, kappa_name, kappa_c, kappa_m, w_c(2), res, stimulus );
std_pred_un = sqrt(-2*log(1-var_pred_un))*180/pi;
std_pred_de = sqrt(-2*log(1-var_pred_de))*180/pi;

% estimate distribution -> error distribution
error_angle = -180:res:180-res;
prob_error_angle_pred_un = zeros(length(adjust1), length(stimulus));
prob_error_angle_pred_de = zeros(length(adjust1), length(stimulus));
for i = 1:length(stimulus)
    [~,n_shift] = min(abs(color180(error_angle-(adjust1(1)-stimulus(i))))); n_shift = n_shift-1;
    prob_error_angle_pred_un(:,i) = circshift(p_adj_theta_T_un(i,:)', n_shift);
    prob_error_angle_pred_de(:,i) = circshift(p_adj_theta_T_de(i,:)', n_shift);
end

%% plot
model_color = [216, 82, 24]/255;

figure(1)
set(gcf,'Position',[0, 0, 1200, 600]);

pmax = max(max(prob_error_angle_pred_un, [], 1), [], 2);
subplot(2,2,1)
hold on
imagesc(stimulus([1,end]), error_angle([1,end]), prob_error_angle_pred_un)

pmax = max([max(prob_error_angle_pred_de, [], 1), pmax], [], 2);
subplot(2,2,3)
hold on
imagesc(stimulus([1,end]), error_angle([1,end]), prob_error_angle_pred_de)

for i = [1,3]
    subplot(2,2,i)
    p0 = plot([0 360], [0 0], 'w--', 'LineWidth',1);
    p0.Color(4) = 0.4;
    box on
    colormap(flipud(gray))
    caxis([0 pmax])
    xlim([0 360])
    ylim([-90 90])
    set(gca,'XTick',0:60:360,'YTick',-90:30:90)
    xlabel('Hue angle (deg)')
    ylabel('Error (deg)')
    set(gca, 'FontSize', 18)
end

subplot(2,2,2)
hold on
for j = 1:length(theta_name)
    p0 = plot([theta_name(j), theta_name(j)], [-20, 20], 'k-', 'LineWidth',1);
    p0.Color(4) = 0.2;
end
p0 = plot([0 360], [0 0], 'k-', 'LineWidth',1);
p0.Color(4) = 0.2;
plot(stimulus, bias_pred_de, 'LineWidth', 2, 'Color', 1-(1-model_color)*1/2)
plot(stimulus, bias_pred_un, 'LineWidth', 2, 'Color', model_color)
xlim([0 360])
ylim([-20 20])
set(gca,'XTick',0:60:360,'YTick',-20:10:20)
xlabel('Hue angle (deg)')
ylabel('Bias (deg)')
set(gca, 'FontSize', 18)

subplot(2,2,4)
hold on
for j = 1:length(theta_name)
    p0 = plot([theta_name(j), theta_name(j)], [0, 30], 'k-', 'LineWidth',1);
    p0.Color(4) = 0.2;
end
plot(stimulus, std_pred_de, 'LineWidth', 2, 'Color', 1-(1-model_color)*1/2)
plot(stimulus, std_pred_un, 'LineWidth', 2, 'Color', model_color)
xlim([0 360])
ylim([0 30])
set(gca,'XTick',0:60:360,'YTick',0:5:30)
xlabel('Hue angle (deg)')
ylabel('SD (deg)')
set(gca, 'FontSize', 18)


