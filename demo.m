% simulate noisy test, noiseless probe
clear all;
load('Coppola_outdoor_indoor_fit_flip_90.mat')

kappa_i = 15.79; % sensory noise
kappa_e = inf; % stimulus noise
kappa_c = 8.29; % category noise
p_c = 0.6; % cardinal probability
kappa_m = 34.63; % motor noise, <=700
w_c = 10.75; % categorical loss weight

res = 0.5; 
adaptor = 0;
stimulus = 0:180; % must be on 0:res:180 grid

[ estimation, p_sti_est ] = ECCategpredict_02_noisysoft( fit_spline_90, nrml_90, kappa_i, kappa_e, kappa_c, p_c, kappa_m, w_c, res, stimulus );

%% Gaussian kernel smoothing
% estimation -> error
[X,Y] = meshgrid([estimation, estimation(1)+180], stimulus);
[x_err0, y_sti0] = meshgrid(-90:res:90-res, stimulus(1:end-1)); % error & stimulus
x_est = NaN(size(x_err0));
for i = 1:size(y_sti0,1)
    x_est(i,:) = circ90(x_err0(i,:) + y_sti0(i,1)); % estimation
end
prob_sti_err = interp2(X, Y, [p_sti_est, p_sti_est(:,1)], x_est, y_sti0);

% Gaussian kernel smoothing
sigma = 5;
err_mid = -89.5:89.5;
sti_mid = 0.5:1:179.5;
[x_err, y_sti] = meshgrid(err_mid, sti_mid); %(-89.5:89.5, 0.5:1:179.5);
prob_sti_err_smooth = zeros(size(x_err));
for i = 1:size(prob_sti_err,1)
    for j = 1:size(prob_sti_err,2)
        prob_sti_err_smooth = prob_sti_err_smooth + prob_sti_err(i,j) * normpdf(sqrt(circ90(y_sti-y_sti0(i,1)).^2+circ90(x_err-x_err0(1,j)).^2),0,sigma);
    end
end
prob_sti_err_smooth = prob_sti_err_smooth./sum(prob_sti_err_smooth,2); % sti * err

figure(1);
set(gcf,'Position',[0, 0, 400, 300]);
imagesc(sti_mid([1,end]), err_mid([1,end]), prob_sti_err_smooth')
colormap(flipud(gray))
xlabel('Test orientation (deg)')
ylabel('Error (deg)')
set(gca, 'XTick',0:45:180, 'YTick',-90:45:90)
set(gca,'YDir','normal')
colorbar
set(gca, 'FontSize', 16)

%% bias, std
bias_pred_smooth = NaN(1,size(prob_sti_err_smooth,1));
var_pred_smooth = NaN(1,size(prob_sti_err_smooth,1));
for i = 1:size(prob_sti_err_smooth,1)
    [ bias_pred_smooth(i), var_pred_smooth(i) ] = circ_mean_var( err_mid, prob_sti_err_smooth(i,:) );
    bias_pred_smooth(i) = circ90(bias_pred_smooth(i));
end
std_pred_smooth = sqrt(-2*log(1-var_pred_smooth))*180/pi/2;

figure(2);
set(gcf,'Position',[0, 0, 800, 300]);
subplot(1,2,1)
hold on
plot(sti_mid, bias_pred_smooth, 'LineWidth', 2, 'Color', [216, 82, 24]/255);
ref_line = plot([0, 180], [0, 0], 'k', 'LineWidth', 2);
ref_line.Color(4) = 0.4;
xlim([0 180])
ylim([-10 10])
set(gca, 'XTick',0:45:180)
xlabel('Test orientation (deg)')
ylabel('Bias (deg)')
set(gca, 'FontSize', 16)

subplot(1,2,2)
plot(sti_mid, std_pred_smooth, 'LineWidth', 2, 'Color', [216, 82, 24]/255);
xlim([0 180])
ylim([0 15])
set(gca, 'XTick',0:45:180)
xlabel('Test orientation (deg)')
ylabel('SD (deg)')
set(gca, 'FontSize', 16)
