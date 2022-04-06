% noisy probe
clear all;
load('Coppola_outdoor_indoor_fit_flip_90.mat')

%% noiseless test
kappa_i = 698.50; % sensory noise
kappa_e = 19.14; % stimulus noise
kappa_c = 4.79; % category noise
p_c = 0.52; % cardinal probability
kappa_m = 38.20; % motor noise, <=700
w_c = 0.52; % categorical loss weight

res = 0.5; 
adaptor = 0;
stimulus = 0:180; % must be on 0:res:180 grid

[ bias_pred, var_pred ] = ECCategpredict_03_noisysoft( fit_spline_90, nrml_90, kappa_i, kappa_e, kappa_c, p_c, kappa_m, w_c, res, stimulus );
std_pred = sqrt(-2*log(1-var_pred))*180/pi/2;

figure(1);
set(gcf,'Position',[0, 0, 800, 300]);
subplot(1,2,1)
hold on
plot(stimulus, bias_pred, 'LineWidth', 2, 'Color', [216, 82, 24]/255);
ref_line = plot([0, 180], [0, 0], 'k', 'LineWidth', 2);
ref_line.Color(4) = 0.4;
xlim([0 180])
% ylim([-10 10])
set(gca, 'XTick',0:45:180)
xlabel('Test orientation (deg)')
ylabel('Bias (deg)')
set(gca, 'FontSize', 16)

subplot(1,2,2)
plot(stimulus, std_pred, 'LineWidth', 2, 'Color', [216, 82, 24]/255);
xlim([0 180])
ylim([0 10])
set(gca, 'XTick',0:45:180)
xlabel('Test orientation (deg)')
ylabel('SD (deg)')
set(gca, 'FontSize', 16)

%% noisy test
kappa_pi = 698.50; % probe sensory noise
kappa_pe = 19.14; % probe stimulus noise

res = 0.5; 
adaptor = 0;
stimulus = 0:180; % must be on 0:res:180 grid

[ bias_pred, var_pred ] = ECCategpredict_05_noisysoft( fit_spline_90, nrml_90, kappa_i, kappa_e, kappa_pi, kappa_pe, kappa_c, p_c, kappa_m, w_c, res, stimulus );
std_pred = sqrt(-2*log(1-var_pred))*180/pi/2;

figure(2);
set(gcf,'Position',[0, 0, 800, 300]);
subplot(1,2,1)
hold on
plot(stimulus, bias_pred, 'LineWidth', 2, 'Color', [216, 82, 24]/255);
ref_line = plot([0, 180], [0, 0], 'k', 'LineWidth', 2);
ref_line.Color(4) = 0.4;
xlim([0 180])
% ylim([-10 10])
set(gca, 'XTick',0:45:180)
xlabel('Test orientation (deg)')
ylabel('Bias (deg)')
set(gca, 'FontSize', 16)

subplot(1,2,2)
plot(stimulus, std_pred, 'LineWidth', 2, 'Color', [216, 82, 24]/255);
xlim([0 180])
ylim([0 15])
set(gca, 'XTick',0:45:180)
xlabel('Test orientation (deg)')
ylabel('SD (deg)')
set(gca, 'FontSize', 16)

