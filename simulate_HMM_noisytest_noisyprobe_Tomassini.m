% simulate model for switching test & probe
% w/ best fitting parameters from Tomassini dataset
clear all;

load('Coppola_1998_analysis/Coppola_outdoor_indoor_fit_flip_90.mat')

n_ext = 2; % number of external noise levels

load('fit_Tomassini_HMM.mat')

stimulus_pred = 0:180;
res = 0.5;

figure(1) 
set(gcf,'Position',[0, 0, 1200, 450]);
FontSize = 12;

%% test: noisy test
bias_pred_test = NaN(n_ext, length(stimulus_pred));
var_pred_test = NaN(n_ext, length(stimulus_pred));
std_pred_test = NaN(n_ext, length(stimulus_pred));
for ii = 1:n_ext
    kappa_e = kappa_e_group(ii);

    [ ~, ~, Bias1, Var1 ] = HMM_noisytest( fit_spline_90, nrml_90, kappa_i, kappa_e, kappa_b, p_c, kappa_m, w_c, res, stimulus_pred );

    bias_pred_test(ii,:) = Bias1;
    var_pred_test(ii,:) = Var1;
    std_pred_test(ii,:) = sqrt(-2*log(1-var_pred_test(ii,:)))*180/pi/2;

    subplot(2, n_ext*2, ii)
    hold on
    plot(stimulus_pred, bias_pred_test(ii,:), 'LineWidth', 2, 'Color', [216, 82, 24]/255)
    refline = plot([0 180], [0 0], 'k', 'LineWidth', 2);
    refline.Color(4) = 0.2;
    xlim([0 180])
    set(gca,'XTick',0:45:180)
    ylim([-6 6])
    set(gca,'YTick',-6:3:6)
    xlabel('Test orientation (deg)')
    ylabel('Bias (deg)')
    set(gca, 'FontSize', FontSize)

    subplot(2, n_ext*2, ii+n_ext*2)
    hold on
    plot(stimulus_pred, std_pred_test(ii,:), 'LineWidth', 2, 'Color', [216, 82, 24]/255)
    xlim([0 180])
    set(gca,'XTick',0:45:180)
    ylim([0 15])
    set(gca,'YTick',0:5:15)
    xlabel('Test orientation (deg)')
    ylabel('SD (deg)')
    set(gca, 'FontSize', FontSize)
end

%% ctrl: noisy probe
bias_pred_ctrl = NaN(n_ext, length(stimulus_pred));
var_pred_ctrl = NaN(n_ext, length(stimulus_pred));
std_pred_ctrl = NaN(n_ext, length(stimulus_pred));
for ii = 1:n_ext
    kappa_e = kappa_e_group(ii);

    [ ~, ~, Bias1, Var1 ] = HMM_noisyprobe( fit_spline_90, nrml_90, kappa_i, kappa_e, kappa_b, p_c, kappa_m, w_c, res, stimulus_pred );

    bias_pred_ctrl(ii,:) = Bias1;
    var_pred_ctrl(ii,:) = Var1;
    std_pred_ctrl(ii,:) = sqrt(-2*log(1-var_pred_ctrl(ii,:)))*180/pi/2;

    subplot(2, n_ext*2, ii+n_ext)
    hold on
    plot(stimulus_pred, bias_pred_ctrl(ii,:), 'LineWidth', 2, 'Color', [216, 82, 24]/255)
    refline = plot([0 180], [0 0], 'k', 'LineWidth', 2);
    refline.Color(4) = 0.2;
    xlim([0 180])
    set(gca,'XTick',0:45:180)
    ylim([-6 6])
    set(gca,'YTick',-6:3:6)
    xlabel('Test orientation (deg)')
    ylabel('Bias (deg)')
    set(gca, 'FontSize', FontSize)

    subplot(2, n_ext*2, ii+n_ext*3)
    hold on
    plot(stimulus_pred, std_pred_ctrl(ii,:), 'LineWidth', 2, 'Color', [216, 82, 24]/255)
    xlim([0 180])
    set(gca,'XTick',0:45:180)
    ylim([0 15])
    set(gca,'YTick',0:5:15)
    xlabel('Test orientation (deg)')
    ylabel('SD (deg)')
    set(gca, 'FontSize', FontSize)
end
