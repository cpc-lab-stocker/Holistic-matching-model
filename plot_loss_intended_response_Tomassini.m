% Fig.8
clear all;

load('Coppola_1998_analysis/Coppola_outdoor_indoor_fit_flip_90.mat')

kappa_e = 15;
kappa_i = inf; 
kappa_b = 12; 
p_c = 0.6;
w_c = 8/9;
kappa_m = 35;
adaptor = 0;
measurement = 90+15;
cardinal = 0;

res = 0.5;

modelColor2 = [0, 113, 188]/255;
modelColor3 = [216, 82, 24]/255;
modelColor23 = [148, 23, 81]/255;

figure(1) 
set(gcf,'Position',[0, 0, 1000, 400]);

%% noisy test
[ x1, L_cat, L_adj, L_tot, estimation0, estimation1 ] = ...
    HMM_noisytest_loss( fit_spline_90, nrml_90, kappa_i, kappa_e, kappa_b, p_c, kappa_m, w_c, res, measurement, cardinal );

figure(1)
subplot(1,2,1)
hold on
ylabel('Feature / Categorical / Total loss')
measureline = plot([measurement, measurement], [0 4], 'k-', 'LineWidth', 1);
measureline.Color(4) = 0.4;
p_adj = plot(x1, L_adj, '-', 'LineWidth', 2, 'Color', modelColor2);
adjline = plot([estimation0, estimation0], [0 4], '-', 'LineWidth', 1, 'Color', modelColor2);
totline = plot([estimation1, estimation1], [0 4], '-', 'LineWidth', 1, 'Color', modelColor3);
p_cat = plot(x1, L_cat, '-', 'LineWidth', 2, 'Color', modelColor23);
ylim([0 4])
set(gca,'YTick',0:4)
p_tot = plot(x1, L_tot, '-', 'LineWidth', 2, 'Color', modelColor3);
xlim([0 180])
set(gca,'XTick', 0:45:180)
xlabel('Response (deg)')
set(gca, 'FontSize', 16)



%% noisy probe
[ x1, L_cat, L_adj, L_tot, estimation0, estimation1 ] = ...
    HMM_noisyprobe_loss( fit_spline_90, nrml_90, kappa_i, kappa_e, kappa_b, p_c, kappa_m, w_c, res, measurement, cardinal );

figure(1)
subplot(1,2,2)
hold on
ylabel('Feature / Categorical / Total loss')
measureline = plot([measurement, measurement], [0 4], 'k-', 'LineWidth', 1);
measureline.Color(4) = 0.4;
p_adj = plot(x1, L_adj, '-', 'LineWidth', 2, 'Color', modelColor2);
adjline = plot([estimation0, estimation0], [0 4], '-', 'LineWidth', 1, 'Color', modelColor2);
totline = plot([estimation1, estimation1], [0 4], '-', 'LineWidth', 1, 'Color', modelColor3);
p_cat = plot(x1, L_cat, '-', 'LineWidth', 2, 'Color', modelColor23);
ylim([0 4])
set(gca,'YTick',0:4)
p_tot = plot(x1, L_tot, '-', 'LineWidth', 2, 'Color', modelColor3);
xlim([0 180])
set(gca,'XTick', 0:45:180)
xlabel('Response (deg)')
set(gca, 'FontSize', 16)




