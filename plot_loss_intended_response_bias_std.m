% Fig.2

clear all;

modelColor2 = [0, 113, 188]/255;
modelColor3 = [216, 82, 24]/255;
modelColor23 = [148, 23, 81]/255;

figure(1)
set(gcf,'Position',[0, 0, 1000, 500]);

%% prior
load('Coppola_1998_analysis/Coppola_outdoor_indoor_fit_flip_90.mat')
n_prior = length(nrml_90);
x = 0:0.2:180;
prior_pred = zeros(1, length(x));
for i_p = 1:n_prior
    prior_pred = prior_pred + fit_spline_90{i_p}(x)'/nrml_90(i_p)/n_prior;
end

figure(1)
subplot(2,3,1)
hold on
plot(x,prior_pred, 'k-', 'LineWidth', 2)
xlim([0 180])
set(gca,'XTick',0:45:180)
ylim([0, 0.022])
set(gca,'YTick',0:0.01:0.02)
xlabel('Orientation (deg)')
ylabel('Probability')
set(gca, 'FontSize', 12)

%% likelihood, posterior, loss
res = 0.5;

measurement = 90+10;
cardinal = 0;

kappa_i = 15; 
kappa_e = inf; 
kappa_b = 8; 
p_c = 0.4; 
w_c = 8/9; 
kappa_m = 35;

[ x1, L_cat, L_adj, L_tot, estimation0, estimation1, likelihood_test, posterior_test ] = ...
    HMM_noisytest_loss( fit_spline_90, nrml_90, kappa_i, kappa_e, kappa_b, p_c, kappa_m, w_c, res, measurement, cardinal );

subplot(2,3,2)
hold on
likearea = area(x1', likelihood_test');
    likearea.FaceColor = [0,0,0];
    likearea.FaceAlpha = 0.2;
    likearea.EdgeColor = 'none';
    likearea.ShowBaseLine = 'off';
postcurve = plot(x1, posterior_test, '-', 'LineWidth', 2, 'Color', modelColor2);
measureline = plot([measurement, measurement], [0 1.1*max(posterior_test)], 'k-', 'LineWidth', 1);
measureline.Color(4) = 0.4;
adjline = plot([estimation0, estimation0], [0 1.1*max(posterior_test)], '-', 'LineWidth', 1, 'Color', modelColor2);
totline = plot([estimation1, estimation1], [0 1.1*max(posterior_test)], '-', 'LineWidth', 1, 'Color', modelColor3);
ylim([0 1.1*max(posterior_test)])
set(gca,'YTick',[])
xlim([0 180])
set(gca,'XTick', 0:45:180)
xlabel('Orientation (deg)')
ylabel('Probability')
set(gca, 'FontSize', 12)

subplot(2,3,5)
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
set(gca, 'FontSize', 12)


%% category structure
stimulus = 0:179;

[p_CW_theta, p_CCW_theta, p_card0_theta, p_card90_theta] = prob_cat_theta(stimulus, kappa_b, p_c);

figure(1)
subplot(2,3,4)
hold on
plot(stimulus, p_card0_theta, 'k-', 'LineWidth', 2)
plot(stimulus, p_card90_theta, 'k-', 'LineWidth', 2)
plot(stimulus, p_CW_theta, 'k-', 'LineWidth', 2)
plot(stimulus, p_CCW_theta, 'k-', 'LineWidth', 2)
xlim([0 180])
set(gca,'XTick',0:90:180)
xlabel('Orientation (deg)')
ylim([0 1.1])
set(gca,'YTick',0:0.5:1)
ylabel('p (C | orientation)')
set(gca, 'FontSize', 12)


%% biad & SD
w_cs = [0, 1/2, 2/3, 4/5, 8/9];
model_Color_2to3 = [...
    0, 113, 188;...
    89, 79, 183;...
    179, 44, 178;...
    198, 64, 101;...
    216, 82, 24]/255;

figure(1)
subplot(2,3,3)
hold on
refline = plot([0 180], [0 0], 'k', 'LineWidth', 1);
refline.Color(4) = 0.2;
subplot(2,3,6)
hold on

for i = 1:length(w_cs)
    w_c = w_cs(i);
    model_Color_i = model_Color_2to3(i,:);
    
    [ ~, ~, bias_pred, var_pred ] = HMM_noisytest( fit_spline_90, nrml_90, kappa_i, kappa_e, kappa_b, p_c, kappa_m, w_c, res, stimulus );
    std_pred = sqrt(-2*log(1-var_pred))*180/pi/2;

    subplot(2,3,3)
    plot(stimulus, bias_pred, 'LineWidth', 2, 'Color', model_Color_i);
    
    subplot(2,3,6)
    plot(stimulus, std_pred, 'LineWidth', 2, 'Color', model_Color_i);
end

figure(1)
subplot(2,3,3)
xlim([0 180])
set(gca,'XTick',0:45:180)
ylim([-9.6 9.6]*1.2)
set(gca,'YTick',-8:4:8)
xlabel('Test orientation (deg)')
ylabel('Bias (deg)')
set(gca, 'FontSize', 12)

figure(1)
subplot(2,3,6)
xlim([0 180])
set(gca,'XTick',0:45:180)
ylim([0 18])
set(gca,'YTick',0:4:16)
xlabel('Test orientation (deg)')
ylabel('SD (deg)')
set(gca, 'FontSize', 12)

