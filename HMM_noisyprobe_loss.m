function [ x1, L_cat, L_adj, L_tot, estimation0, estimation1 ] = HMM_noisyprobe_loss( Fit_prior, nrml_prior, kappa_i, kappa_e, kappa_b, p_c, kappa_m, w_c, res, measurement, cardinal )
% stimulus, measurement, cardinal: on 0:res:180-res grid points

if kappa_m > 700
    error('motor noise too small.')
end

x1 = 0:res:180-res; % physical space, [0,180)
dx1 = ones(1, length(x1))*res;
prior1 = prior_fun2num(x1, Fit_prior, nrml_prior);

x2 = space_nonuni2uni(x1, Fit_prior, nrml_prior, 180); % sensory space, [0,180)
prior2 = 1/180;

%% efficient encoding
p1_m1_theta = zeros(length(x1), length(x1)); % P( external sample m1 | stimulus theta ); m1 * theta
if kappa_e > 700
    p1_m1_theta = eye(length(x1));
else
    for j = 1:length(x1)
        p1_m1_theta(:,j) = circ_vmpdf(2*x1/180*pi,2*x1(j)/180*pi,kappa_e)'/90*pi;
    end
end
p1_m2_m1 = zeros(length(x1), length(x1)); % P( internal sample m2 | external sample m1 ); m2 * m1
if kappa_i > 700
    p1_m2_m1 = eye(length(x1));
else
    for j = 1:length(x1)
        % uniform von mises noise in sensory space, then transformed to physical space
        p1_m2_m1(:,j) = circ_vmpdf(2*x2/180*pi,2*x2(j)/180*pi,kappa_i)'/90*pi.*prior1'/prior2;
    end
end
p1_m2_theta = p1_m2_m1 * (p1_m1_theta .* dx1'); % P( internal sample m2 | stimulus theta ); m2 * theta
p1_m2_theta = p1_m2_theta./(dx1*p1_m2_theta); % normalize

%% category structure
% P( category | theta )
[cat0_CW, cat0_CCW, cat0_card0, cat0_card90] = prob_cat_theta(x1, kappa_b, p_c);

% noisy signal of cardinal position = cardinal
[M_m, I_m] = min(abs(x1-measurement));
[M_c, I_c] = min(abs(x1-cardinal));

cat_CW = circshift(cat0_CW, I_c-1);
cat_CCW = circshift(cat0_CCW, I_c-1);
cat_card0 = circshift(cat0_card0, I_c-1);
cat_card90 = circshift(cat0_card90, I_c-1);

post_csr_mest = p1_m2_theta.*prior1;
post_csr_mest = (post_csr_mest./(post_csr_mest*dx1'))'; % csr * measurement of est

%% matching
post2 = zeros(1,length(x1)); 
post2(I_m) = 1/res;
post2_CW = (post2.*cat_CW)*dx1';  % boundary position * 1
post2_CCW = (post2.*cat_CCW)*dx1';  % boundary position * 1
post2_card0 = (post2.*cat_card0)*dx1';  % boundary position * 1
post2_card90 = (post2.*cat_card90)*dx1';  % boundary position * 1

mean_cos_mest = cosd(2*x1) * (post_csr_mest.*dx1'); % 1 * est
mean_sin_mest = sind(2*x1) * (post_csr_mest.*dx1');

p_adj_CW = cat_CW * (post_csr_mest.*dx1');
p_adj_CCW = cat_CCW * (post_csr_mest.*dx1');
p_adj_card0 = cat_card0 * (post_csr_mest.*dx1');
p_adj_card90 = cat_card90 * (post_csr_mest.*dx1');

L_cat = post2_CW .* (1-p_adj_CW) + post2_CCW .* (1-p_adj_CCW) + post2_card0 .* (1-p_adj_card0) + post2_card90 .* (1-p_adj_card90); % boundary position * estimate
L_adj = 2 - 2 * ( mean_cos_mest*cosd(2*x1(I_m)) + mean_sin_mest*sind(2*x1(I_m)) );
L_tot = L_cat * w_c + L_adj * (1-w_c);

[M, I] = min(L_tot, [], 2);
estimation1 = x1(I);
[M, I] = min(L_adj, [], 2);
estimation0 = x1(I);

end
