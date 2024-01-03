function [ adjust1, p1_adj_theta_T, Bias, Var ] = HMM_noisyprobe( Fit_prior, nrml_prior, kappa_i, kappa_e, kappa_b, p_c, kappa_m, w_c, res, stimulus )
% stimulus: on 0:res:180-res grid points

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

% noisy boudary position
cat_CW = zeros(length(x1),length(x1)); % boundary position * theta(x1)
cat_CCW = zeros(length(x1),length(x1));
cat_card0 = zeros(length(x1),length(x1));
cat_card90 = zeros(length(x1),length(x1));
for i = 1:length(x1)
    cat_CW(i,:) = circshift(cat0_CW, i-1);
    cat_CCW(i,:) = circshift(cat0_CCW, i-1);
    cat_card0(i,:) = circshift(cat0_card0, i-1);
    cat_card90(i,:) = circshift(cat0_card90, i-1);
end


post_csr_mest = p1_m2_theta.*prior1;
post_csr_mest = (post_csr_mest./(post_csr_mest*dx1'))'; % csr * measurement of est

%% matching
estimation1 = zeros(length(x1),length(x1)); % boundary position * stimulus
estimation1_ind = zeros(length(x1),length(x1)); % boundary position * stimulus
for j = 1:length(x1) % each measurement m2
    post2 = zeros(1,length(x1)); 
    post2(j) = 1/res;
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
    
    L_cat_exp = post2_CW .* (1-p_adj_CW) + post2_CCW .* (1-p_adj_CCW) + post2_card0 .* (1-p_adj_card0) + post2_card90 .* (1-p_adj_card90); % boundary position * estimate
    L_adj = 2 - 2 * ( mean_cos_mest*cosd(2*x1(j)) + mean_sin_mest*sind(2*x1(j)) );
    L_tot = L_cat_exp * w_c + L_adj * (1-w_c);
    
    [M, I] = min(L_tot, [], 2);
    estimation1(:,j) = x1(I);
    estimation1_ind(:,j) = I;
        
end

%% response probability
if kappa_b > 700
    p_bound_dx1 = zeros(1,1,length(x1));
    p_bound_dx1(1) = 1;
else
    p_bound_dx1 = circ_vmpdf(2*x1/180*pi,2*x1(1)/180*pi,kappa_b)/90*pi .* dx1;
    p_bound_dx1 = reshape(p_bound_dx1, 1, 1, length(x1));
end

% motor noise
if kappa_m>700
    p_adj = [zeros(length(x1)-1, 1); 1/res];
else
    p_adj = circ_vmpdf(2*x1/180*pi,2*x1(end)/180*pi,kappa_m)'/90*pi;
end
p_adj_circshift = NaN(length(x1), length(x1)); % p( adjust | theta_probe); adj * est
for i = 1:length(x1)
    p_adj_circshift(:, i) = circshift(p_adj, i);
end

% posterior of probe p( theta | m2 ) assuming uniform prior
p_theta_m2_postest = (p1_m2_theta./(p1_m2_theta*dx1'))';

% integrate over motor noise
p_adj_mest_postest = p_adj_circshift * (p_theta_m2_postest .* dx1'); % p( adjust | m2_probe); adj * measurement of probe
p_adj_mest_postest = p_adj_mest_postest./(dx1*p_adj_mest_postest); % normalize

% p( adj | theta, boundary)
p1_adj_theta_bound = NaN(length(x1), length(x1), length(x1));
for j = 1:length(x1) % each boundary position
    p1_adj_theta_bound(:,:,j) = p_adj_mest_postest(:,estimation1_ind(j,:));
end
% integrate over boundary positions
p1_adj_theta = sum(p1_adj_theta_bound.*p_bound_dx1,3);
adjust1 = x1;
    
%% statistics
POE = NaN(size(stimulus)); % mean response
Var = NaN(size(stimulus)); % circular variance of response
p1_adj_theta_return = NaN(size(p1_adj_theta,1), length(stimulus)); % response distribution at stimulus
for i = 1:length(stimulus)
    sti = stimulus(i);
    sti_ind = round(circ180(sti)/res)+1; % index of stimulus in x1
    p1_adj_theta_return(:,i) = p1_adj_theta(:,sti_ind);
    [ POE(i), Var(i) ] = circ_mean_var( x1, p1_adj_theta(:,sti_ind)'.*dx1 );
end
Bias = circ90(POE-stimulus);
p1_adj_theta_T = p1_adj_theta_return';
end
