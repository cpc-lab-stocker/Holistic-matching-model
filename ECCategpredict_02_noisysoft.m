function [ adjust, p_theta_adj_return ] = ECCategpredict_02_noisysoft( Fit_prior, nrml_prior, kappa_i, kappa_e, kappa_b, p_c, kappa_m, w_c, res, stimulus )
%% noisy test, noiseless probe

if kappa_m > 700
    error('motor noise too small.')
end

n_prior = length(nrml_prior);

x1 = circ180(-90:res:90-res); % physical space, [0,180)
x2 = NaN(1, length(x1));
for i = 1:length(x2)
    if n_prior ==1
    x2(i) = circ90(integrate(Fit_prior, x1(i), 0)/nrml_prior*180); % internal space, [-90,90)
    else
        x2(i) = 0;
        for i_p = 1:n_prior
            x2(i) = x2(i) + integrate(Fit_prior{i_p}, x1(i), 0)/nrml_prior(i_p)*180/n_prior;
        end
        x2(i) = circ90(x2(i)); % internal space, [-90,90)
    end
end
if n_prior ==1
    prior0 = Fit_prior(x1)'/nrml_prior;
else
    prior0 = zeros(1, length(x1));
    for i_p = 1:n_prior
        prior0 = prior0 + Fit_prior{i_p}(x1)'/nrml_prior(i_p)/n_prior;
    end
end
dx1 = ones(1, length(x1))*res;
prior2 = 1/180;
f1 = prior0/prior2;

%% likelihood
p_m2_m1 = zeros(length(x2), length(x2)); % same m1 in same column
p_m1_theta = zeros(length(x1), length(x1)); % same theta in same column
if kappa_e > 700
    p_m1_theta = eye(length(x2));
else
    for j = 1:length(x2)
        p_m1_theta(:,j) = circ_vmpdf(2*x1/180*pi,2*x1(j)/180*pi,kappa_e)'/90*pi .* dx1';
    end
end
if kappa_i > 700
    p_m2_m1 = eye(length(x2));
else
    for j = 1:length(x2)
        p_m2_m1(:,j) = circ_vmpdf(2*x2/180*pi,2*x2(j)/180*pi,kappa_i)'/90*pi.*f1';
    end
end
p_m2_theta = p_m2_m1 * p_m1_theta; % same theta in same column, m2 in x1
p_m2_theta = p_m2_theta./(dx1*p_m2_theta);

%% category
[M0,I0] = min(x1);
[M179,I179] = max(x1);
cond_chop0_card0 = zeros(1,length(x1));
if kappa_b > 700
    for i = 1:length(x1)
        if x1(i) == 0
            cond_chop0_card0(i) = p_c;
        end
    end
else
    for i = 1:length(x1)
        cond_chop0_card0(i) = circ_vmpdf( 2*x1(i)/180*pi, 0/180*pi, kappa_b )/circ_vmpdf( 0, 0/180*pi, kappa_b ) * p_c; 
    end
end
cond_chop0_card90 = zeros(1,length(x1));
if kappa_b > 700
    for i = 1:length(x1)
        if x1(i) == 90
            cond_chop0_card90(i) = p_c;
        end
    end
else
    for i = 1:length(x1)
        cond_chop0_card90(i) = circ_vmpdf( 2*x1(i)/180*pi, 2*90/180*pi, kappa_b )/circ_vmpdf( 2*90/180*pi, 2*90/180*pi, kappa_b ) * p_c; 
    end
end
cond_chop0_CW = zeros(1,length(x1));
if kappa_b > 700
    for i = 1:length(x1)
        if x1(i) == 0 || x1(i) == 90
            cond_chop0_CW(i) = 0.5;
        elseif x1(i)<90 && x1(i) > 0
            cond_chop0_CW(i) = 1;
        end
    end
else
    F = @(x)circ_vmpdf( x, 0, kappa_b )-circ_vmpdf( x, pi, kappa_b );
    for i = 1:length(x1)
        cond_chop0_CW(i) = integral(F,0,2*x1(i)/180*pi)+0.5-integral(F,0,0); 
    end
end
cond_chop0_CW = cond_chop0_CW .* (1-cond_chop0_card0-cond_chop0_card90);
cond_chop0_CCW = zeros(1,length(x1));
if kappa_b > 700
    for i = 1:length(x1)
        if x1(i) == 90 || x1(i) == 180
            cond_chop0_CCW(i) = 0.5;
        elseif x1(i)<180 && x1(i) > 90
            cond_chop0_CCW(i) = 1;
        end
    end
else
    F = @(x)circ_vmpdf( x, pi, kappa_b )-circ_vmpdf( x, 2*pi, kappa_b );
    for i = 1:length(x1)
        cond_chop0_CCW(i) = integral(F,0,2*x1(i)/180*pi)+0.5-integral(F,0,pi); 
    end
end
cond_chop0_CCW = cond_chop0_CCW .* (1-cond_chop0_card0-cond_chop0_card90);

cond_chop0_tot = cond_chop0_CW + cond_chop0_CCW + cond_chop0_card0 + cond_chop0_card90;
cond_chop0_CW = cond_chop0_CW./cond_chop0_tot;
cond_chop0_CCW = cond_chop0_CCW./cond_chop0_tot;
cond_chop0_card0 = cond_chop0_card0./cond_chop0_tot;
cond_chop0_card90 = cond_chop0_card90./cond_chop0_tot;

cond_chop_CW = zeros(length(x1),length(x1));
cond_chop_CCW = zeros(length(x1),length(x1));
cond_chop_card0 = zeros(length(x1),length(x1));
cond_chop_card90 = zeros(length(x1),length(x1));
for i = 1:length(x1)
    cond_chop_CW(i,:) = circshift(cond_chop0_CW, i-I0);
    cond_chop_CCW(i,:) = circshift(cond_chop0_CCW, i-I0);
    cond_chop_card0(i,:) = circshift(cond_chop0_card0, i-I0);
    cond_chop_card90(i,:) = circshift(cond_chop0_card90, i-I0);
end

%% matching
estimation = zeros(length(x1),length(x1)); % boundary position * measurement
for j = 1:length(x2)
    post = p_m2_theta(j,:).*f1; %%%%%%%%
    post = post./(post*dx1');
    post_CW = (post.*cond_chop_CW)*dx1';  % boundary position * 1
    post_CCW = (post.*cond_chop_CCW)*dx1';  % boundary position * 1
    post_card0 = (post.*cond_chop_card0)*dx1';  % boundary position * 1
    post_card90 = (post.*cond_chop_card90)*dx1';  % boundary position * 1
    
    mean_cos_x = cosd(2*x1) * (post.*dx1)';
    mean_sin_x = sind(2*x1) * (post.*dx1)';
    
    L_cat_exp = (post_CW .* (1-cond_chop_CW) + post_CCW .* (1-cond_chop_CCW) + post_card0 .* (1-cond_chop_card0) + post_card90 .* (1-cond_chop_card90) ) * w_c; % boundary position * estimate
    L_adj = 2 - 2 * ( mean_cos_x*cosd(2*x1) + mean_sin_x*sind(2*x1) );
    L_tot = L_cat_exp + L_adj;
    
    [M, I] = min(L_tot, [], 2);
    estimation(:,j) = x1(I);
        
end

%% response distribution
p_est_theta_bound_dest = repmat(p_m2_theta.*dx1', 1, 1, length(x1));

% integrate over motor noise: est -> adj
p_adj_theta_bound = NaN(size(p_est_theta_bound_dest));
for j = 1:length(x1) % bound
    p_adj_est = circ_vmpdf(2*x1'/180*pi,2*estimation(j,:)/180*pi,kappa_m)/90*pi;
    p_adj_theta_bound(:,:,j) = p_adj_est * p_est_theta_bound_dest(:,:,j);
end

% integrate over categorical noise
if kappa_b > 700
    p_bound_dx1 = zeros(1,1,length(x1));
    p_bound_dx1(I0) = 1;
else
    p_bound_dx1 = circ_vmpdf(2*x1/180*pi,2*x1(I0)/180*pi,kappa_b)/90*pi .* dx1;
    p_bound_dx1 = reshape(p_bound_dx1, 1, 1, length(x1));
end
p_adj_theta = sum(p_adj_theta_bound.*p_bound_dx1,3);
adjust = x1;
    
%% prob density of estimation1 given stimulus(theta)
p_adj_theta_return = NaN(size(p_adj_theta,1), length(stimulus));
for i = 1:length(stimulus)
    sti = stimulus(i);
    sti_ind = round(circ180(sti+90)/res)+1; % x1(sti_ind) = sti
    p_adj_theta_return(:,i) = p_adj_theta(:,sti_ind);
end
adjust = circ90(adjust); 
p_theta_adj_return = p_adj_theta_return';



end
