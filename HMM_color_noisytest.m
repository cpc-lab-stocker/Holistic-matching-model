function [ Bias1, Var1, adjust1, p1_adj_theta_T ] = HMM_color_noisytest( color_num, prior_num, kappa_i, kappa_e, theta_b, kappa_b, kappa_c, kappa_m, w_c, res, stimulus )

if kappa_m > 700
    error('motor noise too small.')
end

color_num = [color_num(end)-360, color_num, color_num(1)+360];
prior_num = [prior_num(end), prior_num, prior_num(1)];

x1 = 0:res:360-res; % physical space, [0,360)
prior0 = interp1(color_num, prior_num, x1);
x2 = cumtrapz(x1, prior0) * 360; % sensory space, [0,360)
dx1 = ones(1, length(x1))*res;
prior2 = 1/360;
f1 = prior0/prior2;

%% efficient encoding
p1_m1_theta = zeros(length(x1), length(x1)); % P( external sample m1 | stimulus theta ); m1 * theta
if kappa_e > 700
    p1_m1_theta = eye(length(x1));
else
    for j = 1:length(x1)
        p1_m1_theta(:,j) = circ_vmpdf(x1/180*pi,x1(j)/180*pi,kappa_e)'/180*pi;
    end
end
p1_m2_m1 = zeros(length(x1), length(x1)); % P( internal sample m2 | external sample m1 ); m2 * m1
if kappa_i > 700
    p1_m2_m1 = eye(length(x1));
else
    for j = 1:length(x1)
        % uniform von mises noise in sensory space, then transformed to physical space
        p1_m2_m1(:,j) = circ_vmpdf(x2/180*pi,x2(j)/180*pi,kappa_i)'/180*pi.*f1';
    end
end
p1_m2_theta = p1_m2_m1 * (p1_m1_theta.*dx1'); % P( internal sample m2 | stimulus theta ); m2 * theta
p1_m2_theta = p1_m2_theta./(dx1*p1_m2_theta); % normalize

%% category structure
n_cat = length(theta_b);
p_bound = circ_vmpdf(x1*pi/180, 0, kappa_b)*pi/180; 
p_bound = p_bound/sum(p_bound)/res;
prob_cat_theta = zeros(n_cat, length(x1), length(x1)); % P( category | theta ); catgory * theta(x1) * boundary jitter position 

prob0_cat_theta = zeros(n_cat, length(x1));
for j = 1:n_cat
    F2 = @(x)circ_vmpdf(x, theta_b(j)*pi/180, kappa_c); 
    if j == 1
        F1 = @(x)circ_vmpdf(x, theta_b(n_cat)*pi/180, kappa_c); 
        prob = @(x) integral(F1, theta_b(n_cat)*pi/180-3*pi(), x) - integral(F2, theta_b(j)*pi/180-pi(), x); 
    else
        F1 = @(x)circ_vmpdf(x, theta_b(j-1)*pi/180, kappa_c); 
        prob = @(x) integral(F1, theta_b(j-1)*pi/180-pi(), x) - integral(F2, theta_b(j)*pi/180-pi(), x); 
    end
    for i = 1:length(x1)
        prob0_cat_theta(j,i) = prob(x1(i)*pi/180);
    end
end
for i = 1:length(x1) % boundary jitter position
    prob_cat_theta(:,:,i) = circshift(prob0_cat_theta, i-1, 2);
end

%% matching
estimation1_I = zeros(length(x1),length(x1)); % boundary position * measurement
post2 = p1_m2_theta.*prior0;
post2 = post2./(post2*dx1'); % m2 * theta
mean_cos_x = cosd(x1) * (post2.*dx1)'; % 1 * m2
mean_sin_x = sind(x1) * (post2.*dx1)';
L_adj = 2 - 2 * ( mean_cos_x'.*cosd(x1) + mean_sin_x'.*sind(x1) ); % m2 * estimate

for i = 1:length(x1) % each boundary position
    post2_cat = prob_cat_theta(:,:,i) * (post2.*dx1)';  % n_cat * m2
    L_cat_exp = post2_cat'*(1-prob_cat_theta(:,:,i)); % m2 * estimate
    L_tot = L_cat_exp * w_c + L_adj * (1-w_c);

    [M, I] = min(L_tot, [], 2); % m2 * 1
    estimation1_I(i,:) = I';
end

%% response probability
p_adj_x1end = circ_vmpdf(x1'/180*pi,x1(end)/180*pi,kappa_m)/180*pi;
p1_adj_theta_bound = NaN(length(x1), length(x1), length(x1));
for j = 1:length(x1) % each boundary position
    p_adj_est = zeros(length(x1), length(x1));
    for i = 1:length(x1) % m2
        p_adj_est(:,i) = circshift(p_adj_x1end, estimation1_I(j,i));
    end
    p1_adj_theta_bound(:,:,j) = p_adj_est * (p1_m2_theta.*dx1');
end
p_bound_dx1 = reshape(p_bound*res, 1, 1, length(x1));
p1_adj_theta = sum(p1_adj_theta_bound.*p_bound_dx1,3);
adjust1 = x1;

%% statistics
POE1 = NaN(size(stimulus));
Var1 = NaN(size(stimulus));
p1_adj_theta_return = NaN(size(p1_adj_theta,1), length(stimulus));
for i = 1:length(stimulus)
    sti = stimulus(i);
    sti_ind = round(color360(sti)/res)+1; % x1
    [ POE1(i), Var1(i) ] = color_mean_var( adjust1, p1_adj_theta(:,sti_ind)'.*dx1 );
    p1_adj_theta_return(:,i) = p1_adj_theta(:,sti_ind);
end
Bias1 = color180(POE1-stimulus);
p1_adj_theta_T = p1_adj_theta_return';



end
