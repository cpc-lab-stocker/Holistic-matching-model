function [cat0_CW, cat0_CCW, cat0_card0, cat0_card90] = prob_cat_theta(x, kappa_b, p_c)
%P(C|theta) from category parameters
cat0_card0 = zeros(1,length(x));
if kappa_b > 700
    for i = 1:length(x)
        if x(i) == 0
            cat0_card0(i) = p_c;
        end
    end
else
    for i = 1:length(x)
        cat0_card0(i) = circ_vmpdf( 2*x(i)/180*pi, 0/180*pi, kappa_b )/circ_vmpdf( 0, 0/180*pi, kappa_b ) * p_c; 
    end
end
cat0_card90 = zeros(1,length(x));
if kappa_b > 700
    for i = 1:length(x)
        if x(i) == 90
            cat0_card90(i) = p_c;
        end
    end
else
    for i = 1:length(x)
        cat0_card90(i) = circ_vmpdf( 2*x(i)/180*pi, 2*90/180*pi, kappa_b )/circ_vmpdf( 2*90/180*pi, 2*90/180*pi, kappa_b ) * p_c; 
    end
end
cat0_CW = zeros(1,length(x));
if kappa_b > 700
    for i = 1:length(x)
        if x(i) == 0 || x(i) == 90
            cat0_CW(i) = 0.5;
        elseif x(i)<90 && x(i) > 0
            cat0_CW(i) = 1;
        end
    end
else
    F = @(x)circ_vmpdf( x, 0, kappa_b )-circ_vmpdf( x, pi, kappa_b );
    for i = 1:length(x)
        cat0_CW(i) = integral(F,0,2*x(i)/180*pi)+0.5-integral(F,0,0); 
    end
end
cat0_CW = cat0_CW .* (1-cat0_card0-cat0_card90);
cat0_CCW = zeros(1,length(x));
if kappa_b > 700
    for i = 1:length(x)
        if x(i) == 90 || x(i) == 180
            cat0_CCW(i) = 0.5;
        elseif x(i)<180 && x(i) > 90
            cat0_CCW(i) = 1;
        end
    end
else
    F = @(x)circ_vmpdf( x, pi, kappa_b )-circ_vmpdf( x, 2*pi, kappa_b );
    for i = 1:length(x)
        cat0_CCW(i) = integral(F,0,2*x(i)/180*pi)+0.5-integral(F,0,pi); 
    end
end
cat0_CCW = cat0_CCW .* (1-cat0_card0-cat0_card90);

cat0_tot = cat0_CW + cat0_CCW + cat0_card0 + cat0_card90;
cat0_CW = cat0_CW./cat0_tot;
cat0_CCW = cat0_CCW./cat0_tot;
cat0_card0 = cat0_card0./cat0_tot;
cat0_card90 = cat0_card90./cat0_tot;
end

