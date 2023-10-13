function [x2] = space_nonuni2uni(x1, fun_prior, nrml_prior, range)
%change to a space where the prior is uniform
%   x1: grid point in non-uniform space
%   fun_prior: prior in functional form; cell(1,n_prior)
%   nrml_prior: normalizing factor for each functional prior
%   range: range of uniform space
%   x2: grid point in uniform space
x2 = zeros(1, length(x1));
n_prior = length(nrml_prior);
if n_prior == 1
    for i = 1:length(x1)
        x2(i) = integrate(fun_prior, x1(i), 0)/nrml_prior*range;
    end
else
    for i = 1:length(x1)
        for i_p = 1:n_prior
            x2(i) = x2(i) + integrate(fun_prior{i_p}, x1(i), 0)/nrml_prior(i_p)*range/n_prior;
        end
    end
end
end

