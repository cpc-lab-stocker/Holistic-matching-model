function [num_prior] = prior_fun2num(x, fun_prior, nrml_prior)
%Convert prior from functional form to numerical form. num_prior(x) = sum(fun_prior(x)./nrml_prior)/n_prior
%   x: grid point
%   fun_prior: prior in functional form; cell(1,n_prior)
%   nrml_prior: normalizing factor for each functional prior
n_prior = length(nrml_prior);
if n_prior == 1
    num_prior = Fit_prior(x1)'/nrml_prior;
else
    num_prior = zeros(1, length(x));
    for i_p = 1:n_prior
        num_prior = num_prior + fun_prior{i_p}(x)'/nrml_prior(i_p)/n_prior;
    end
end
end

