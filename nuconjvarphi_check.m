function [y,f, nu_conj] = nuconjvarphi_check(x_m,A,a,b)
% Description: this function do numerical integration for intergral w.r.t
% upper bound density (nu_conjugate)
    theta_evaluate = [0,-x_m];
    nu_conj = get_nu_conjbds(theta_evaluate,A);
    f = gampdf(x_m,a,b);
    y = nu_conj*f;
end