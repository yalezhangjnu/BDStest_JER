function y = integrand_nu_conj(x_m,A,a,b)
% Description: this function do numerical integration for intergral w.r.t
% upper bound density (nu_conjugate)
    theta_evaluate = [0,-x_m];
    nu = get_nu_conjbds(theta_evaluate,A);
    f = gampdf(x_m,a,b);
    y = nu*f;
end