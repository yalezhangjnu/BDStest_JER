function y = integrand_nu(x_m,A,a,b)
% Description: this function do numerical integration for intergral we
% are interested in
% inputs:
% x_m: the evaluation points
% A:   the events
% a:   the range of first parameter
% b:   the range of second parameter
    theta_evaluate = [0,-x_m];
    nu = get_nubds(theta_evaluate,A);
    f = gampdf(x_m,a,b);
    y = nu*f;
end