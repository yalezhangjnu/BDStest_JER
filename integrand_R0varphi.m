function y = integrand_R0varphi(x_m,bds_test,outcomes, a,b)
% Description: this function do numerical integration for intergral we are interested \
% inputs:
% x_m:         the evaluation points
% bds_test:    the bdstesat from algorithm inner optimization
% outocome:    the single outcome of events 
% a:           the range of first parameter
% b:           the range of second parameter
theta_evaluate = [0,-x_m];
L = length(outcomes);
nu_conj = NaN(L,1);
for i=1:L
    A=outcomes{i};
    nu_conj(i,1) = get_nu_conjbds(theta_evaluate,A);
end
R0=sum(nu_conj.*bds_test');  % here for all outcomes instead of a single outcome

f = gampdf(x_m,a,b);
y = R0*f;
end