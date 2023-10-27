function Qf=obj_func(x, As, J1, J2, a, b, theta_0s, theta_1s, mu1_alt,zeta)
% This is the function for our optimization objective function
% inputs:
% x:         the optimizer
% As:        the number of events we consider in outcome space S. Notice
%            in two player entry game we only
% J1         the number of points for shape parameter
% J2         the number of points for scale parameter. Notice It decides
% a          shape parameter evaluate points
% b          scale parameter evaluate points
% theta_0s:  the null hypothesis we consider
% theta_1s:  the alternative hypothesis we consider
% mu1_alt:   the alternative prior mu1 as weightiied uniform distribution
% zeta:      the size control parameter

%%
% Step 0: calculate the eta based on input values
eta=tauvtoeta_trans(x);  % return the value of eta and the dimension of eta is 1 by J (double check in calculation)
% Step 1: choose l and compute the least favorable pairs
[q_0,q_1]=inneropt_LFpairs(eta, As, J1, J2, a, b, theta_0s, theta_1s, mu1_alt);
% Step 2: calculate the objective function Q value given x
Qf=-Q_tauv(x, q_0, q_1,zeta); %notice that we need to give a negative sign and try to use the fmincon find the minimum

