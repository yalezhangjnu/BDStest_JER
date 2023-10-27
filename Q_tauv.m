function [Qf, tau]=Q_tauv(x, q_0, q_1, zeta)
% description: this function calculate the value of objective function
% inputs:
% x:         the inputs vector (tau, v_1,v_2,...,v_{l-1}), the first
%            element is tau
% bds_test:  the likelihood ratio test
% zeta: size control parameter
tau=x(1);
bds_test= q_1./q_0 > tau./(zeta.*(1-tau));
R0=sum(q_0.*bds_test);              % Expectation of bds_test under q_0
R1=zeta*(1-sum(q_1.*bds_test));     % Expectation of bds_test under q_1
Qf=tau*R0+(1-tau)*zeta*R1;
end