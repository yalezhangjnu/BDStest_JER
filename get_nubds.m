function nu=get_nubds(theta,A)
% computes the theoretical lower bound or belief function
% theta: structure parameters, in this case the theta change to 6 by 1 and
%        we only take value theta(5) and theta(6),revised maybe needed in
%        the future for covariates case
% instead of 2 by 1
% mu:    the mean of error term
% sigma: the sd of error term
% A    : a column vector representing the event
% 
Phi = @(x)normcdf(x,0,1);
if isempty(A)
    nu=0;      % empty set with zero probability
elseif A == 0
    nu=Phi(0)*Phi(0); % pr(u1<0)*pr(u2<0)
elseif A == 11
    nu=Phi(theta(1))*Phi(theta(2)); % pr(u1>-theta(1))*pr(u2>-theta(2))
elseif A == 01
    nu=(1/2-Phi(theta(1)))*Phi(theta(2))+Phi(0)*Phi(0); % pr(0<u1<-theta1)*pr(u2>-theta2)+pr(u1<0)*pr(u2>0)
elseif A == 10
    nu=Phi(theta(1))*(1-Phi(theta(2)))+(1/2-Phi(theta(1)))*1/2;  %pr(u1>-theta(1))*pr(u2<-theta(2))+pr(0<u1<-theta(1))pr(u2<0)
elseif isequal(A, [0; 11])
    nu=Phi(0)*Phi(0)+Phi(theta(1))*Phi(theta(2)); % there is only a unique distribution for (0,0) and (1,1)
elseif isequal(A, [0;01])
    nu=Phi(0)*Phi(0)+(1/2-Phi(theta(1)))*Phi(theta(2))+Phi(0)*Phi(0);  %
elseif isequal(A, [0;10])
    nu=Phi(0)*Phi(0)+Phi(theta(1))*(1-Phi(theta(2)))+(1/2-Phi(theta(1)))*1/2;
else
    disp('error:wrong event')
end
end