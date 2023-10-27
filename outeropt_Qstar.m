% This is a function file to calculation Qstar
function [Qstar,CV, WAP, Lambda,R0,R1,bds_test,exitflag]=outeropt_Qstar(zeta, ns, upper, lower, J1, J2, a, b,As)
J=J1*J2;
%% Setup parameter we are interested in
[delta1,delta2] = meshgrid(linspace(lower,upper,ns+1)',linspace(lower,upper,ns+1)');   % ns+1: number of grid points on each axis
deltas          = [reshape(delta1,[],1),reshape(delta2,[],1)];                         % Reshape the all possible alternatives
theta_0s        = deltas(deltas(:,1)==0,:);     % THETA_0: delta^1 = 0, delta^2 <= 0
theta_1s        = deltas(deltas(:,1)~=0,:);     % THETA_1: delta^1 < 0, delta^2 <= 0 
delta1_null     = theta_0s(:,1);                % only setup delta_1 equals to zero, 
delta1_alt      = theta_1s(:,1);                % setup alternative hypothesis (not equal to zero same as the null hypothesis)
La              = length(delta1_alt);           % The length of grid points in alternative hypothesis
Ln              = length(delta1_null);          % The length of grid points in null hypothesis
mu1_alt         = repmat(1/La, 1, La)';         % Fixed alternative prior mu_1 as weighted uniform distribution. 
%% Compute the the optimizer
% part 1: setup fimincon initial value
rng(0);
temp=  rand(J-1,1);      % The dimension minus 1 because the first elements of vector is tau
v   =  temp/(sum(temp));
tau =  sum(v);
tauv_initial=[tau;v];    % The dimension of this vector is v
% part 2: setup optimization option
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','StepTolerance',1e-10);
% part 3 : setup objective function 
objectf=@(x)obj_func(x, As, J1, J2, a, b, theta_0s, theta_1s, mu1_alt,zeta);
% part 4 : setup upper and lower bound
lb=zeros(J,1);
ub=ones(J,1);
% part 5: setup inequalitis and equalities
A=[-1 ones(1,J-1)];                  
B=0;
Aeq=[];
Beq=[];
nonlcons=[];
% part 5: solve optimization problem
[x,~,exitflag]=fmincon(objectf,tauv_initial,A, B,Aeq, Beq,lb,ub,nonlcons,options);
%% compute the optimal BDS-test
opteta = tauvtoeta_trans(x);
opttau = x(1);
[q_0opt, q_1opt]=inneropt_LFpairs(opteta, As, J1, J2, a, b, theta_0s, theta_1s, mu1_alt);
bds_test= q_1opt./q_0opt > opttau./(zeta.*(1-opttau));
Lambda=q_1opt./q_0opt;
CV=opttau./(zeta.*(1-opttau));  % Calculate the critical value
R1=(1-opttau)*zeta*(1-sum(q_1opt.*bds_test)); 
WAP=sum(q_1opt.*bds_test);% BDS weighted power
R0=opttau*sum(q_0opt.*bds_test);
Qstar=R0+R1;  % Calculate the maximum BDS risk
end