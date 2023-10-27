% File Name:   main_tab1fig1rep.m
% Description: This code try to replicateat the result of paper Figure1 and Table1 
%              outer optimization w.r.t. fixed value of zeta=0.065
%              To replicate figure 1 and Table 1, pleae change the value of
%              J1, J2 to J1=J2=3, J1=J2=4, and J1=J2=5 to get the full table
%              results.
% Author:      Yi & Hiro
% Date:        July,13,2023
% Refernces:   KAIDO AND ZHANG (2023):JER
%% Basic Setup
%S=100;
%Timer=NaN(S,1);
%for s=1:S   % this loop is used to compute the total time we used for it.
tic          % setup the start time of timer
zeta=0.065;   % Fix value of zeta at current stage
ns=10;       % number of grid for alternative paramter space for approximate
upper = 0;   % the upper bound for parameter space
lower = -5;  % the lower bound for approximate parameter space, originally it goes to negative infinity
J1 = 3;      % number of cases we consider for gamma distribution with first shape parameter k
J2 = 3;      % number of cases we consider for gamma distribution with second scale parameter theta
J=J1*J2;     % the total number of basis function we are using, J is the parameter l we defined in the paper
a = linspace(1,25,J1);         % shape parameter evaluate points
b = linspace(0.1,1.5,J2);      % scale parameter evaluate points
As             ={                               % we only use seven events with consideration of both upper bound and lower bound
    00;
    11;
    01;
    10;
    [00;11];
    [00;01];
    [00;10];
    };                                          % This is a 7 by 1 vector
%% Setup parameter we are interested in
[delta1,delta2] = meshgrid(linspace(lower,upper,ns+1)',linspace(lower,upper,ns+1)');   % ns+1: number of grid points on each axis
deltas          = [reshape(delta1,[],1),reshape(delta2,[],1)];                         % Reshape the all possible alternatives
theta_0s        = deltas(deltas(:,1)==0,:);     % THETA_0: delta^1 = 0, delta^2 <= 0
theta_1s        = deltas(deltas(:,1)~=0,:);     % THETA_1: delta^1 < 0, delta^2 <= 0 
delta1_null     = theta_0s(:,1);                % only setup delta_1 equals to zero, 
delta1_alt      = theta_1s(:,1);                % setup alternative hypothesis (not equal to zero same as the null hypothesis)
La              = length(delta1_alt);           % The length of grid points in alternative hypothesis
Ln              = length(delta1_null);          % The length of grid points in null hypothesis
mu1_alt         = repmat(1/La, 1, La)';         % Fixed alternative prior mu_1 as weighted uniform distribution. Take lower bound to value for approximation
%% Compute the the optimizer
% part 1: setup fimincon initial value
rng(0);
temp=  rand(J-1,1);      % The dimension minus 1 because the first elements of vector is tau
v   =  temp/(sum(temp));
tau =  sum(v);
tauv_initial=[tau;v];    % The dimension of this vector is v
% part 2: setup optimization option
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
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
[x,fval]=fmincon(objectf,tauv_initial,A, B,Aeq, Beq,lb,ub,nonlcons,options);
%% compute the optimal BDS-test
opteta = tauvtoeta_trans(x);
opttau = x(1);
[q_0opt, q_1opt]=inneropt_LFpairs(opteta, As, J1, J2, a, b, theta_0s, theta_1s, mu1_alt); % This is the LFP reported in the table 1
bds_test= q_1opt./q_0opt > opttau./(zeta.*(1-opttau));
toc
elapsedTime=toc;
%Timer(s,1)=elapsedTime;
%end
%Atimer=mean(Timer);
%% plot the LF prior in Figure 1 for different cases
M = 100;
thetagrid = linspace(-4,0,M);
varphi = zeros(J,M); 
i=1;
for j=1:J1       % index for shape parameter
    for k=1:J2   % index for scale parameter
        varphi(i,:) = gampdf(-thetagrid,a(j),b(k));
        i=i+1;
    end
end
mu0 = opteta'*varphi;
plot(thetagrid,mu0)
xlabel('$\theta^{(2)}$','FontSize', 18, 'Interpreter','latex')
ylabel('least-favorable prior','FontSize', 18)
%% Calculate all outputs in the Table 1
Lambda=q_1opt./q_0opt;          % The likelihood ratio-based test based on LFP in Table 1
CV=opttau./(zeta.*(1-opttau));  % Calculate the critical value in Table 1
Qstar=R0+R1;                    % Calculate the maximum BDS risk in Table 1
R1=(1-opttau)*zeta*(1-sum(q_1opt.*bds_test)); 
WAP=sum(q_1opt.*bds_test);      % Calculate the BDS-WAP in Table 1
R0=opttau*sum(q_0opt.*bds_test);% Calculate the size in Table 1
