function [Q_0, Q_1]=inneropt_LFpairs(eta, As, J1, J2, a, b, theta_0s, theta_1s, mu1_alt)
% Description: this function is compute the inner optimization Algorithm 2
% in JER paper KAIDO and Zhang (2023)
% output: 
% Q_0:       least favorable pairs under H_0
% Q_1:       least favorable pairs under H_1
% inputs:
% eta: the vector that comes from initialval or optimizer from the
%            outer optimization. optimizer=(tau, v_1,...,v_{l-1}) notice
%            the dimension of the optimizer
%            1 by l
% As:        the number of events we consider in outcome space S. Notice
%            in two player entry game we only
% J1         the number of points for shape parameter
% J2         the number of points for scale parameter. Notice It decides
% a          shape parameter evaluate points
% b          scale parameter evaluate points
% mu1_alt:   the alternative prior mu1 as weightiied uniform distribution
%% Approximate the finite parameter space using sieve method
L=length(As);  % the total number of events
Ln=length(theta_0s);
La=length(theta_1s);
%J = J1*J2;     % Total number basis function we are hiring in gamma distribution case
kappa_A=NaN(L,2);                   % construct a vector to store the lower bound of belief function under H0 and H1
kappa_conj_A=NaN(L,2);              % construct a vector to store the upper bound of belief function under H0 and H1
%% Approximate the intergral over a space
nu0f_A=NaN(L,J1,J2);        % store the approximate value based on lower bound of belief function based on different combination of J1 and J2
nu0conjf_A=NaN(L,J1,J2);    % store the approximate value based on upper bound of belief function based on different combiantion of J1 and J2
for i=1:L            % index for events fixed event i
    for j=1:J1       % index for shape parameter
        for k=1:J2   % index for scale parameter
            nu0f_A(i,j,k)     = integral(@(x)integrand_nu(x,As{i},a(j),b(k)),0,Inf);
            nu0conjf_A(i,j,k) = integral(@(x)integrand_nu_conj(x,As{i},a(j),b(k)),0,Inf);
        end
    end
end
varphinu0f_j=shape_integrand(nu0f_A);
varphiconj0f_j=shape_integrand(nu0conjf_A);
%% calculate the eta based on the new definition
for i=1:L          % index for events
    kappa_A(i,1)=varphinu0f_j(i,:)*eta;          % caluclate the integration multiply eta for as kappa_0l
    kappa_conj_A(i,1)=varphiconj0f_j(i,:)*eta;   % calculate the intergration multiply eta for as upper bound 
end
nu_alt = zeros(La,1);               % store lower bound belief function under alter theta
nu_alt_conj = zeros(La,1);          % store upper bound belief function under alter theta
for i=1:L            % index for each events
    A=As{i};         % take corresponding event
    for j=1:La       % index for each valuate points
        theta_evaluate=theta_1s(j,:);                 % take alternative hypothesis parameter value
        nu_alt(j,1)=get_nubds(theta_evaluate, A);     % Calculate the alternative parameter lower bound for event A
        nu_alt_conj(j,1)=get_nu_conjbds(theta_evaluate, A);  % Calculate the alternative parameter upper bound for event A
    end
    kappa_A(i,2)=nu_alt'*mu1_alt;                        %  put uniform prior mu1 in the paper for alternative hypothesis lower bound
    kappa_conj_A(i,2)=nu_alt_conj'*mu1_alt;              %  put uniform prior mu1 in the paper for alternative hypothesis upper bound
end
%% SETUP CONVEX PROBLEM FOR THIS CASE
%%
single_cum=0;
for j=1:length(As)
    single_cum=single_cum+(length(As{j})==1);
end
k=2*single_cum;      % This is the general case for number of solvers we
%would like to compute
%% setup convex problem
%%% Start to construct convex problem we need to solve
cvx_begin quiet
cvx_precision high
cvx_solver sedumi
variable x(k); %  there are n elements in this vector
%%% x(1) is p_0(0,0), x(2) is p_0(1,1)
%%% x(3) is p_0(0,1), x(4) is p_0(1,0)
%%% x(5) is p_1(0,0), x(6) is p_1(1,1)
%%% x(7) is p_1(0,1), x(8) is p_1(1,0)

%%% There are three parts of inequalities we need to consider

%%% Part 1: For each singleton event, there exists upper and lower
%%% bound for probability.
A1_0=[diag(ones(4,1)),zeros(4,4);
    diag(-ones(4,1)),zeros(4,4)];
b1_0=[kappa_conj_A(1:4,1);
    -kappa_A(1:4,1)];


A1_1=[zeros(4,4), diag(ones(4,1));
    zeros(4,4), diag(-ones(4,1))];
b1_1=[kappa_conj_A(1:4,2);
    -kappa_A(1:4,2)];
%%% Part 2: For nonsingleton event, there exists upper and lower
%%% bound for probability

% This matrix is for upper and lower bound of [00,11],[00,01],[00,10] under
% null hypothesis
A2_0=[1 1 0 0 0 0 0 0;
    1 0 1 0 0 0 0 0;
    1 0 0 1 0 0 0 0;
    -1 -1 0 0 0 0 0 0;
    -1 0 -1 0 0 0 0 0;
    -1 0 0 -1 0 0 0 0];
b2_0=[kappa_conj_A(5:end,1);
   -kappa_A(5:end,1)];

% This matrix is for upper and lower bound of [00,11],[00,01],[00,10] under
% alternative hypothesis

A2_1=[0 0 0 0 1 1 0 0;
    0 0 0 0 1 0 1 0;
    0 0 0 0 1 0 0 1;
    0 0 0 0 -1 -1 0 0;
    0 0 0 0 -1 0 -1 0;
    0 0 0 0 -1 0 0 -1];
b2_1=[kappa_conj_A(5:end,2);
    -kappa_A(5:end,2)];

%%% Part 3: This part has nothing to do with constrains, instead it comes
%%% from the fact that probability sum up equal to 1.

A3=[ones(1,4), zeros(1,4);
    -ones(1,4), zeros(1,4);
    zeros(1,4), ones(1,4);
    zeros(1,4), -ones(1,4)];
b3=[1;-1;1;-1];

%%% Combine 3 parts together we have constrain matrix A and bound vector b
A=[A1_0;A1_1;A2_0;A2_1;A3];
ub=[b1_0;b1_1;b2_0;b2_1;b3];

p=[x(1),x(2),x(3),x(4)];
q=[x(5),x(6),x(7),x(8)];

A*x<=ub;

minimize(sum(rel_entr(p+q,p)))

cvx_end

Q_0=p;
Q_1=q;
end