% File name: main_plotQstar
% Description: This file is used to solve problem in page 14 of the paper.
%              To find the value of zeta, we may need to plot this graph to
%              make sure the optimization is not local

%% Basic Setup
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
    };
% We need a grid for different value of zetas
zeta_bar= 0.1;                                 % The upper bar of zeta space
M       = 10;% The number of grid we are going to use for the graph
thegrid_zeta=linspace(0.01,zeta_bar, M);% Line space to generate the total number of points we are interested
Qsim=zeros(length(thegrid_zeta),1);
CVsim  = zeros(length(thegrid_zeta),1);
WAPsim = zeros(length(thegrid_zeta),1);
R0sim  = zeros(length(thegrid_zeta),1);
R1sim  = zeros(length(thegrid_zeta),1);
exitflagsim = zeros(length(thegrid_zeta),1);
Lambdasim = zeros(length(thegrid_zeta),4);
bds_testsim = zeros(length(thegrid_zeta),4);
for i=1:length(thegrid_zeta)
    X = ['The gird point number ',num2str(i),' is running.'];
    disp(X)
    zeta=thegrid_zeta(i);
    try
        [Qstar,CV, WAP, Lambda,R0,R1,bds_test,exitflag]=outeropt_Qstar(zeta, ns, upper, lower, J1, J2, a, b,As);
        Qsim(i,1) = Qstar;
        CVsim(i)  = CV;
        WAPsim(i) = WAP;
        R0sim(i)  = R0;
        R1sim(i)  = R1;
        exitflagsim(i) = exitflag;
        Lambdasim(i,:) = Lambda;
        bds_testsim(i,:) = bds_test;
    catch
        exitflagsim(i) = 0;
    end
end
plot(thegrid_zeta(exitflagsim==1),Qsim(exitflagsim==1))
xlabel('$\zeta$','FontSize', 18, 'Interpreter','latex')
ylabel('Q^*(zeta)','FontSize', 18)