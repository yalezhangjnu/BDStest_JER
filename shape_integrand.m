function phi_j = shape_integrand(belief_A)
% Description: this function do numerical integration for intergral we
% are interested in
% inputs:
% belieff_A : the numerical integration of belief function and paramters with three dimension (# of event * dimension of shape parameter * dimension of scale paramter )
L=size(belief_A, 1);
J1=size(belief_A,2);
J2=size(belief_A,3);
J=J1*J2;
phi_j=NaN(L,J);
for i=1:L
    phi_j(i,:)=reshape(belief_A(i,:,:),1,J);
end