function nu_conj= get_nu_conjbds(theta, A)
Phi = @(x)normcdf(x,0,1);
if isempty(A)
    nu_conj=0;
elseif A==0
    nu_conj=Phi(0)*Phi(0); % pr(u1<0)*pr(u2<0)
elseif A==11
    nu_conj=Phi(theta(1))*Phi(theta(2)); % pr(u1>-theta(1))*pr(u2>-theta(2))
elseif A==01
    nu_conj=Phi(-theta(1))*(1-Phi(0)); %use lower probability add up the area which multiple equl exists
elseif A==10
    nu_conj=(1-Phi(0))*Phi(-theta(2));
elseif isequal(A, [0; 11])
    nu_conj=Phi(0)*Phi(0)+Phi(theta(1))*Phi(theta(2));
elseif isequal(A, [0; 01])
    nu_conj=Phi(0)*Phi(0)+Phi(-theta(1))*(1-Phi(0));
elseif isequal(A, [0;10])
    nu_conj=Phi(0)*Phi(0)+(1-Phi(0))*Phi(-theta(2));
end
end