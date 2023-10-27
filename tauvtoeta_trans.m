function eta=tauvtoeta_trans(tauv)
% Description: this function transform the optimizer in outer optimization
%              to eta weights we need in inner optimization steps.
% inputs:
% tauv=(tau, v_1,...,v_{l-1})' dimension is J by 1
J=size(tauv,1);
eta=NaN(J,1);                  % the dimension is J by 1
for i=2:J
    eta(i-1,1)=tauv(i,1)/tauv(1,1);
end
eta(J,1)=tauv(1,1)-sum(tauv(2:end));
end