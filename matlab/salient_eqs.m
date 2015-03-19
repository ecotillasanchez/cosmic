function [err,J] = salient_eqs(Ea_delta,Va,Xd,Xq,Pg,Qg)
% these are the two algebraic equations for the simplified
% salient pole generator model from Bergen and Vittal

Ea = Ea_delta(1);
delta = Ea_delta(2);

err = zeros(2,1);
err(1) = (Ea.*Va./Xd).*sin(delta) +...
         Va.^2/2*(1./Xq-1./Xd).*sin(2*delta) - Pg; % eq. 6.45, B&V
err(2) = (Ea.*Va./Xd).*cos(delta) - ...
         Va.^2.*(cos(delta).^2/Xd + sin(delta).^2/Xq) - Qg; % eq. 6.46, B&V

J = zeros(2,2);
%dP_dEa
J(1,1) = Va./Xd.*sin(delta);
%dP_ddelta
J(1,2) = (Ea.*Va./Xd).*cos(delta) + Va.^2*(1./Xq-1./Xd).*cos(2*delta);
%dQ_dEa
J(2,1) = Va./Xd.*cos(delta);
%dQ_ddelta
J(2,2) = -Ea.*Va./Xd.*sin(delta) -...
         Va.^2.*(-sin(2*delta)/Xd + sin(2*delta)/Xq);
