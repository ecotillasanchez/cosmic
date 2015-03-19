function status = verify_mac(delta_m, Ea, Eap, Xd, Xdp, Xq, Va, Pe, Qe)
% usage: status = verify_mac(delta_m, Ea, Eap, Xd, Xdp, Xq, Va, Pe, Qe)
status = true;

P1 = (Ea.*Va./Xd).*sin(delta_m) +...
    Va.^2/2*(1./Xq-1./Xd).*sin(2*delta_m); % eq. 6.45, B&V
Q1 = (Ea.*Va./Xd).*cos(delta_m) - ...
    Va.^2.*(cos(delta_m).^2/Xd + sin(delta_m).^2/Xq); % eq. 6.46, B&V

P2 = (Eap.*Va./Xdp).*sin(delta_m) +...
    Va.^2/2*(1./Xq-1./Xdp).*sin(2*delta_m); % eq. 6.45, B&V
Q2 = (Eap.*Va./Xdp).*cos(delta_m) - ...
    Va.^2.*(cos(delta_m).^2/Xdp + sin(delta_m).^2/Xq); % eq. 6.46, B&V

if abs(P1 - P2) > 1e-6
    disp('P part of the equations do not match');
    fprintf('P1 = %g ... P2 = %g \n', P1, P2);
    status = false;
end

if abs(P1 - Pe) > 1e-6
    disp('P part of the equations do not match');
    fprintf('P1 = %g ... Pe = %g \n', P1, P2);
    status = false;
end

if abs(Q1 - Q2) > 1e-6
    disp('Q part of the equations do not match');
    fprintf('Q1 = %g ... Q2 = %g \n', Q1, Q2);
    status = false;
end

if abs(Q1 - Qe) > 1e-6
    disp('Q part of the equations do not match');
    fprintf('Q1 = %g ... Qe = %g \n', Q1, Q2);
    status = false;
end



if ~status
    keyboard
end
