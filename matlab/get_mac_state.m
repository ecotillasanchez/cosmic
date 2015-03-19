function [mac,exc,gov] = get_mac_state(ps,mode)
% usage: [mac,exc,gov] = get_mac_state(ps,mode)
%  ps - ps structure typical
%  mode is one of the following:
%   'classical', Ea behind Xd
%   'salient', Salient pole machine

% NOTE that this function assumes that there is one mac entry per gen, in
% the same sequence.

% grab some data
C  = psconstants;
j  = 1i;
Pg = ps.gen(:,C.ge.Pg)/ps.baseMVA;
Qg = ps.gen(:,C.ge.Qg)/ps.baseMVA;
G  = ps.bus_i(ps.gen(:,1));
theta_g     = ps.bus(G,C.bu.Vang)*pi/180;
Va          = ps.bus(G,C.bu.Vmag).*exp(1i*theta_g);
Sa          = Pg + j*Qg; 
Ia          = conj(Sa./Va);
phi         = angle(Sa);        % power factor angle for the generator

if abs(phi) > pi/2
    error('strange angle for the power factor');
end
Xd  = ps.mac(:,C.ma.Xd);
Xdp = ps.mac(:,C.ma.Xdp);
Xq  = ps.mac(:,C.ma.Xq);
mac = ps.mac;
exc = ps.exc;
gov = ps.gov;
ng = size(mac,1);
omega_0 = 2*pi*ps.frequency;

switch mode
    case 'classical'
        % find Pm in per unit
        mac(:,C.ma.Pm) = Pg; % at nominal frequency
        % find delta
        delta_m_theta = atan2( Pg , ( Qg + Va.^2 ./ Xd ) );
        delta = delta_m_theta; % changed from delta_m_theta + theta_g;
        mac(:,C.ma.delta) = delta;
        % find Ea
        Ea_mag = Pg .* Xd ./ ( Va .* sin(delta) );
        mac(:,C.ma.Ea) = Ea_mag;
        % delta_dot = 0
        mac(:,C.ma.omega) = omega_0;
    case 'salient'
        if any(Xq==0)
            error('we need Xq values to do salient pole model');
        end
        
        aprime = Va + j.*Xq.*Ia;  % + r.*Ia
        delta = angle(aprime);
        delta_m = delta - theta_g;
        Vaq = abs(Va).*cos(delta_m);
        Iad_mag = abs(Ia).*sin(delta-angle(Ia));
        Eap_mag = Vaq + Xdp.*Iad_mag; % + r.*Iad
        Ea_mag = (Xd - Xdp).*Iad_mag + Eap_mag;
        % verify the machine
        for i=1:ng
            if ~verify_mac(delta_m(i), Ea_mag(i), Eap_mag(i), ...
                    Xd(i), Xdp(i), Xq(i), abs(Va(i)), Pg(i), Qg(i))
                error('machine %d salient model is not consistent...\n',i);
            end
        end
        % save the results
        mac(:,C.ma.Ea)      = Ea_mag;
        mac(:,C.ma.Eap)     = Eap_mag;
        mac(:,C.ma.delta_m) = delta_m;
        mac(:,C.ma.delta)   = delta;
        mac(:,C.ma.omega)   = omega_0;
        mac(:,C.ma.Pm0)     = Pg; % at nominal frequency
        exc(:,C.ex.Vref)    = Ea_mag./ ps.exc(:,C.ex.Ke) + abs(Va);
        exc(:,C.ex.E1)      = Ea_mag./ ps.exc(:,C.ex.Ke);
        exc(:,C.ex.Efd)     = Ea_mag;
        gov(:,C.go.Pref)    = Pg;
        gov(:,C.go.P3)      = 0;  % ? I think it is right 
    
    otherwise
        error('not a valid mode for get_mac_state');
end
