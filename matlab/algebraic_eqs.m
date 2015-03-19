function [g,dg_dx,dg_dy] = algebraic_eqs(t,x,y,ps,opt) %#ok<INUSL>
% usage: [g,dg_dx,dg_dy] = algebraic_eqs(t,x,y,ps,opt)
% algebraic equation that describes the network
%
% inputs:
%   t   -> time in seconds
%   x   -> [delta delta_dot Pm Eap] for each machine
%   y   -> [Vmag Theta]
%   ps  -> power system structure
%   opt -> options structure
%
% outputs: 
%   g(1) -> Pg, real power
%   g(2) -> Qg, reactive power
%
% Equations and notation are based on Bergen & Vittal, 2000

% constants
C           = psconstants;
n           = size(ps.bus,1);
ng          = size(ps.mac,1);
m           = size(ps.branch,1);
n_sh        = size(ps.shunt,1);
ix          = get_indices(n,ng,m,n_sh,opt);

angle_ref = opt.sim.angle_ref;                 % angle reference: 0:delta_sys,1:delta_coi
COI_weight = opt.sim.COI_weight;               % weight of center of inertia

if COI_weight
    weight      = ps.mac(:,C.ma.M);
else
    weight     = ps.gen(:,C.ge.mBase);
end

% extract differential variables
deltas      = x(ix.x.delta);
Eaps        = x(ix.x.Eap);

% extract algebraic variables and fix the slack bus angle
mac_buses   = ps.mac(:,C.mac.gen);
mac_bus_i   = ps.bus_i(mac_buses);
is_slack    = ps.bus(:,C.bu.type)==C.REF;

Vmags       = y(ix.y.Vmag);
Thetas      = y(ix.y.theta);
Theta_slack = ps.bus(is_slack,C.bu.Vang)*pi/180;

mac_Vmags   = Vmags(mac_bus_i);
mac_Thetas  = Thetas(mac_bus_i);

% machine angles, relative to the bus angles:
if ~angle_ref
    delta_sys = y(ix.y.delta_sys);
    delta_m = deltas + delta_sys - mac_Thetas;
else
    delta_coi = sum(weight.*deltas)/sum(weight);
    delta_m = deltas - delta_coi - mac_Thetas;
end

% initialize output
g = zeros(length(y),1);

% extract parameters from ps
Xdps        = ps.mac(:,C.ma.Xdp);
Xqs         = ps.mac(:,C.ma.Xq);

% % calculate power from bus i to bus j 
% V       = Vmags.*exp(1i*Thetas);
% S_trans = V.*conj(Ybus*V);

% sin/cos version of S_trans
[I,K,y_IK] = find(ps.Ybus);
Vmag_I = Vmags(I);
Vmag_K = Vmags(K);
Vmag_IK = Vmag_I .* Vmag_K;
theta_IK = Thetas(I) - Thetas(K);
sin_theta_IK = sin(theta_IK);
cos_theta_IK = cos(theta_IK);
g_IK = real(y_IK);
b_IK = imag(y_IK);
P_IK = Vmag_IK .* (g_IK.*cos_theta_IK + b_IK.*sin_theta_IK);
Q_IK = Vmag_IK .* (g_IK.*sin_theta_IK - b_IK.*cos_theta_IK);
P_trans = sparse(I,1,P_IK,n,1);
Q_trans = sparse(I,1,Q_IK,n,1);

% calculate load
load_locs = ps.bus_i(ps.shunt(:,1));
Vd = Vmags(load_locs);
Pd_base = ps.shunt(:,C.sh.P) .* ps.shunt(:,C.sh.factor)/ps.baseMVA;
Qd_base = ps.shunt(:,C.sh.Q) .* ps.shunt(:,C.sh.factor)/ps.baseMVA;
% compute the load that comes from the const impedance factor
Pd_const_Z = Pd_base .* Vd.^2 .* ps.shunt(:,C.sh.frac_Z);
Qd_const_Z = Qd_base .* Vd.^2 .* ps.shunt(:,C.sh.frac_Z);
% compute the load that comes from the const current factor
Pd_const_I = Pd_base .* Vd .* (1-(ps.shunt(:,C.sh.frac_Z)+ps.shunt(:,C.sh.frac_S)+ps.shunt(:,C.sh.frac_E)));
Qd_const_I = Qd_base .* Vd .* (1-(ps.shunt(:,C.sh.frac_Z)+ps.shunt(:,C.sh.frac_S)+ps.shunt(:,C.sh.frac_E)));
% compute the const S load
Pd_const_S = Pd_base .* ps.shunt(:,C.sh.frac_S);
Qd_const_S = Qd_base .* ps.shunt(:,C.sh.frac_S);
% compute the exponential load
Pd_const_E = Pd_base .* Vd.^ps.shunt(:,C.sh.gamma) .* ps.shunt(:,C.sh.frac_E);
Qd_const_E = Qd_base .* Vd.^ps.shunt(:,C.sh.gamma) .* ps.shunt(:,C.sh.frac_E);

Pd_bus = sparse(ps.bus_i(ps.shunt(:,1)),1, Pd_const_Z + Pd_const_I + Pd_const_S + Pd_const_E, size(ps.bus,1),1);
Qd_bus = sparse(ps.bus_i(ps.shunt(:,1)),1, Qd_const_Z + Qd_const_I + Qd_const_S + Qd_const_E, size(ps.bus,1),1);

% calculate P and Q for the machines
Pgen = (Eaps.*mac_Vmags./Xdps).*sin(delta_m) +...
    (mac_Vmags.^2./2).*(1./Xqs-1./Xdps).*sin(2*delta_m); % eq. 7.81, Bergen & Vittal
Qgen = (Eaps.*mac_Vmags./Xdps).*cos(delta_m) - ...
    mac_Vmags.^2.*(cos(delta_m).^2./Xdps + sin(delta_m).^2./Xqs); % modified from eq. 6.46, Bergen & Vittal

mac_P = sparse(mac_bus_i,1,Pgen,n,1);
mac_Q = sparse(mac_bus_i,1,Qgen,n,1);

% output
g(ix.g.P)       = P_trans - mac_P + Pd_bus;
g(ix.g.Q)       = Q_trans - mac_Q + Qd_bus;
if ~angle_ref
    g(ix.g.slack)   = Thetas(is_slack) - Theta_slack;
end

% output dg_dx and dg_dy if requested
if nargout>1
    % real power & reactive power
    if ~angle_ref
        dPg_ddelta_values = (Eaps.*mac_Vmags./Xdps) .* cos(delta_m) + ...
            mac_Vmags.^2.*(1./Xqs-1./Xdps).*cos(2*delta_m);
        dQg_ddelta_values = (Eaps.*mac_Vmags./Xdps) .* (-sin(delta_m)) -....
            mac_Vmags.^2.*(-sin(2*delta_m)./Xdps + sin(2*delta_m)./Xqs);
    else
        dPg_ddelta_values = ((Eaps.*mac_Vmags./Xdps) .* cos(delta_m) + ...
            mac_Vmags.^2.*(1./Xqs-1./Xdps).*cos(2*delta_m))*ones(1,ng).*((eye(ng)-weight./sum(weight)*ones(1,ng))');
        dPg_ddelta_values = dPg_ddelta_values';
        dPg_ddelta_values = dPg_ddelta_values(:);
        
        dQg_ddelta_values = ((Eaps.*mac_Vmags./Xdps) .* (-sin(delta_m)) -...
            mac_Vmags.^2.*(-sin(2*delta_m)./Xdps + sin(2*delta_m)./Xqs))*ones(1,ng).*((eye(ng)-weight./sum(weight)*ones(1,ng))');
        dQg_ddelta_values = dQg_ddelta_values';
        dQg_ddelta_values = dQg_ddelta_values(:);
    end
    dPg_dEa_values    = (mac_Vmags./Xdps) .* sin(delta_m);
    dQg_dEa_values    = (mac_Vmags./Xdps) .* cos(delta_m);
    
    % assemble dg_dx
    dg_dx = sparse(ix.ny,ix.nx);    
    Pg_rows = ix.g.P(mac_bus_i);
    Qg_rows = ix.g.Q(mac_bus_i);
    if ~angle_ref
        Pg_rows_loc = Pg_rows;
        delta_loc = ix.x.delta;
        Qg_rows_loc = Qg_rows;
    else       
        [Pg_row,~] = meshgrid(Pg_rows, ix.x.delta);
        Pg_rows_loc = Pg_row(:);
        delta_loc = ix.COI.delta;
        [Qg_row,~] = meshgrid(Qg_rows, ix.x.delta);
        Qg_rows_loc = Qg_row(:);
    end
    % dPg_ddelta
    dg_dx = dg_dx + sparse(Pg_rows_loc,delta_loc,-dPg_ddelta_values,ix.ny,ix.nx);
    % dPg_dEa
    dg_dx = dg_dx + sparse(Pg_rows,ix.x.Eap,-dPg_dEa_values,ix.ny,ix.nx);
    % dQg_ddelta
    dg_dx = dg_dx + sparse(Qg_rows_loc,delta_loc,-dQg_ddelta_values,ix.ny,ix.nx);
    % dQg_dEa
    dg_dx = dg_dx + sparse(Qg_rows,ix.x.Eap,-dQg_dEa_values,ix.ny,ix.nx);
end
if nargout>2
    % assumes that the Ybus in ps is updated based on the discrete state of the system
    % find dP_IK_dVmag_I and K
    dP_IK_dVmag_I = Vmag_K .* (g_IK.*cos_theta_IK + b_IK.*sin_theta_IK);
    dP_IK_dVmag_K = Vmag_I .* (g_IK.*cos_theta_IK + b_IK.*sin_theta_IK);
    % find dQ_IK_dVmag_I and K
    dQ_IK_dVmag_I = Vmag_K .* (g_IK.*sin_theta_IK - b_IK.*cos_theta_IK);
    dQ_IK_dVmag_K = Vmag_I .* (g_IK.*sin_theta_IK - b_IK.*cos_theta_IK);
    % find dP_IK_dtheta_I and K
    dP_IK_dtheta_I = Vmag_IK .* (-g_IK.*sin_theta_IK + b_IK.*cos_theta_IK);
    dP_IK_dtheta_K = Vmag_IK .* (+g_IK.*sin_theta_IK - b_IK.*cos_theta_IK);
    % find dQ_IK_dtheta_I and K
    dQ_IK_dtheta_I = Vmag_IK .* (+g_IK.*cos_theta_IK + b_IK.*sin_theta_IK);
    dQ_IK_dtheta_K = Vmag_IK .* (-g_IK.*cos_theta_IK - b_IK.*sin_theta_IK);
    
    % assemble dg_dy
    dg_dy = sparse(ix.ny,ix.ny);
    % assemble dP_IK_dVmag_I and K
    dg_dy = dg_dy + sparse(ix.g.P(I),ix.y.Vmag(I),dP_IK_dVmag_I,ix.ny,ix.ny);
    dg_dy = dg_dy + sparse(ix.g.P(I),ix.y.Vmag(K),dP_IK_dVmag_K,ix.ny,ix.ny);
    % assemble dQ_IK_dVmag_I and K
    dg_dy = dg_dy + sparse(ix.g.Q(I),ix.y.Vmag(I),dQ_IK_dVmag_I,ix.ny,ix.ny);
    dg_dy = dg_dy + sparse(ix.g.Q(I),ix.y.Vmag(K),dQ_IK_dVmag_K,ix.ny,ix.ny);
    % assemble dP_IK_dtheta_I and K
    dg_dy = dg_dy + sparse(ix.g.P(I),ix.y.theta(I),dP_IK_dtheta_I,ix.ny,ix.ny);
    dg_dy = dg_dy + sparse(ix.g.P(I),ix.y.theta(K),dP_IK_dtheta_K,ix.ny,ix.ny);
    % assemble dQ_IK_dtheta_I and K
    dg_dy = dg_dy + sparse(ix.g.Q(I),ix.y.theta(I),dQ_IK_dtheta_I,ix.ny,ix.ny);
    dg_dy = dg_dy + sparse(ix.g.Q(I),ix.y.theta(K),dQ_IK_dtheta_K,ix.ny,ix.ny);
    
    % put the voltages back into the generator derivatives
    % differentiate Pgen = (Eaps.*mac_Vmags./Xdps).*sin(delta_m);
    dPgen_dVgen = Eaps ./ Xdps .* sin(delta_m) +...
        mac_Vmags.*(1./Xqs-1./Xdps).*sin(2*delta_m);
    % differentiate Qgen = (Eaps.*mac_Vmags./Xdps).*cos(delta_m) - (mac_Vmags.^2)./Xdps;
    dQgen_dVgen = Eaps ./ Xdps .* cos(delta_m) - ...
        2*mac_Vmags.*(cos(delta_m).^2./Xdps + sin(delta_m).^2./Xqs);
    % stick the stuff above into dg_dy
    dg_dy = dg_dy + sparse(ix.g.P(mac_bus_i),ix.y.Vmag(mac_bus_i),-dPgen_dVgen,ix.ny,ix.ny);
    dg_dy = dg_dy + sparse(ix.g.Q(mac_bus_i),ix.y.Vmag(mac_bus_i),-dQgen_dVgen,ix.ny,ix.ny);
    
    % now do the generator P/Q with respect to the generator bus angles
    dPg_dtheta_values = (Eaps.*mac_Vmags./Xdps) .* cos(delta_m) + ...
                         mac_Vmags.^2.*(1./Xqs-1./Xdps).*cos(2*delta_m);
    dQg_dtheta_values = (Eaps.*mac_Vmags./Xdps) .* (-sin(delta_m)) -...
                         mac_Vmags.^2.*(-sin(2*delta_m)./Xdps + sin(2*delta_m)./Xqs);
    dg_dy = dg_dy + sparse(ix.g.P(mac_bus_i),ix.y.theta(mac_bus_i),+dPg_dtheta_values,ix.ny,ix.ny);
    dg_dy = dg_dy + sparse(ix.g.Q(mac_bus_i),ix.y.theta(mac_bus_i),+dQg_dtheta_values,ix.ny,ix.ny);
    
    % fix the derivatives with [Z]IPE contributions
    dg_dy = dg_dy + sparse(ix.g.P(load_locs),ix.y.Vmag(load_locs), 2*Pd_base.*Vd.*ps.shunt(:,C.sh.frac_Z),ix.ny,ix.ny);
    dg_dy = dg_dy + sparse(ix.g.Q(load_locs),ix.y.Vmag(load_locs), 2*Qd_base.*Vd.*ps.shunt(:,C.sh.frac_Z),ix.ny,ix.ny);
    % fix the derivatives with Z[I]PE contributions
    dg_dy = dg_dy + sparse(ix.g.P(load_locs),ix.y.Vmag(load_locs), Pd_base.*(1-(ps.shunt(:,C.sh.frac_Z)+ps.shunt(:,C.sh.frac_S)+ps.shunt(:,C.sh.frac_E))),ix.ny,ix.ny);
    dg_dy = dg_dy + sparse(ix.g.Q(load_locs),ix.y.Vmag(load_locs), Qd_base.*(1-(ps.shunt(:,C.sh.frac_Z)+ps.shunt(:,C.sh.frac_S)+ps.shunt(:,C.sh.frac_E))),ix.ny,ix.ny);
    % fix the derivatives with ZIP[E] contributions
    dg_dy = dg_dy + sparse(ix.g.P(load_locs),ix.y.Vmag(load_locs), ps.shunt(:,C.sh.gamma).*Pd_base.*Vd.^(ps.shunt(:,C.sh.gamma)-1).*ps.shunt(:,C.sh.frac_E),ix.ny,ix.ny);
    dg_dy = dg_dy + sparse(ix.g.Q(load_locs),ix.y.Vmag(load_locs), ps.shunt(:,C.sh.gamma).*Qd_base.*Vd.^(ps.shunt(:,C.sh.gamma)-1).*ps.shunt(:,C.sh.frac_E),ix.ny,ix.ny);
    
    if ~angle_ref
        % fix the derivatives for the slack bus / delta_sys
        dg_dy = dg_dy + sparse(ix.y.delta_sys,ix.y.theta(is_slack),1,ix.ny,ix.ny);
        dg_dy = dg_dy + sparse(ix.g.P(mac_bus_i),ix.y.delta_sys,-dPg_ddelta_values,ix.ny,ix.ny);
        dg_dy = dg_dy + sparse(ix.g.Q(mac_bus_i),ix.y.delta_sys,-dQg_ddelta_values,ix.ny,ix.ny);
    end
end