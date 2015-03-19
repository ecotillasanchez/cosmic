function [f,df_dx,df_dy] = differential_eqs(t,x,y,ps,opt) %#ok<INUSL>
% usage: [f,df_dx,df_dy] = differential_eqs(t,x,y,ps,opt)
% differential equations that model the elements of the power system
%
% inputs:
%   t   -> time in seconds
%   x   -> [delta omega Pm Eap E1 Efd] for each machine
%   y   -> [Vmag Theta]
%   ps  -> power system structure
%   opt -> options structure
%
% outputs: 
%   f(1) -> ddelta_dt
%   f(2) -> domega_dt
%   f(3) -> dPm_dt
%   f(4) -> dEap_dt
%   f(5) -> dE1_dt
%   f(6) -> dEfd_dt
%   f(7) -> dP3_dt

% constants
C           = psconstants;
n           = size(ps.bus,1);
ng          = size(ps.mac,1);
m           = size(ps.branch,1);
n_sh        = size(ps.shunt,1);
ix          = get_indices(n,ng,m,n_sh,opt);
gc          = opt.sim.gen_control; 
angle_ref = opt.sim.angle_ref;                 % angle reference: 0:delta_sys,1:delta_coi
COI_weight = opt.sim.COI_weight;               % weight of center of inertia

% extract parameters from ps
Xds     = ps.mac(:,C.ma.Xd);
Xdps    = ps.mac(:,C.ma.Xdp);
Xqs     = ps.mac(:,C.ma.Xq);
Td0ps   = ps.mac(:,C.ma.Td0p);
Ds      = ps.mac(:,C.ma.D);
Ms      = ps.mac(:,C.ma.M);
omega_0 = 2*pi*ps.frequency;

if COI_weight
    weight      = Ms;
else
    weight     = ps.gen(:,C.ge.mBase);
end

% extract differential variables
deltas      = x(ix.x.delta);
omegas_pu   = x(ix.x.omega_pu);
Pms         = x(ix.x.Pm);
Eaps        = x(ix.x.Eap);
Efds        = x(ix.x.Efd);
E1s         = x(ix.x.E1);
P3s         = x(ix.x.P3);

% extract algebraic variables and fix the slack bus angle
mac_buses   = ps.mac(:,C.mac.gen);
mac_bus_i   = ps.bus_i(mac_buses);
Vmags       = y(ix.y.Vmag);
Thetas      = y(ix.y.theta);
mac_Vmags   = Vmags(mac_bus_i);
mac_Thetas  = Thetas(mac_bus_i);

% machine angles, relative to the bus angles:
if ~angle_ref
    delta_sys = y(ix.y.delta_sys);
    delta_m = deltas + delta_sys - mac_Thetas;
else
    delta_coi = sum(weight.*deltas)/sum(weight);
    delta_m = deltas - delta_coi - mac_Thetas;
    omegas_coi = sum(weight.*omegas_pu,1)/sum(weight,1);
end


% calculate Pe: Eq. 7.81 from Bergen & Vittal
Pes = (Eaps.*mac_Vmags./Xdps).*sin(delta_m) + mac_Vmags.^2./2.*(1./Xqs-1./Xdps).*sin(2*delta_m);

% initialize output
f = zeros(length(x),1);

% calculate swing equations
if ~angle_ref
    f(ix.f.delta_dot) = omega_0.*(omegas_pu-1);
    f(ix.f.omega_dot) = (Pms - Pes - Ds.*(omegas_pu-1))./Ms;
else
    f(ix.f.delta_dot) = omega_0.*(omegas_pu-omegas_coi);
    f(ix.f.omega_dot) = (Pms - Pes - Ds.*(omegas_pu-omegas_coi))./Ms;
end
f(ix.f.Eap_dot)   = -Eaps.*Xds./(Td0ps.*Xdps)+(Xds./Xdps-1).*mac_Vmags.*cos(delta_m)./Td0ps+Efds./Td0ps; % Eq. 7.75 from Bergen & Vittal   

% calculate governor and exciter equations
if gc
    [f(ix.f.Pm_dot),f(ix.f.P3_dot),df_dx_gov]               = governor_eqs_modified([Pms';P3s'],omegas_pu,ps);
    [f(ix.f.Efd_dot),f(ix.f.E1_dot),df_dx_exc,df_dy_exc]   	= exciter_eqs([Efds';E1s'],mac_Vmags,ps);
end

% output df_dx and df_dy if requested
if nargout>1
    % build df_dx
    if ~angle_ref
        dPg_ddelta = (Eaps.*mac_Vmags./Xdps).*cos(delta_m) + mac_Vmags.^2.*(1./Xqs-1./Xdps).*cos(2*delta_m);
        dFswing_ddelta_values = -dPg_ddelta./Ms;
        dFswing_domega_values = -Ds./Ms;
        dFdelta_dot_domega    = omega_0;
        dEap_ddelta_values    = -(Xds./Xdps-1).*mac_Vmags.*sin(delta_m)./Td0ps; 
    else
        dPg_ddelta = ((Eaps.*mac_Vmags./Xdps).*cos(delta_m) + mac_Vmags.^2.*(1./Xqs-1./Xdps).*cos(2*delta_m))*ones(1,ng).*((eye(ng)-weight./sum(weight)*ones(1,ng))');
        dFswing_ddelta_values = -dPg_ddelta./(Ms*ones(1,ng));
        dFswing_ddelta_values = dFswing_ddelta_values';
        dFswing_ddelta_values = dFswing_ddelta_values(:);
        
        dFswing_domega_values = -Ds*ones(1,ng).*(eye(ng)-weight./sum(weight)*ones(1,ng))'./(Ms*ones(1,ng));
        dFswing_domega_values = dFswing_domega_values';
        dFswing_domega_values = dFswing_domega_values(:);
        
        dFdelta_dot_domega    = omega_0.*(eye(ng)-weight./sum(weight)*ones(1,ng))';
        dFdelta_dot_domega    = dFdelta_dot_domega';
        dFdelta_dot_domega    = dFdelta_dot_domega(:);
        
        dEap_ddelta_values    = -(Xds./Xdps-1).*mac_Vmags.*sin(delta_m)./Td0ps*ones(1,ng).*((eye(ng)-weight./sum(weight)*ones(1,ng))');
        dEap_ddelta_values    = dEap_ddelta_values';
        dEap_ddelta_values    = dEap_ddelta_values(:);
    end

    dFswing_dPm_values    = 1./Ms;
    dFswing_dEa_values    = -sin(delta_m).*mac_Vmags./(Ms.*Xdps);
    dEap_dEap_values      = -Xds./(Td0ps.*Xdps);    
    dEap_dEfd_values      = 1./Td0ps;
    if gc
        dE1_dE1_values         = df_dx_exc(:,1);
        dEfd_dE1_values        = df_dx_exc(:,2);
        dEfd_dEfd_values       = df_dx_exc(:,3);
        dPm_domegas_values     = df_dx_gov(:,1);
        dPm_dPm_values         = df_dx_gov(:,2);
        dPm_dP3_values         = df_dx_gov(:,3);
        dP3_dP3_values         = df_dx_gov(:,4);
        dP3_domegas_values     = df_dx_gov(:,5);
    end

    if ~angle_ref
        omega_dot_loc = ix.f.omega_dot;
        delta_loc = ix.x.delta;
        delta_dot_loc = ix.f.delta_dot;
        omega_pu_loc = ix.x.omega_pu;
        Eap_dot_loc = ix.f.Eap_dot;
    else
        omega_dot_loc = ix.COI.omega_dot;
        delta_loc = ix.COI.delta;
        delta_dot_loc = ix.COI.delta_dot;
        omega_pu_loc = ix.COI.omega_pu;
        Eap_dot_loc = ix.COI.Eap_dot;
    end
    % assemble df_dx
    df_dx = sparse(ix.nx,ix.nx);
    % dFswing_ddelta
    df_dx = df_dx + sparse(omega_dot_loc,delta_loc,dFswing_ddelta_values,ix.nx,ix.nx);
    % dFswing_domega
    df_dx = df_dx + sparse(omega_dot_loc,omega_pu_loc,dFswing_domega_values,ix.nx,ix.nx);
    % dFswing_dPm
    df_dx = df_dx + sparse(ix.f.omega_dot,ix.x.Pm,dFswing_dPm_values,ix.nx,ix.nx);
    % dFswing_dEa
    df_dx = df_dx + sparse(ix.f.omega_dot,ix.x.Eap,dFswing_dEa_values,ix.nx,ix.nx);
    % dFdelta_dot_domega
    df_dx = df_dx + sparse(delta_dot_loc,omega_pu_loc,dFdelta_dot_domega,ix.nx,ix.nx);
	% dEap_dot_dEap
    df_dx = df_dx + sparse(ix.f.Eap_dot,ix.x.Eap,dEap_dEap_values,ix.nx,ix.nx);
    % dEap_dot_ddelta
    df_dx = df_dx + sparse(Eap_dot_loc,delta_loc,dEap_ddelta_values,ix.nx,ix.nx);
    % dEap_dot_dEfd
    df_dx = df_dx + sparse(ix.f.Eap_dot,ix.x.Efd,dEap_dEfd_values,ix.nx,ix.nx);
    if gc
        % dE1_dot_dE1
        df_dx = df_dx + sparse(ix.f.E1_dot,ix.x.E1,dE1_dE1_values,ix.nx,ix.nx);
        % dEfd_dot_dE1
        df_dx = df_dx + sparse(ix.f.Efd_dot,ix.x.E1,dEfd_dE1_values,ix.nx,ix.nx);
        % dEfd_dot_dEfd
        df_dx = df_dx + sparse(ix.f.Efd_dot,ix.x.Efd,dEfd_dEfd_values,ix.nx,ix.nx);
        % dPm_domega
        df_dx = df_dx + sparse(ix.f.Pm_dot,ix.x.omega_pu,dPm_domegas_values,ix.nx,ix.nx);
        % dPm_dot_dPm
        df_dx = df_dx + sparse(ix.f.Pm_dot,ix.x.Pm,dPm_dPm_values,ix.nx,ix.nx);
        % dPm_dP3
        df_dx = df_dx + sparse(ix.f.Pm_dot,ix.x.P3,dPm_dP3_values,ix.nx,ix.nx);
        % dP3_dot_dP3
        df_dx = df_dx + sparse(ix.f.P3_dot,ix.x.P3,dP3_dP3_values,ix.nx,ix.nx);
        % dP3_dot_domegas
        df_dx = df_dx + sparse(ix.f.P3_dot,ix.x.omega_pu,dP3_domegas_values,ix.nx,ix.nx);
    end
end
if nargout>2
    % build df_dy (change in f wrt the algebraic variables)
    dPg_dVmag = Eaps.*sin(delta_m)./Xdps + mac_Vmags.*(1./Xqs-1./Xdps).*sin(2*delta_m);
    dFswing_dVmag_values        = -(dPg_dVmag)./Ms;
	dEap_dVmag_values           = (Xds./Xdps-1).*cos(delta_m)./Td0ps;
    if gc
        dE1_dVmag_values            = df_dy_exc(:,1);
        dEfd_dVmag_values           = df_dy_exc(:,2);
    end
    
	% assemble df_dy
    cols = ix.y.Vmag(mac_bus_i);
    df_dy = sparse(ix.f.omega_dot,cols,dFswing_dVmag_values,ix.nx,ix.ny);
    
    if ~angle_ref
        dFswing_dtheta_values = dFswing_ddelta_values;
        cols = ix.y.delta_sys;
        % domega_dot_ddelta_sys
        df_dy = df_dy + sparse(ix.f.omega_dot,cols,dFswing_ddelta_values,ix.nx,ix.ny);
        %dEap_dot_ddelta_sys
        cols = ix.y.delta_sys;
        df_dy = df_dy + sparse(ix.f.Eap_dot,cols,dEap_ddelta_values,ix.nx,ix.ny);
    else
        dFswing_dtheta_values = -((Eaps.*mac_Vmags./Xdps).*cos(delta_m) + mac_Vmags.^2.*(1./Xqs-1./Xdps).*cos(2*delta_m))./Ms;
    end
    cols = ix.y.theta(mac_bus_i);
    df_dy = df_dy + sparse(ix.f.omega_dot,cols,-dFswing_dtheta_values,ix.nx,ix.ny);
	
	% dEap_dot_dVmag
    cols  = ix.y.Vmag(mac_bus_i);
    df_dy = df_dy + sparse(ix.f.Eap_dot,cols,dEap_dVmag_values,ix.nx,ix.ny);
    % dEap_dot_dtheta
    cols = ix.y.theta(mac_bus_i);
    dEap_dtheta_values = -(Xds./Xdps-1).*mac_Vmags.*sin(delta_m)./Td0ps;
    df_dy = df_dy + sparse(ix.f.Eap_dot,cols,-dEap_dtheta_values,ix.nx,ix.ny);

    if gc
        % dE1_dot_dVmag
        cols  = ix.y.Vmag(mac_bus_i);
        df_dy = df_dy + sparse(ix.f.E1_dot,cols,dE1_dVmag_values,ix.nx,ix.ny);
        % dEfd_dot_dVmag
        cols  = ix.y.Vmag(mac_bus_i);
        df_dy = df_dy + sparse(ix.f.Efd_dot,cols,dEfd_dVmag_values,ix.nx,ix.ny);
    end
end
