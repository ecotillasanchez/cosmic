function [x,y] = get_xy_rec(ps,opt)
% usage: [x,y] = get_xy_rec(ps)
% this function gets the x and y vectors from the data in ps for
% rectangular formulation

% get some constants
C  = psconstants;
n  = size(ps.bus,1);
ng = size(ps.gen,1);
m  = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices_rec(n,ng,m,n_sh,opt);

angle_ref = opt.sim.angle_ref;                 % angle reference: 0:delta_sys,1:delta_coi
COI_weight = opt.sim.COI_weight;               % weight of center of inertia

if COI_weight
    weight      = ps.mac(:,C.ma.M);
else
    weight     = ps.gen(:,C.ge.mBase);
end
% get data from the ps structure
omega_0     = 2*pi*ps.frequency;
if ~angle_ref
    mac_bus_i   = ps.bus_i(ps.mac(:,1));    
    delta_sys   = ps.bus(1,C.bu.delta_sys);     % in case we were in islanded state
    Vr          = ps.bus(:,C.bu.Vr);
    Vi          = ps.bus(:,C.bu.Vi);
    theta       = ps.bus(:,C.bu.Vang)*pi/180;   % it is called in delta
    delta       = ps.mac(:,C.ma.delta_m) - delta_sys + theta(mac_bus_i);
else
    delta       = ps.mac(:,C.ma.delta);
    delta_coi = sum(weight.*delta)/sum(weight);          % Center of inertia
    theta       = ps.bus(:,C.bu.Vang)*pi/180-delta_coi;  % use COI as reference
    V           = ps.bus(:,C.bu.Vmag).*exp(1i.*theta);
    Vr          = real(V);
    Vi          = imag(V);
end

omega_pu    = ps.mac(:,C.ma.omega)/omega_0;
Pm          = ps.mac(:,C.ma.Pm);
Eap         = ps.mac(:,C.ma.Eap);
E1 			= ps.exc(:,C.ex.E1);
Efd 		= ps.exc(:,C.ex.Efd);
P3  		= ps.gov(:,C.go.P3);

% build x
x = zeros(ix.nx,1);
x(ix.x.delta)    = delta;
x(ix.x.omega_pu) = omega_pu;
x(ix.x.Pm)       = Pm;
x(ix.x.Eap)      = Eap;
x(ix.x.E1) 		 = E1;
x(ix.x.Efd) 	 = Efd;
x(ix.x.P3) 	     = P3;

% build y
y = zeros(ix.ny,1);
y(ix.y.Vr)          = Vr;
y(ix.y.Vi)          = Vi;
if ~angle_ref
    y(ix.y.delta_sys)   = delta_sys;
end
