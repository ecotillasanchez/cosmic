function [x,y] = get_xy(ps,opt)
% usage: [x,y] = get_xy(ps)
% this function gets the x and y vectors from the data in ps

% get some constants
C  = psconstants;
n  = size(ps.bus,1);
ng = size(ps.gen,1);
m  = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices(n,ng,m,n_sh,opt);

angle_ref = opt.sim.angle_ref;                 % angle reference: 0:delta_sys,1:delta_coi
COI_weight = opt.sim.COI_weight;               % weight of center of inertia

if COI_weight
    weight      = ps.mac(:,C.ma.M);
else
    weight     = ps.gen(:,C.ge.mBase);
end

% get data from the ps structure
omega_0     = 2*pi*ps.frequency;
Vmag        = ps.bus(:,C.bu.Vmag);

if ~angle_ref
    mac_bus_i   = ps.bus_i(ps.mac(:,1));
    delta_sys   = ps.bus(1,C.bu.delta_sys);     % in case we were in islanded state
    theta       = ps.bus(:,C.bu.Vang)*pi/180;
    delta       = ps.mac(:,C.ma.delta_m) - delta_sys + theta(mac_bus_i);
else
    delta       = ps.mac(:,C.ma.delta);
    delta_coi   = sum(weight.*delta)/sum(weight);        % Center of inertia
    theta       = ps.bus(:,C.bu.Vang)*pi/180-delta_coi;  % use COI as reference
end

omega_pu    = ps.mac(:,C.ma.omega)/omega_0;
Pm          = ps.mac(:,C.ma.Pm);
Eap         = ps.mac(:,C.ma.Eap);
E1 			= ps.exc(:,C.ex.E1);
Efd 		= ps.exc(:,C.ex.Efd);
P3  		= ps.gov(:,C.go.P3);
Temperature = ps.relay(ix.re.temp,C.re.state_a);

% build x
x = zeros(ix.nx,1);
x(ix.x.delta)    = delta;
x(ix.x.omega_pu) = omega_pu;
x(ix.x.Pm)       = Pm;
x(ix.x.Eap)      = Eap;
x(ix.x.E1) 		 = E1;
x(ix.x.Efd) 	 = Efd;
x(ix.x.P3) 	     = P3;
x(ix.x.temp)     = Temperature;

% build y
y = zeros(ix.ny,1);
y(ix.y.Vmag)        = Vmag;
y(ix.y.theta)       = theta;
if ~angle_ref
    y(ix.y.delta_sys)   = delta_sys;
end
