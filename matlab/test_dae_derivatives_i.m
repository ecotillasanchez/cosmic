%% driver that tests the f and g equations for 9-bus case
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
end

% select data case to test derivatives
ps = updateps(case9_ps);
%load case2383_mod_ps_dyn;

% initialize the case
opt = psoptions;
opt.nr.use_fsolve = true;
opt.verbose = true;
opt.sim.gen_control = 1;        % 0 = generator without exciter and governor, 1 = generator with exciter and governor
opt.sim.angle_ref  = 0;         % 0 = delta_sys, 1 = center of inertia---delta_coi
                                % Now, center of inertia doesn't work when having islanding
opt.sim.COI_weight = 0;  
opt.sim.dyn_load   = 0;         % 0 = ZIPE load; 1 = dynamic load

ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);

% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
% initialize relays
ps.relay    				= get_relays(ps,'all');
ps.shunt                    = get_load_state(ps,'ZIPE_only'); % dynamic or ZIPE_only

% build indices
n = size(ps.bus,1);
ng = size(ps.gen,1);
m = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices(n,ng,m,n_sh,opt);

% build x and y
[x0,y0] = get_xy(ps,opt);
keyboard
xy0     = [x0;y0];
xyp0    = zeros(size([x0;y0]));

%% check derivatives for DAE equations
%% check dF_dxy

disp('Checking dF_dxy for x0');
F_handle = @(xy) differential_algebraic_i([],xy,xyp0,ps,opt);
checkDerivatives(F_handle,[],xy0);

xy1 = xy0+randn(size(xy0))*0.1;
disp('Checking dF_dxy for perturbed x0');
F_handle = @(xy) differential_algebraic_i([],xy,xyp0,ps,opt);
checkDerivatives(F_handle,[],xy1);

return
