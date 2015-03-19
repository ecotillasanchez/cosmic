%% driver that tests the f and g equations
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
end

% select data case to test steady state
%ps = updateps(case9_ps);
load case2383_mod_ps_dyn;

% set some options
opt = psoptions;
opt.sim.dt_default = 0.1;
opt.nr.use_fsolve = true;
opt.verbose = true;
% initialize the case
ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);

% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
% initialize relays
opt.sim.uvls_limit = 0.9;
ps.relay                    = get_relays(ps,'all',opt);

% build indices
n  = size(ps.bus,1);
ng = size(ps.gen,1);
m  = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices_i(n,ng,m,n_sh);

% build xy
[xy0] = get_xy_i(ps);
xyp0  = zeros(size(xy0));

% check derivatives for DAE equations
disp('checking derivatives...');
F_handle = @(xy) differential_algebraic_i([],xy,xyp0,ps,opt);
checkDerivatives(F_handle,[],xy0);

% prepare the inputs for the ode solver
fn_F   = @(t,xy,xyp) differential_algebraic_i(t,xy,xyp,ps,opt);
options = odeset(   'Jacobian', @(t,xy,xyp) get_jacobian_i(t,xy,xyp,ps,opt), ...
                    'Stats','on', ... 
                    'NormControl','off');

tic
[tout,xyout] = ode15i(fn_F,[0 500],xy0,xyp0,options);
plot_dyn_results_i(tout,xyout,ix);
toc;
