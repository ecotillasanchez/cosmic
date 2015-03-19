%% driver that tests the f and g equations
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
end

% select data case to test steady state
ps = case9_ps;
%load case3120sp_ps_dyn;

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
ix   = get_indices(n,ng,m,n_sh);

% build x and y
[x0,y0] = get_xy(ps);
xy0     = [x0;y0];

% check derivatives for differential equations
disp('checking differential equations...');
f_handle = @(x) differential_eqs([],x,y0,ps,opt);
checkDerivatives(f_handle,[],x0);
% check f0
f0 = f_handle(x0);
max_f0 = max(abs(f0));
fprintf('max abs f = %e\n',max_f0);

% check derivatives for algebraic equations
disp('checking algebraic equations...');
g_handle = @(y) algebraic_eqs_only([],x0,y,ps,opt);
checkDerivatives(g_handle,[],y0);
% check g0
g0 = g_handle(y0);
max_g0 = max(abs(g0));
fprintf('max abs g = %e\n',max_g0);

% prepare the event handle function for relays
event_handle = @(t,xy) endo_event(t,xy,ix,ps);

% prepare the inputs for the ode solver
fn_fg   = @(t,xy) differential_algebraic(t,xy,ix.nx,ps,opt);
mass_matrix = sparse(1:ix.nx,1:ix.nx,1,ix.nx+ix.ny,ix.nx+ix.ny);
options = odeset(   'Mass',mass_matrix, ...
                    'MassSingular','yes', ...
                    'Jacobian', @(t,xy) get_jacobian(t,xy,ix.nx,ps,opt), ...
                    'Events', event_handle, ...
                    'Stats','on', ... 
                    'NormControl','off');

tic
odeout = ode15s(fn_fg,[0 500],xy0,options);
plot_dyn_results(odeout.x,odeout.y,ix);
toc;
