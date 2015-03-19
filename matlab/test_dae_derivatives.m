%% driver that tests the f and g equations for 9-bus case
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
end

% select data case
%load case2383_mod_ps_dyn;
ps  = updateps(case68_ps);

% initialize the case
opt = psoptions;
opt.nr.use_fsolve = true;
opt.verbose = true;
opt.sim.gen_control = 1;        % 0 = generator without exciter and governor, 1 = generator with exciter and governor
opt.sim.angle_ref = 0;          % 0 = delta_sys, 1 = center of inertia---delta_coi
                                % Now, center of inertia doesn't work when having islanding
opt.sim.COI_weight = 0;  

ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);

% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
% initialize relays
ps.relay    				= get_relays(ps,'all');

% build indices
n = size(ps.bus,1);
ng = size(ps.gen,1);
m = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices(n,ng,m,n_sh,opt);

% build x and y
[x0,y0] = get_xy(ps,opt);

%% check derivatives for differential equations
%% check df_dx
x = x0; y = y0;
disp('Checking df_dx for x0');
[f0,df_dx,df_dy] = differential_eqs(0,x,y,ps,opt);
f_handle = @(x) differential_eqs(0,x,y,ps,opt);
checkDerivatives(f_handle,df_dx,x);
% repeat from a perturbed starting point
disp('Checking df_dx for perturbed x0');
x = x0+randn(size(x0))*0.1;
[f0,df_dx,df_dy] = differential_eqs(0,x,y,ps,opt);
checkDerivatives(f_handle,df_dx,x);

%% check df_dy
x = x0; y = y0;
disp('Checking df_dy for y0');
[f0,df_dx,df_dy] = differential_eqs(0,x,y,ps,opt);
f_handle = @(y) differential_eqs(0,x,y,ps,opt);
checkDerivatives(f_handle,df_dy,y);
% repeat from a perturbed starting point
disp('Checking df_dy for perturbed y0');
y = y0+randn(size(y0))*0.1;
[f0,df_dx,df_dy] = differential_eqs(0,x,y,ps,opt);
checkDerivatives(f_handle,df_dy,y);

%% check dg_dx
x = x0; y = y0;
disp('Checking dg_dx for x0');
[f0,dg_dx,dg_dy] = algebraic_eqs(0,x,y,ps,opt);
g_handle = @(x) algebraic_eqs(0,x,y,ps,opt);
checkDerivatives(g_handle,dg_dx,x);
% repeat from a perturbed starting point
disp('Checking dg_dx for perturbed x0');
x = x0+randn(size(x0))*0.1;
[f0,dg_dx,dg_dy] = algebraic_eqs(0,x,y,ps,opt);
checkDerivatives(g_handle,dg_dx,x);

%% check dg_dy
x = x0; y = y0;
disp('Checking dg_dy for y0');
[f0,dg_dx,dg_dy] = algebraic_eqs(0,x,y,ps,opt);
g_handle = @(y) algebraic_eqs(0,x,y,ps,opt);
checkDerivatives(g_handle,dg_dy,y);
% repeat from a perturbed starting point
disp('Checking dg_dy for perturbed y0');
y = y0+randn(size(y0))*0.1;
[f0,dg_dx,dg_dy] = algebraic_eqs(0,x,y,ps,opt);
checkDerivatives(g_handle,dg_dy,y);

return
