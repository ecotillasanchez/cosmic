%%  test discrete change
clear all; close all; clc;
addpath('../data');
C  = psconstants;
opt = psoptions;
ps = case39_ps;

ps = newpf(ps);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
% build indices
n  = size(ps.bus,1);
ng = size(ps.gen,1);
m  = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices(n,ng,m,n_sh,opt);

% build x and y
[x0,y0] = get_xy(ps,opt);
xy0     = [x0;y0];

algebraic_eqs_only([],x0,y0,ps,opt)

% inject a perturbation to the algebraic variables and resolve
y1 = y0 + 0.01*randn(size(y0));
%x0(1) = x0(1) + 0.001;
y1p = solve_algebraic([],x0,y1,ps,opt);

return
% implement the discrete event and update ps structure
ps.branch(5,C.br.status) = 0;
ps = runPowerFlow(ps);
ps.Ybus = getYbus(ps,false);
ps.mac  = get_mac_state(ps,'salient');
[x1,y1] = get_xy(ps);
xy1     = [x1;y1];
algebraic_eqs_only([],x1,y1,ps)

g_handle = @(y) algebraic_eqs_only([],x1,y,ps,opt);
checkDerivatives(g_handle,[],y1);
return

ps.mac(:,C.mac.delta)   = x1(ix.x.delta) - y1(ix.y.theta(mac_bus_i));
ps.bus(:,C.bus.Vmag)    = y1(ix.y.Vmag); 
ps.bus(:,C.bus.Vang)    = y1(ix.y.theta)/pi*180;

% solve algebraic vars and integrate within-change time slot
y1p = solve_algebraic(t1,x1,y1,ps);