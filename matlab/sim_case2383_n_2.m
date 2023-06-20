function [outputs,ps] = sim_case2383_n_2(a,b) 
%% simulate 2383-bus polish case
C = psconstants;
clear t_delay t_prev_check num_ls

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

% simulation time
t_max = 100;

% select data case to simulate
load case2383_mod_ps_dyn ps;
ps.branch(:,C.br.tap)       = 1;
ps.shunt(:,C.sh.factor)     = 1;
ps.shunt(:,C.sh.status)     = 1;
ps.shunt(:,C.sh.frac_S)     = 0;
ps.shunt(:,C.sh.frac_E)     = 1;
ps.shunt(:,C.sh.frac_Z)     = 0;
ps.shunt(:,C.sh.gamma)      = 0.08;

% fix the negative load in data, which causes problem for load shedding
if any(ps.shunt(:,C.sh.P)<0)
    neg_load = find(ps.shunt(:,C.sh.P)<0);
    ps.shunt(neg_load,C.sh.P:C.sh.Q) = ps.shunt(neg_load,C.sh.P:C.sh.Q)*(-1);
end
% fix the negative active generation for genartors
if any(ps.gen(:,C.ge.P)<0)
    neg_gen = ps.gen(:,C.ge.P)<0;
    ps.gen(neg_gen,C.ge.P) = 0;
end

% set some options
opt = psoptions2383;
opt.sim.integration_scheme = 1;
opt.sim.dt_default = 0.1;
opt.nr.use_fsolve = true;
opt.verbose = true;
opt.sim.writelog = true;
% opt.pf.linesearch = 'cubic_spline';
opt.sim.gen_control = 1;        % 0 = generator without exciter and governor, 1 = generator with exciter and governor
opt.sim.angle_ref = 0;          % 0 = delta_sys, 1 = center of inertia---delta_coi
                                % Center of inertia doesn't work when having islanding
opt.sim.COI_weight = 1;         % 1 = machine inertia, 0 = machine MVA base(Powerworld)
opt.sim.uvls_tdelay_ini = 0.5;  % 1 sec delay for uvls relay.
opt.sim.ufls_tdelay_ini = 0.5;  % 1 sec delay for ufls relay.
opt.sim.dist_tdelay_ini = 0.5;  % 1 sec delay for dist relay.
% Don't forget to change this value (opt.sim.time_delay_ini) in solve_dae.m

% initialize the case
ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
ps = update_load_freq_source(ps);

% build the machine variables
[ps.mac,ps.exc,ps.gov]      = get_mac_state(ps,'salient');
% initialize relays
ps.relay                    = get_relays(ps,'all',opt);

% initialize global variables
global t_delay t_prev_check dist2threshold state_a
n    = size(ps.bus,1);
ng   = size(ps.mac,1);
m    = size(ps.branch,1);
n_sh = size(ps.shunt,1);
ix   = get_indices(n,ng,m,n_sh,opt);
t_delay = inf(size(ps.relay,1),1);
t_delay([ix.re.uvls])= opt.sim.uvls_tdelay_ini;
t_delay([ix.re.ufls])= opt.sim.ufls_tdelay_ini;
t_delay([ix.re.dist])= opt.sim.dist_tdelay_ini;
t_prev_check = nan(size(ps.relay,1),1);
dist2threshold = inf(size(ix.re.oc,2)*2,1);
state_a = zeros(size(ix.re.oc,2)*2,1);

%% build an event matrix
event = zeros(3,C.ev.cols);
% start
event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
% trip a branch
event(2,[C.ev.time C.ev.type]) = [5 C.ev.trip_branch];
event(2,C.ev.branch_loc) = a;
% trip a branch
event(3,[C.ev.time C.ev.type]) = [5 C.ev.trip_branch];
event(3,C.ev.branch_loc) = b;
% set the end time
event(4,[C.ev.time C.ev.type]) = [t_max C.ev.finish];

%% run the simulation
[outputs,ps] = simgrid(ps,event,'sim_case2383',opt);

%% print the results
% fname = outputs.outfilename;
% [t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd] = read_outfile(fname,ps,opt);
% omega_0 = 2*pi*ps.frequency;
% omega_pu = omega / omega_0;
% 
% % delete *.csv; delete trace*;
% 
% figure(1); clf;
% subplot(4,1,1); hold on; 
% nl = size(omega_pu,2); colorset = varycolor(nl);
% set(gca, 'ColorOrder', colorset);
% plot(t',omega_pu');
% ylabel('omega pu');
% %legend(cellstr(num2str((1:nl)', 'omega_%d'))); legend boxoff; 
% 
% subplot(4,1,2); hold on; 
% nl = size(theta,2); colorset = varycolor(nl);
% set(gca, 'ColorOrder', colorset);
% plot(t',theta');
% ylabel('theta');
% %legend(cellstr(num2str((1:nl)', 'theta_%d'))); legend boxoff;
% 
% subplot(4,1,3); hold on; 
% nl = size(Vmag,2); colorset = varycolor(nl);
% set(gca, 'ColorOrder', colorset);
% plot(t',Vmag');
% ylabel('Vmag');
% xlabel('time');
% %legend(cellstr(num2str((1:nl)', 'Vmag_%d'))); legend boxoff;
% 
% subplot(4,1,4); hold on; 
% nl = size(delta,2); colorset = varycolor(nl);
% set(gca, 'ColorOrder', colorset);
% plot(t',delta');
% ylabel('Delta');
% xlabel('time');
% %legend(cellstr(num2str((1:nl)', 'Delta_%d'))); legend boxoff;
