%% simulate 2383-bus polish case
clear all; close all; clc; C = psconstants;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

% simulation time
t_max = 30;

% select data case to simulate
load case2383_mod_ps_dyn;
ps.branch(:,C.br.tap)       = 1;
ps.shunt(:,C.sh.factor)     = 1;
ps.shunt(:,C.sh.status)     = 1;
ps.shunt(:,C.sh.frac_S)     = 0;
ps.shunt(:,C.sh.frac_E)     = 1;
ps.shunt(:,C.sh.frac_Z)     = 0;
ps.shunt(:,C.sh.gamma)      = 0.08;

% to differentiate the line MVA ratings
rateB_rateA                     = ps.branch(:,C.br.rateB)./ps.branch(:,C.br.rateA);
rateC_rateA                     = ps.branch(:,C.br.rateC)./ps.branch(:,C.br.rateA);
ps.branch(rateB_rateA==1,C.br.rateB)    = 1.1 * ps.branch(rateB_rateA==1,C.br.rateA);
ps.branch(rateC_rateA==1,C.br.rateC)    = 1.5 * ps.branch(rateC_rateA==1,C.br.rateA);

ps = updateps(ps);

% set the outage to simulate
branch_outages = [15,106]; % [15, 106];

% initialize the case
opt = psoptions2383;
opt.sim.integration_scheme = 1;
opt.sim.dt_default = 1/30;
opt.sim.writelog = true;
opt.nr.use_fsolve = true;
opt.verbose = true;
opt.sim.gen_control = 1;        % 0 = generator without exciter and governor, 1 = generator with exciter and governor
opt.sim.angle_ref = 0;          % 0 = delta_sys, 1 = center of inertia---delta_coi
                                % Center of inertia doesn't work when having islanding
opt.sim.COI_weight = 1;         % 1 = machine inertia, 0 = machine MVA base(Powerworld)
opt.sim.uvls_tdelay_ini = 0.5;  % 1 sec delay for uvls relay.
opt.sim.ufls_tdelay_ini = 0.5;  % 1 sec delay for ufls relay.
opt.sim.dist_tdelay_ini = 0.5;  % 1 sec delay for dist relay.
opt.sim.temp_tdelay_ini = 0;    % 0 sec delay for temp relay.
% Don't forget to change this value (opt.sim.time_delay_ini) in solve_dae.m

% initialize the case
ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
ps = update_load_freq_source(ps);
% build the machine variables
[ps.mac,ps.exc,ps.gov]      = get_mac_state(ps,'salient');
% initialize relays
ps.relay                    = get_relays(ps,'all',opt);

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
t_delay([ix.re.temp])= opt.sim.temp_tdelay_ini;
t_prev_check = nan(size(ps.relay,1),1);
dist2threshold = inf(size(ix.re.oc,2)*2,1);
state_a = zeros(size(ix.re.oc,2)*2,1);

%% build an event matrix with branch outages
n_event = length(branch_outages)+2;
event = zeros(n_event,C.ev.cols);
% start
event(1,[C.ev.time C.ev.type]) = [0 C.ev.start];
% trip branches
for i = 1:length(branch_outages)
    event(i+1,[C.ev.time C.ev.type]) = [10 C.ev.trip_branch];
    event(i+1,C.ev.branch_loc) = branch_outages(i);
end
% set the end time
event(end,[C.ev.time C.ev.type]) = [t_max C.ev.finish];

%% run the simulation
[outputs,ps] = simgrid(ps,event,'sim_case2383',opt);

%% print the results
if opt.sim.writelog 
    fname = outputs.outfilename;
    [t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd,P3,Temperature] = read_outfile(fname,ps,opt);
    omega_0 = 2*pi*ps.frequency;
    omega_pu = omega / omega_0;
   
end