%% simulate RTS-96 case with topology known attack
clearvars; clear global; close all; clc; C = psconstants;

% load('mid_size_events_2.mat');rec_events = rec;
% load('rec_events.mat')
load rts96_2;
load_before_disturb = ps.shunt;

total_losses = NaN(1,3);
% for itr = 1:size(rec_events,2)
itr = 1;
total_losses(itr,1) = 0;

% do not touch path if we are deploying code
if ~(ismcc || isdeployed)
    addpath('../data');
    addpath('../numerics');
end

% simulation time
t_max = 130;
load rts96_2;

% load RTS_Demand;
% % 
% ps.shunt(:,C.sh.P)= RTS_Demand(:,2);
% select data case to simulate

ps = updateps(ps);
ps.branch(:,6:8) =ps.branch(:,6:8).*3;
% ps = replicate_case(ps,2);          
ps = unify_generators(ps); 
ps.branch(:,C.br.tap)       = 1;
ps.shunt(:,C.sh.factor)     = 1;
ps.shunt(:,C.sh.status)     = 1;
% ps.shunt(:,C.sh.frac_S)     = 1;
% ps.shunt(:,C.sh.frac_E)     = 0;
% ps.shunt(:,C.sh.frac_Z)     = 1;
% ps.shunt(:,C.sh.gamma)      = 0.08;

[ZIP_Coff_P_Q_Final] = ZIP_Coeff();
for qq=3:3:size(ps.shunt,1)
    % real power for residential
    ps.shunt(qq,C.sh.frac_S)     = ZIP_Coff_P_Q_Final.Pp_2(1);
    ps.shunt(qq,C.sh.frac_Z)     = ZIP_Coff_P_Q_Final.Zp_2(1);
    % reactive power for residential
    ps.shunt(qq,C.sh.frac_Q_Z)     = ZIP_Coff_P_Q_Final.Zq_2(1);
    ps.shunt(qq,C.sh.frac_Q_I)     = ZIP_Coff_P_Q_Final.Iq_2(1);
    ps.shunt(qq,C.sh.frac_Q_S)     = ZIP_Coff_P_Q_Final.Pq_2(1);
    
    if qq+1<=size(ps.shunt,1)
    % real power for Commercial
        ps.shunt(qq+1,C.sh.frac_S)     = ZIP_Coff_P_Q_Final.Pp_2(2);
        ps.shunt(qq+1,C.sh.frac_Z)     = ZIP_Coff_P_Q_Final.Zp_2(2);
    % reactive power for Commercial
        ps.shunt(qq+1,C.sh.frac_Q_Z)     = ZIP_Coff_P_Q_Final.Zq_2(2);
        ps.shunt(qq+1,C.sh.frac_Q_I)     = ZIP_Coff_P_Q_Final.Iq_2(2);
        ps.shunt(qq+1,C.sh.frac_Q_S)     = ZIP_Coff_P_Q_Final.Pq_2(2);
    end
    % real power for Industrial
     if qq+2<=size(ps.shunt,1)
        ps.shunt(qq+2,C.sh.frac_S)     = ZIP_Coff_P_Q_Final.Pp_2(3);
        ps.shunt(qq+2,C.sh.frac_Z)     = ZIP_Coff_P_Q_Final.Zp_2(3);
    % reactive power for Industrial
        ps.shunt(qq+2,C.sh.frac_Q_Z)     = ZIP_Coff_P_Q_Final.Zq_2(3);
        ps.shunt(qq+2,C.sh.frac_Q_I)     = ZIP_Coff_P_Q_Final.Iq_2(3);
        ps.shunt(qq+2,C.sh.frac_Q_S)     = ZIP_Coff_P_Q_Final.Pq_2(3);
     end
   
end

lembda_list=[0.1:0.1:1.2];

lembda=lembda_list(10);
% lembda=1.23;
ps.shunt(:,C.sh.P)     = ps.shunt(:,C.sh.P).*lembda;

limit_gen = find(ps.gen(:,C.ge.Pg)*lembda>=ps.gen(:,C.ge.Pmax));
for g = 1:size(limit_gen,1)
    ps.gen(limit_gen(g),C.ge.Pg)    =   ps.gen(limit_gen(g),C.ge.Pmax);
    ps.gov(limit_gen(g),C.go.Pmax)  =   ps.gov(limit_gen(g),C.go.Pmax).*1.1;
end

limit_gen_min = find(ps.gen(:,C.ge.Pg)*lembda<ps.gen(:,C.ge.Pmax));
for gg = 1:size(limit_gen_min,1)
    ps.gen(limit_gen_min(gg),C.ge.Pg)    =   ps.gen(limit_gen_min(gg),C.ge.Pg).*lembda;
    ps.gov(limit_gen_min(gg),C.go.Pmax)  = ps.gov(limit_gen_min(gg),C.go.Pmax).*lembda.*1.1;
end

% to differentiate the line MVA ratings
rateB_rateA                     = ps.branch(:,C.br.rateB)./ps.branch(:,C.br.rateA);
rateC_rateA                     = ps.branch(:,C.br.rateC)./ps.branch(:,C.br.rateA);
ps.branch(rateB_rateA==1,C.br.rateB)    = 1.1 * ps.branch(rateB_rateA==1,C.br.rateA);
ps.branch(rateC_rateA==1,C.br.rateC)    = 1.5 * ps.branch(rateC_rateA==1,C.br.rateA);

% set some options
opt = psoptions;
opt.sim.integration_scheme = 1;
opt.sim.dt_default = 1/10;
opt.nr.use_fsolve = true;
opt.sim.var_step = true;        % fixed (false) or variable (true) integration step size
% opt.pf.linesearch = 'cubic_spline';
opt.verbose = true;
opt.sim.gen_control = 1;        % 0 = generator without exciter and governor, 1 = generator with exciter and governor
opt.sim.angle_ref = 0;          % 0 = delta_sys, 1 = center of inertia---delta_coi
                                % Center of inertia doesn't work when having islanding
opt.sim.COI_weight = 0;         % 1 = machine inertia, 0 = machine MVA base(Powerworld)
opt.sim.oc_limit = 1.75;
opt.sim.uvls_tdelay_ini = 0.2;  % 1 sec delay for uvls relay.
opt.sim.ufls_tdelay_ini = 0.2;  % 1 sec delay for ufls relay.
opt.sim.dist_tdelay_ini = 0.5;  % 1 sec delay for dist relay.
opt.sim.temp_tdelay_ini = 0.1;    % 0 sec delay for temp relay.
opt.sim.oc_tdelay_ini = 0;  % 1 sec delay for oc relay.

% other relays' settings
opt.sim.uvls_limit = 0.87;     % 0.87 threshold at which under voltage load shedding occurs
opt.sim.uvls_delta = 0.25;     % the fraction of load that is shed during simulation if under voltage
opt.sim.ufls_limit = 0.985;    % 0.975 threshold at which under frequency load shedding occurs
opt.sim.ufls_delta = 0.25;     % load shedding fraction
opt.sim.zone1_distance = 0.9;  % zone 1 default distance (90% of line impedance)
% Don't forget to change this value (opt.sim.time_delay_ini) in solve_dae.m

opt.sim.ofgs_delta = 0.25;      % the fraction of generation that is shed during simulation if over frequency
opt.sim.ofgs_limit = 1.010;    % 1.02 threshold at which over frequency generation shedding occurs
opt.sim.ofgs_tdelay_ini = 0.2;  % 1 sec delay for uvls relay.

% initialize the case
ps = newpf(ps,opt);
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
ps = update_load_freq_source(ps);
% build the machine variables
[ps.mac,ps.exc,ps.gov] 		= get_mac_state(ps,'salient');
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
t_delay([ix.re.uvls]) = opt.sim.uvls_tdelay_ini;
t_delay([ix.re.ufls]) = opt.sim.ufls_tdelay_ini;
% t_delay([ix.re.dist]) = opt.sim.dist_tdelay_ini;
% t_delay([ix.re.temp]) = opt.sim.temp_tdelay_ini;
% t_delay([ix.re.ofgs])= opt.sim.ofgs_tdelay_ini;
% t_delay([ix.re.oc])= opt.sim.oc_tdelay_ini;

t_prev_check = nan(size(ps.relay,1),1);
dist2threshold = inf(size(ix.re.oc,2)*2,1);
state_a = zeros(size(ix.re.oc,2)*2,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   build an event matrix      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%%


load_attack = [35:51];
size_Att1 = size(load_attack,2);

freq1 = 9;
freq2 = 9;
freq3 = 9;
frac_A = NaN(3,1);
frac_A(1,1)=0.1;
frac_A(2,1)=0.1;
frac_A(3,1)=0.31;
[event1,k,step,event_1,event_2,event_3] = event_different_frequency(C,load_attack,t_max,frac_A,freq1,freq2,freq3);
area2_loads = find(event1(:,6)>=1 & event1(:,6)<35);
event1(area2_loads,:)=[];

%% build an event matrix
load_attack = [35:51];
size_Att2 = size(load_attack,2);
% time_step=[1,20,10,5,(1/0.3),(1/0.6),(1/0.4),(1/0.5),(1/0.8),(1/0.7),(1/0.9)];
% freq1 = rec_events(itr).freq(1);
% freq2 = rec_events(itr).freq(2);
% freq3 = rec_events(itr).freq(3);
freq1 = 9;
freq2 = 9;
freq3 = 9;
frac_A(1,1)=0.1;
frac_A(2,1)=0.1;
frac_A(3,1)=0.31;
[event,k,step,event_1,event_2,event_3] = event_different_frequency(C,load_attack,t_max,frac_A,freq1,freq2,freq3);
event(find(event(:,8) == frac_A(2,1) | event(:,8) == -frac_A(2,1) ),:)=[];

% event(find(event(:,8) == frac_A(1,1) | event(:,8) == -frac_A(1,1) ),:)=[];
% event = [event(1:69,:); event1(138:end-1,:); event1(end,:)];
% event(end,1) = 40; 

% to oscillate one load 
% oscillating_load=event(find(event(:,6)==49),:); % determine the load
% event = [event(1,:); oscillating_load; event(end,:)];% integrate it to event matrix


% event(512:end-1,:) = [];
% event(end,1) = event(end-1,1)+1;
%% Dr.kim Idea
%%%%%%%%%%%%%%%% Class A extra range %%%%%%%%%%%%

%% deleting area1 & 2
% area2_loads = find(event(:,6)>=1 & event(:,6)<35);
% event(area2_loads,:)=[];

%% deleting area 1 &3
% area2_loads = find(event(:,6)>=1 & event(:,6)>34 );
% event(area2_loads,:)=[];
% 
% area2_loads = find(event(:,6)>=1 & event(:,6)<18 );
% event(area2_loads,:)=[];
%% deleting area 2 & 3
% area2_loads = find(event(:,6)>=1 & event(:,6)>17 );
% event(area2_loads,:)=[];

% area3_loads = find(event(:,6)>=17 & event(:,6)<34);
% area2_loads = find(event(:,6)<28  & event(:,6)~=0 );
% area3_loads = find(event(:,6)<=15);
% event(area3_loads,8) = event(area3_loads,8).*-1;

% area2_loads = find(event(:,6)>=1 & event(:,6)<18 );

% area3_loads = find(event(:,6)==34);
% event(area3_loads,6) = 51;
% event(area3_loads,8) = event(area3_loads,8).*-1;
% event(area3_loads,8)=event(area3_loads,8).*(0.25/0.3);
% area3_loads = find(event(:,6)>37  & event(:,6)<40  & event(:,6)~=0 );
% event(area3_loads,8)=event(area3_loads,8).*(0.20/0.3);




% event(13:end-1,:)=[];
% event([167:end-1],:)=[];
% area3_loads = find(event(:,6)>=1 & event(:,6)<17);
% event(area3_loads,:)=[];
% event([14 15 16 17 18 19],:)=[];
% event1(53:end-1,1) = event1(53:end-1,1)-5;
% event = [event(1:154,:); event1(53:end-1,:); event1(end,:)];
% event = [event(1:154,:); event1(53:69,:); event1(end,:)];
% event(138:end-1,:)=[];
% event([155:end-1],:)=[];

% Exp-B
% event1(19:35,1) = event1(19:35,1)+10;

% event = [event(1:154,:); event1(19:35,:); event1(end,:)];
% event = [event(1:120,:); event1(19:35,:); event1(end,:)];
% event(69:end-1,:)=[];
% Exp-C

% event = [event(1:137,:); event1(104:end-1,:); event1(end,:)];

% Carter work
% event1(2:end,1) = event1(2:end,1)+60;
% 
% event = [event(1:end-1,:); event1(2:86,:); event1(end,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% New way for attack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % event(:,8) = 0;
% E = find(event(:,6)==36);
% E =E(E>138);
% E = event(E,:);
% E(:,8) = (E(:,8)./0.25).*0.05;
% event(138:end-1,:) = [];
% event = [event(1:end-1,:);E;event(end,:)];
% % event(end,1) = 30;

attack_event = find(event(:,6)<=17 & event(:,8)>0 & event(:,1)>=10);
ps.time_area1 = unique(event(attack_event,1));

attack_event = find(event(:,6)>17 & event(:,6)<35 & event(:,8)>0 & event(:,1)>=10);
ps.time_area2 = unique(event(attack_event,1));

attack_event = find(event(:,6)>=35 & event(:,8)>0 & event(:,1)>=10);
ps.time_area3 = unique(event(attack_event,1));


attack_event = find(event(:,6)<=17 & event(:,8)<0 & event(:,1)>=10);
ps.time_area1_back = unique(event(attack_event,1));

attack_event = find(event(:,6)>17 & event(:,6)<35 & event(:,8)<0 & event(:,1)>=10);
ps.time_area2_back = unique(event(attack_event,1));

attack_event = find(event(:,6)>=35 & event(:,8)<0 & event(:,1)>=10);
ps.time_area3_back = unique(event(attack_event,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Old way for attack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % delete_load = find(event(:,6)>17 & event(:,6)<35 );
% attack_event = find(event(:,6)==37 & event(:,8)>0 & event(:,1)>10);
% event(attack_event,2:10) = 0;
% event(attack_event,C.ev.type) = C.ev.em_control;
% event(attack_event,C.ev.shunt_loc) =  NaN;
% event(attack_event,C.ev.change_by) =  1;
% 
% attack_event = find(event(:,6)==37 & event(:,8)<0 & event(:,1)>10);
% event(attack_event,2:10) = 0;
% event(attack_event,C.ev.type) = 23;
% event(attack_event,C.ev.shunt_loc) =  NaN;
% event(attack_event,C.ev.change_by) =  1;
% 
% delete_load = find(event(:,6)>34 & event(:,1)>=10);
% event(delete_load,:) = [];
% 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t=0;
% counter =1;
% while t == 0
%     if rec(itr).outputs.event_record(counter,2) == 11
%         t = rec(itr).outputs.event_record(counter,1) ;
%     end
%     counter = counter + 1;   
%     
% end
% event(2,[C.ev.time C.ev.type]) = [3 C.ev.trip_branch];
% event(2,C.ev.branch_loc) = 103;
% 
% event(3,[C.ev.time C.ev.type]) = [3 C.ev.trip_branch];
% event(3,C.ev.branch_loc) = 98;
% 
% event(4,[C.ev.time C.ev.type]) = [3 C.ev.trip_branch];
% event(4,C.ev.branch_loc) = 97;
% 
% event(5,[C.ev.time C.ev.type]) = [3 C.ev.trip_branch];
% event(5,C.ev.branch_loc) = 99;

% event(6,[C.ev.time C.ev.type]) = [3 C.ev.trip_branch];
% event(6,C.ev.branch_loc) = 97;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%   adding  EC       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EC_time = [24:20:40];
% EC_time = [18.5];
% for i = 1:size(EC_time,2)
%     EC_event =zeros(1,10);
%     EC_event_less_t  = find(event(:,C.ev.time)<EC_time(1,i));
%     EC_event_more_t = find(event(:,C.ev.time)>=EC_time(1,i));
%     EC_event(1,C.ev.time)      =  EC_time(1,i);
%     EC_event(1,C.ev.type)      =  C.ev.em_control;
%     EC_event(1,C.ev.shunt_loc) =  NaN;
% %     EC_event(1,C.ev.quantity)  =  0;
%     EC_event(1,C.ev.change_by) =  1;
%     event = [event(EC_event_less_t,:);EC_event;event(EC_event_more_t,:)];
% 
% end

% area2_loads = find(event(:,6)~=0 & event(:,2)~=22 & event(:,2)~=1);
% event(area2_loads,:)=[];

% 2;
% 
% load_attack_percent = sum(ps.shunt(load_attack,2))*0.73/sum(ps.shunt(:,2));

% event_1([30:33],8)=event_1([30:33],8).*-1;
% event_2([29 33:34],8)=event_2([29 33:334],8).*-1;
% event = [event_3([10:18],:);event_1([30:33],:);event_2([29 33:34],:);event_3([27:35],:);];

    ps.original  =  ps.shunt(:,C.sh.P);
    ps.target = [119];
    ps.itr = 0;
    ps.f_int = 103.1693;
    ps.count1 = 0;
    ps.count2 = 0;
%%
% Z =NaN(73,7);
% for i = 1:73
%     Z(i,1) = ps.bus(i,1);
%     Z1 = find(ps.branch(:,1) == ps.bus(i,1));
%     Z1 = ps.branch(Z1,2);
%     Z2 = find(ps.branch(:,2) == ps.bus(i,1));
%     Z2 = ps.branch(Z2,1);
%     Z3 = [];
%     Z3 = [Z1;Z2];
%     for j = 1:size(Z3,1)
%         Z3(j) = find(ps.bus(:,1)== Z3(j));
%     end
%     Z3 = sort(unique(Z3));
%     Z(i,2:size(Z3,1)+1) = Z3';
%     
% end

ps.goal = 1;
%% run the simulation
blockout = 0;
k = 1;
% % [outputs,ps,blockout,event] = simgrid_attack_toplogy_RTS(ps,event,'sim_caseRTS_96_EC',opt,k,blockout,step,event_1,event_2,event_3,load_attack);
% [outputs,ps,blockout,event] = simgrid_dynamic_attack_RTS(ps,event,'sim_caseRTS_96_EC',opt,blockout,event_1,event_2,event_3,load_attack);
[outputs,ps,blockout] = simgrid_EC(ps,event,'sim_caseRTS_96_EC',opt,k,blockout,step,event_1,event_2,event_3,load_attack);
% [outputs,ps,blockout,event,t_simulation,adaptive_summary] = simgrid_adaptive_control(ps,event,'sim_caseRTS_96_EC',opt,k,blockout,event_1,event_2,event_3,load_attack);



total_losses(itr,2) = outputs.demand_lost;
total_losses(itr,3) = sum(ps.branch(:,C.br.status));
% loads(itr).after_distur = ps.shunt(:,[1 2]);
% end
%% print the results
fname = outputs.outfilename;
% [t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd,P3,Temperature,If,It] = read_outfile(fname,ps,opt);
[t,delta,omega,Pm,Eap,Vmag,theta,E1,Efd,P3,Temperature,If,It,Pf,Qf,P_inj] = read_outfile_attack(fname,ps,opt,outputs);
omega_0 = 2*pi*ps.frequency;
omega_pu = omega / omega_0;

%%
