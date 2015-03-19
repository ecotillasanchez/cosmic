function [ps,t_out,X,Y] = superset_odeout(ps,ps1,ps2,t_out1,t_out2,X1,Y1,X2,Y2,opt)
% usage: [ps,t_out,X,Y] = superset_odeout(ps,ps1,ps2,t_out1,t_out2,X1,Y1,X2,Y2,opt)

C = psconstants;
angle_ref = opt.sim.angle_ref;                 % angle reference: 0:delta_sys,1:delta_coi

% get number of buses, machines and branches
n   = size(ps.bus,1);   n_macs  = size(ps.mac,1);   m   = size(ps.branch,1);
n1  = size(ps1.bus,1);  n_macs1 = size(ps1.mac,1);  m1  = size(ps1.branch,1);
n2  = size(ps2.bus,1);  n_macs2 = size(ps2.mac,1);  m2  = size(ps2.branch,1);

nxsmac = 7;				% number of differential variables per machine 

% store upstream the ps information from each island
buses1      = ismember(ps.bus(:,1),ps1.bus(:,1));       buses2  = ~buses1;
gens1       = ismember(ps.gen(:,1),ps1.gen(:,1));       gens2   = ~gens1;
shunts1     = ismember(ps.shunt(:,1),ps1.shunt(:,1));   shunts2 = ~shunts1;
relx_bu1    = ismember(ps.relay(:,C.re.bus_loc),ps1.relay(ps1.relay(:,C.re.bus_loc)~=0,C.re.bus_loc));
relx_br1    = ismember(ps.relay(:,C.re.branch_loc),ps1.relay(ps1.relay(:,C.re.branch_loc)~=0,C.re.branch_loc));
relx_ge1    = ismember(ps.relay(:,C.re.gen_loc),ps1.relay(ps1.relay(:,C.re.gen_loc)~=0,C.re.gen_loc));
relx_sh1    = ismember(ps.relay(:,C.re.shunt_loc),ps1.relay(ps1.relay(:,C.re.shunt_loc)~=0,C.re.shunt_loc));
relx1       = relx_bu1 | relx_br1 | relx_ge1 | relx_sh1; 
relx_bu2    = ismember(ps.relay(:,C.re.bus_loc),ps2.relay(ps2.relay(:,C.re.bus_loc)~=0,C.re.bus_loc));
relx_br2    = ismember(ps.relay(:,C.re.branch_loc),ps2.relay(ps2.relay(:,C.re.branch_loc)~=0,C.re.branch_loc));
relx_ge2    = ismember(ps.relay(:,C.re.gen_loc),ps2.relay(ps2.relay(:,C.re.gen_loc)~=0,C.re.gen_loc));
relx_sh2    = ismember(ps.relay(:,C.re.shunt_loc),ps2.relay(ps2.relay(:,C.re.shunt_loc)~=0,C.re.shunt_loc));
relx2       = relx_bu2 | relx_br2 | relx_ge2 | relx_sh2; 

if ~angle_ref
ps.bus(buses1,C.bu.delta_sys)       = ps1.bus(:,C.bu.delta_sys);
ps.bus(buses2,C.bu.delta_sys)       = ps2.bus(:,C.bu.delta_sys);
end
ps.bus(buses1,C.bu.Vmag:C.bu.Vang)  = ps1.bus(:,C.bu.Vmag:C.bu.Vang);
ps.bus(buses2,C.bu.Vmag:C.bu.Vang)  = ps2.bus(:,C.bu.Vmag:C.bu.Vang);
ps.mac(gens1,C.ma.Eap:C.ma.omega)   = ps1.mac(:,C.ma.Eap:C.ma.omega);
ps.mac(gens2,C.ma.Eap:C.ma.omega)   = ps2.mac(:,C.ma.Eap:C.ma.omega);
ps.exc(gens1,C.ex.Efd:C.ex.E1)      = ps1.exc(:,C.ex.Efd:C.ex.E1);
ps.exc(gens2,C.ex.Efd:C.ex.E1)      = ps2.exc(:,C.ex.Efd:C.ex.E1);
ps.relay(relx1,C.re.setting1:C.re.timer_start) = ps1.relay(:,C.re.setting1:C.re.timer_start);
ps.relay(relx2,C.re.setting1:C.re.timer_start) = ps2.relay(:,C.re.setting1:C.re.timer_start);

% synchronize time vectors from the islands by uniform resampling at PMU time
if all(isnan(X1(:)))
    t_out1  = [t_out1 t_out2(end)];
    X1      = [X1 nan(size(X1,1),1)];
    Y1      = [Y1 nan(size(Y1,1),1)];
elseif all(isnan(X2(:)))
    t_out2  = [t_out2 t_out1(end)];
    X2      = [X2 nan(size(X2,1),1)];
    Y2      = [Y2 nan(size(Y2,1),1)];
end

% to avoid losing data when two time range are way differernt, during
% synchronizing time series.
if t_out1(end)>t_out2(end)
    nt2 = length(t_out2);
    t_out2 = [t_out2(1:end-1),t_out2(end):opt.sim.dt_default:t_out1(end)];
    nt2_new = length(t_out2);
    X2 = [X2,nan(size(X2,1),nt2_new-nt2)];
    Y2 = [Y2,nan(size(Y2,1),nt2_new-nt2)];
    ps2.shunt(:,C.sh.factor) = 0;      % where there is a numerical issue like this, declare it blackout.
elseif t_out2(end)>t_out1(end)
    nt1 = length(t_out1);
    t_out1 = [t_out1(1:end-1),t_out1(end):opt.sim.dt_default:t_out2(end)];
    nt1_new = length(t_out1);
    X1 = [X1,nan(size(X1,1),nt1_new-nt1)];
    Y1 = [Y1,nan(size(Y1,1),nt1_new-nt1)];
    ps1.shunt(:,C.sh.factor) = 0;      % where there is a numerical issue like this, declare it blackout.
end

ps.shunt(shunts1,C.sh.factor)       = ps1.shunt(:,C.sh.factor);
ps.shunt(shunts2,C.sh.factor)       = ps2.shunt(:,C.sh.factor);

ts1x        = timeseries(X1',t_out1); ts1y  = timeseries(Y1',t_out1);
ts2x        = timeseries(X2',t_out2); ts2y  = timeseries(Y2',t_out2);
[ts1x, ts2x] = synchronize(ts1x,ts2x,'Uniform','Interval',opt.sim.dt_default);
[ts1y, ts2y] = synchronize(ts1y,ts2y,'Uniform','Interval',opt.sim.dt_default);
X1          = ts1x.Data';   Y1  = ts1y.Data';    
X2          = ts2x.Data';   Y2  = ts2y.Data';
t_out       = ts1x.time';
tk          = length(t_out);

% merge the machine portion of X vectors
X_macs  = zeros(nxsmac*n_macs,tk);
if n_macs1>0
    for i = 1:n_macs1
        j = find(ps.mac(:,1) == ps1.mac(i,1));
        ind     = nxsmac*j-(nxsmac-1) : 1 : nxsmac*j ;
        ind1    = nxsmac*i-(nxsmac-1) : 1 : nxsmac*i ;
        X_macs(ind,:) = X1(ind1,:);
    end
end

if n_macs2>0
    for i = 1:n_macs2,
        j = find(ps.mac(:,1) == ps2.mac(i,1));
        ind     = nxsmac*j-(nxsmac-1) : 1 : nxsmac*j ;
        ind2    = nxsmac*i-(nxsmac-1) : 1 : nxsmac*i ;
        X_macs(ind,:) = X2(ind2,:);
    end
end

% merge the temperature portion of X vectors
X_br    = zeros(m,tk);
pars1   = false;
pars2   = false;

if m1>0
    for i = 1:m1
        j = find(ismember(ps.branch(:,1:2),ps1.branch(i,1:2),'rows'));
        if length(j) > 1
            if length(j) > 2, error('three parallel branches?'); end
            % we have a pair of parallel branches
            if pars1
                % we are on the second one
                j = j(2);
                pars1 = false;
            else
                % we are on the first one
                j = j(1);
                pars1 = true;
            end
        end
        X_br(j,:)   = X1(nxsmac*n_macs1 + i,:);
    end
end

if m2>0
    for i = 1:m2
        j = find(ismember(ps.branch(:,1:2),ps2.branch(i,1:2),'rows'));
        if length(j) > 1
            if length(j) > 2, error('three parallel branches?'); end
            % we have a pair of parallel branches
            if pars2
                % we are on the second one
                j = j(2);
                pars2 = false;
            else
                % we are on the first one
                j = j(1);
                pars2 = true;
            end
        end
        X_br(j,:)   = X2(nxsmac*n_macs2 + i,:);
    end
end
X = [X_macs; X_br];
% to avoid issues of synchronzing data the in the next level, since X_macs
% might be empty
if isempty(X)
    X = [X;nan(1,size(X,2))];
end

% merge the algebraic variables
if ~angle_ref
    Y   = zeros(2*n + 1,tk);
else
    Y   = zeros(2*n,tk);
end
for i = 1:n1
    j = find(ps.bus(:,1) == ps1.bus(i,1));
    Y((2*j)-1:2*j,:) = Y1((2*i)-1:2*i,:);
end
for i = 1:n2
    j = find(ps.bus(:,1) == ps2.bus(i,1));
    Y((2*j)-1:2*j,:) = Y2((2*i)-1:2*i,:);
end

% merge event record
ps.event_record = unique([ps.event_record; ps1.event_record; ps2.event_record],'rows');