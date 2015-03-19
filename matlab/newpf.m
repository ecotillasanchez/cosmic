function [ps,converged,islands] = newpf(ps,opt)
% experimental power flow formulation
% that allows for multiple slack buses (participation factors) and
% uses nrsolve to obtain the solution
% usage: [ps,status,islands] = newpf(ps)
% status is a binary status (1: success, 0: failure)

% TODO
% -work on cascading power flow

%% constants
j = 1i;

C = psconstants;

%% process the inputs
if nargin<2
    opt = psoptions;
end
nr_opts = opt;
nr_opts.nr.verbose = opt.verbose;
nr_opts.nr.linesearch = opt.pf.linesearch;
% update
if opt.pf.update
    ps = updateps(ps);
end

%% initialize outputs
converged  = true;
islands = struct;

%% extract some basic data
% bus data
n = size(ps.bus,1);
[Ybus,Yf,Yt,~,~] = getYbus(ps,false);
V0 = getV(ps);
Vmag0  = abs(V0);
theta0 = angle(V0);
Vmag = Vmag0;
theta = theta0;
% get the gen bus injections
nGen    = size(ps.gen,1);
G       = ps.bus_i(ps.gen(:,1));
status  = ps.gen(:,C.ge.status);
S_gen   = status.*(ps.gen(:,C.ge.P) + j*ps.gen(:,C.ge.Q))/ps.baseMVA;
Sg      = sparse(G,1,S_gen,n,1);

% get the load bus injections with a ZIPE model
D               = ps.bus_i(ps.shunt(:,1));
S_load_base     = (ps.shunt(:,C.sh.P) + j*ps.shunt(:,C.sh.Q)).*ps.shunt(:,C.sh.factor)/ps.baseMVA;
S_load_P        = S_load_base.*ps.shunt(:,C.sh.frac_S);
Sd              = sparse(D,3,S_load_P,n,5);
S_load_Z        = S_load_base.*ps.shunt(:,C.sh.frac_Z);
Sd              = Sd + sparse(D,1,S_load_Z,n,5);
S_load_I        = S_load_base.*(1-(ps.shunt(:,C.sh.frac_Z)+ps.shunt(:,C.sh.frac_S)+ps.shunt(:,C.sh.frac_E)));
Sd              = Sd + sparse(D,2,S_load_I,n,5);
S_load_E        = S_load_base.*ps.shunt(:,C.sh.frac_E);
Sd              = Sd + sparse(D,4,S_load_E,n,5);
S_load_E_gamma  = ps.shunt(:,C.sh.gamma);
Sd              = Sd + sparse(D,5,S_load_E_gamma,n,5);

% check for islands
br_status = ps.branch(:,C.br.status)==1;
if opt.pf.CalcIslands
    [subGridNos,nSubGrids] = findSubGraphs(ps.bus(:,1),ps.branch(br_status,1:2));
else
    nSubGrids = 1;
    subGridNos = true(n,1);
end

%% run the power flow for each island/sub-grid
% the output of this is to calculate Vmag and theta for each bus
for grid_no = 1:nSubGrids
    % find the buses/types in this subgrid
    bus_subset = (grid_no==subGridNos);
    n_sub = sum(bus_subset);
    [pv,pq,ref,pf,part_fact] = getBusTypes(ps,bus_subset);
    if opt.pf.CascadingPowerFlow
        % make all of the generators pf buses
        pf(pv) = true;      %#ok<NASGU>
        pv = false(n,1);
    end
    if all(abs(Vmag(bus_subset))<1e-12) || all(pq(bus_subset)) || ~any(ref(bus_subset))
        Vmag(bus_subset) = 0;
        theta(bus_subset) = 0;
        if opt.verbose
            fprintf('Subgrid %d of %d (%d buses) is black\n',grid_no,nSubGrids,n_sub);
        end
        continue;
    end
    % print something
    if opt.verbose && nSubGrids>1
        fprintf('Computing power flow for subgrid %d of %d (%d buses)\n',grid_no,nSubGrids,n_sub);
    end
    % run power flow on the subgrid
    [Vmag,theta,success] = subgridPowerFlow(bus_subset,Ybus,Vmag,theta,Sg,Sd,pq,pv,ref,part_fact,nr_opts);
    % cascading power flow
    while ~success && opt.pf.CascadingPowerFlow && sum(ref)>0
        error('this part of the code needs revising in order to incorporate the zipe load model');
        % the odds are that we arrived here because of voltage collapse
        % therefore find the shunts at low-voltage buses
        ls_rate     = opt.pf.load_shed_rate; %#ok<UNRCH>
        low_voltage = Vmag<ps.bus(:,C.bu.Vmin)&bus_subset;
        if ~success && ~any(low_voltage)
            fprintf('Blackout in subgrid %d\n',grid_no);
            success = false;
            break
        end
        low_voltage_shunts = ismember(ps.shunt(:,1),ps.bus(low_voltage,1));
        % reduce load by load shed rate
        change = false;
        for sh=find(low_voltage_shunts)'
            factor_old = ps.shunt(sh,C.sh.factor);
            if factor_old>0
                change = true;
                factor_new = max(0,factor_old-ls_rate);
                ps.shunt(sh,C.sh.factor) = factor_new;
                Smag = abs(ps.shunt(sh,C.sh.P));
                % print something
                bu_no = ps.shunt(sh,1);
                fprintf('Reducing load at bus %4d by %.2f%% (from %7.2f to %7.2f MVA)\n',...
                    bu_no,ls_rate*100,factor_old*Smag,factor_new*Smag);
            end
        end
        if ~change
            fprintf('Blackout in subgrid %d of %d (%d buses)\n',grid_no,nSubGrids,n_sub);
            success = false;
            keyboard
            break;
        end
        % get the new Sd
        factor = ps.shunt(:,C.sh.frac_S).*ps.shunt(:,C.sh.factor);
        S_load_P = factor.*S_load_base;
        Sd = sparse(D,1,S_load_P,n,1);
        % get the revised Yshunt
        factor = ps.shunt(:,C.sh.frac_Y).*ps.shunt(:,C.sh.factor);
        y_shunt = factor.*S_load_base;
        Yshunt = sparse(D,D,y_shunt,n,n);
        % check to see if all of the loads in this area are shut down
        shunt_subset = ismember(ps.shunt(:,1),ps.bus(bus_subset,1));
        if all(ps.shunt(shunt_subset,C.sh.factor)<=0)
            fprintf('Blackout in subgrid %d of %d (%d buses)\n',grid_no,nSubGrids,n_sub);
            success = false;
            keyboard
            break;
        end
        % re-run power flow on the subgrid
        Vmag(bus_subset)  = Vmag0(bus_subset);
        theta(bus_subset) = theta0(bus_subset);
        [Vmag,theta,success] = ...
            subgridPowerFlow(bus_subset,Ybus+Yshunt,Vmag,theta,Sg-Sd,pq,pv,ref,part_fact,nr_opts);
        if ~success
            warning('Debug code here');
            figure(3); hold on; %DEBUG
            plot(sort(Vmag));   %DEBUG
        end
    end
    % set the subgrid to blackout state
    if ~success
        Vmag(bus_subset) = 0;
        theta(bus_subset) = 0;
    end
    % save info about this island
    islands(grid_no).status = success;
    islands(grid_no).bus_set = ps.bus(bus_subset,1);
    % set the status variable
    converged = converged & success;
end

%% save the results in the ps structure
Vmag(Vmag<eps) = 0;
V = Vmag.*exp(j*theta);
% bus voltage results
ps.bus(:,C.bu.Vmag) = Vmag;
ps.bus(:,C.bu.Vang) = theta*180/pi;

% branch results
If = Yf*V; % branch status is accounted for in Yf
It = Yt*V; % branch status is accounted for in Yt
F = ps.bus_i(ps.branch(:,1));
T = ps.bus_i(ps.branch(:,2));
Sf = V(F) .* conj(If);
St = V(T) .* conj(It);
ps.branch(:,C.br.Imag_f) = abs(If);
ps.branch(:,C.br.Imag_t) = abs(It);
ps.branch(:,C.br.Pf) = real(Sf) * ps.baseMVA;
ps.branch(:,C.br.Qf) = imag(Sf) * ps.baseMVA;
ps.branch(:,C.br.Pt) = real(St) * ps.baseMVA;
ps.branch(:,C.br.Qt) = imag(St) * ps.baseMVA;

% calculate resulting ZIPE load after the powerflow
zipe_cols       = size(Sd,2);
if zipe_cols == 1
    S_zipe = Sd;
elseif zipe_cols == 5
    S_Z = Sd(:,1) .* Vmag.^2;
    S_I = Sd(:,2) .* Vmag;
    S_P = Sd(:,3);
    S_E = Sd(:,4) .* Vmag.^Sd(:,5);
    S_zipe = S_Z + S_I + S_P + S_E;
else
    error('zipe load model matrix is not the right size');
end

% calculate generator outputs after the powerflow
Sg = V.*conj(Ybus*V) + S_zipe;
Pg = real(Sg);
Qg = imag(Sg);

%  assign to gen structure
for gi = 1:nGen
    gen_bus = ps.gen(gi,1);
    gen_bus_i = ps.bus_i(gen_bus);
    if Vmag(gen_bus_i) == 0
        ps.gen(gi,C.ge.Qg) = 0;
        ps.gen(gi,C.ge.Pg) = 0;
        if opt.pf.CascadingPowerFlow
            ps.gen(gi,C.ge.status) = 0;
        end
    else
        gens_at_bus = (ps.gen(:,1)==gen_bus);
        if sum(gens_at_bus)==1
            Pg_ratio = 1;
        else
            Pg_ratio = ps.gen(gi,C.ge.Pmax)/sum(ps.gen(gens_at_bus,C.ge.Pmax));
        end
        ps.gen(gi,C.ge.Qg) = Qg(gen_bus_i) * ps.baseMVA * Pg_ratio;
        ps.gen(gi,C.ge.Pg) = Pg(gen_bus_i) * ps.baseMVA * Pg_ratio;
    end
end

% shut off failed shunts
if opt.pf.CascadingPowerFlow
    for sh = 1:size(ps.shunt,1)
        sh_bus_i = ps.bus_i(ps.shunt(sh,1));
        if Vmag(sh_bus_i) == 0
            ps.shunt(sh,C.sh.factor) = 0;
        end
    end
end

return
