function [Vmag,theta,success] = subgridPowerFlow(bus_subset,Ybus,Vmag,theta,Sg,Sd,pq,pv,ref,part_fact,nr_opts)
% usage: [Vmag,theta,success] = subgridPowerFlow(bus_subset,Ybus,Vmag,theta,Sg,Sd,pq,pv,ref,part_fact,nr_opts)
% power flow for a subgrid defined by bus_subset
n_sub = sum(bus_subset);
if sum(ref)==0
    Vmag(bus_subset)  = 0;
    theta(bus_subset) = 0;
    success = false;
elseif n_sub==1 % we have one generator bus
    theta(bus_subset) = 0;
    success = true;
else
    % build the decision vector
    npq = sum(pq & bus_subset);
    ix = struct;
    ix.theta = 1:(n_sub-1);
    ix.Vmag  = (1:npq) + max(ix.theta);
    ix.rho = 1 + max(ix.Vmag); % generator ramping variable
    if isempty(ix.rho) % degenerate case
        ix.rho = 1 + max(ix.theta);
    end
    nx = max(ix.rho);
    x = zeros(nx,1);
    x(ix.theta) = theta(~ref & bus_subset);
    x(ix.Vmag)  = Vmag(pq);
    x(ix.rho)   = 0;
    if any(x(ix.Vmag))==0
        x(ix.Vmag) = 1;
    end
    % subset the variables for this subgrid
    Ybus_sub    = Ybus(bus_subset,bus_subset);
    Vmag_sub    = Vmag(bus_subset);
    Sg_sub      = Sg(bus_subset);
    Sd_sub      = Sd(bus_subset,:);
    pq_sub      = pq(bus_subset);
    pv_sub      = pv(bus_subset);
    ref_sub     = ref(bus_subset);
    if part_fact(ref)==0
        part_fact(ref) = 1;
    end
    % set up a virtual function to solve
    g = @(newx)mismatch(newx,Ybus_sub,Vmag_sub,Sg_sub,Sd_sub,pq_sub,pv_sub,ref_sub,part_fact);
    % try to solve the power flow problem
    [x,flag] = nrsolve(g,x,nr_opts);
    
    % if this failed re-try with a flat start
    if flag~=1
        if nr_opts.nr.verbose
           disp('Re-trying power-flow with a flat start');
        end
        x(ix.theta) = 0;
        x(ix.Vmag) = 1;
        [x,flag] = nrsolve(g,x,nr_opts);
    end
    if flag==-1
        disp('Found a power flow solution that is optimal, but not feasible.');
    end
    success = (flag==1);
    % save the results
    Vmag(pq) = x(ix.Vmag);
    theta(~ref & bus_subset) = x(ix.theta);
    theta(ref) = 0;
end