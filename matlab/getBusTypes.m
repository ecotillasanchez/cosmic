function [is_pv,is_pq,is_ref,is_pf,PartFact] = getBusTypes(ps,bus_subset)
% Use this function to figure out which buses are types PV, PQ, or
% Ref/Swing/Slack
% usage: [pv,pq,ref,pf] = getBusTypes(ps,bus_subset)
%  ps is the power system data structure
%  bus_subset (optional) is a vector indicating which nodes are in the
%  current subgrid requested. If specified only these nodes will be dealt
%  with.

C = psconstants;
%% check the inputs
n = size(ps.bus,1);
gen_type = ps.gen(:,C.ge.type);
is_gen_on = ps.gen(:,C.ge.status);
if nargin<2 || isempty(bus_subset)
    bus_subset = true(n,1);
end

%% calculate generation at the buses
Pg = ps.gen(:,C.ge.P) .* is_gen_on;
Pg_bus = sparse(ps.bus_i(ps.gen(:,1)),1,Pg,n,1);

%% find the reference bus for the nodes in bus_subset
% first check for an infinite bus
ref_bus_set = find(ps.bus(:,C.bu.type)==C.INF);
if isempty(ref_bus_set)
    is_ref_gen = (gen_type==C.REF & is_gen_on);
    ref_bus_set = setdiff(ps.bus_i( ps.gen(is_ref_gen,1) ),find(~bus_subset));
end

n_ref = length(ref_bus_set);
if n_ref==0 % find a reference bus if we don't have one
    % choose the largest gen as the reference bus
    [p,bus_ix] = max( Pg_bus .* bus_subset );
    if p>1e-6
        ref = bus_ix;
    else
        ref = [];
    end
elseif n_ref==1
    ref = ref_bus_set;
elseif n_ref>1
    [~,ix] = max( Pg_bus(ref_bus_set) );
    ref = ref_bus_set(ix);
end

is_ref = false(n,1);
is_ref(ref) = true;

%% find the pv buses
pv = ps.bus_i( ps.gen(gen_type==C.PV & is_gen_on,1) );
is_pv  = false(n,1);
is_pv(pv) = true;
is_pv(~bus_subset) = false;
is_pv(is_ref) = false;

%% find the participation factor generators
pf = ps.bus_i( ps.gen(gen_type==C.PF & is_gen_on,1) );
is_pf = false(n,1);
is_pf(pf) = true;
is_pf(~bus_subset) = false;
is_pf(is_ref) = false;

%% find the pq buses
is_pq = bus_subset & ~is_ref & ~is_pv & ~is_pf;

%% calculate participation factors
if nargout>4
    PartFact = zeros(n,1);
    gen_bus_i = ps.bus_i( ps.gen(:,1) );
    PartFact(gen_bus_i) = ps.gen(:,C.ge.part_fact);
    if all(PartFact==0)
        PartFact(ref) = 1;
    end
end

%{
if all(is_pq(bus_subset))
    warning('All pq buses');
elseif isempty(ref)
    warning('did not find a reference bus');
end
%}

return
