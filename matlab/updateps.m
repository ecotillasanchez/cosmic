function ps = updateps(ps)
% check for and fix irregular data in a power system data structure
% see "psconstants" for a description of this structure

C = psconstants;

%% make sure that the matrices are big enough
% bus
[n,bu_cols] = size(ps.bus);
if bu_cols < C.bu.cols
    ps.bus = addcolumns(ps.bus,C.bu.cols);
end
% branch
[m,br_cols] = size(ps.branch);
if br_cols < C.br.cols
    ps.branch = addcolumns(ps.branch,C.br.cols);
end
% gen
if size(ps.gen,2) < C.ge.cols
    ps.gen = addcolumns(ps.gen,C.ge.cols);
end
% shunt
if ~isfield(ps,'shunt')
    ps.shunt = [];
end
if size(ps.shunt,2) < C.sh.cols
    ps.shunt = addcolumns(ps.shunt,C.sh.cols);
end

%% make sure that the system has a base frequency
if ~isfield(ps,'frequency')
    ps.frequency = 60;
end

%% check the bus data and update bus_i
max_bus_no  = max(ps.bus(:,1));
ps.bus_i    = sparse(ps.bus(:,1),1,(1:n)',max_bus_no,1);
no_baseKV   = ps.bus(:,C.bu.baseKV)==0;
if any(no_baseKV)
    ps.bus(no_baseKV,C.bu.baseKV) = C.bu.baseKV_default;
end

%% set generator voltages
gen_status = ps.gen(:,C.ge.status);
is_v_fixed = (ps.gen(:,C.ge.type)==C.PV | ps.gen(:,C.ge.type)==C.REF) & gen_status;
% put the generator voltage set points into the bus matrix
ps.bus(ps.bus_i(ps.gen(is_v_fixed,1)),C.bu.Vmag) = ps.gen(is_v_fixed,C.ge.Vsp);
% shut off generators that are off
ps.gen(~gen_status,C.ge.P:C.ge.Q) = 0;

%% check the generator types
gen_set = (ps.gen(:,C.ge.type)==0);
if any( gen_set )
    bus_ix = ps.bus_i(ps.gen(gen_set,1));
    ps.gen(gen_set,C.ge.type) = ps.bus(bus_ix,C.bu.type);
end
% check that generators have participation factors
ramp_set = (ps.gen(:,C.ge.type)==C.REF | ps.gen(:,C.ge.type)==C.PF) & ps.gen(:,C.ge.P)>10;
fix_set = ramp_set &  ps.gen(:,C.ge.part_fact)==0;
ps.gen(fix_set,C.ge.part_fact) = 1;
% add generator id numbers
needs_id = ps.gen(:,C.ge.id) <= 0;
if any( needs_id )
    max_id = max(max(ps.gen(:,C.ge.id)),0);
    ps.gen(needs_id,C.ge.id) = max_id + (1:sum(needs_id))';
end
ng = size(ps.gen,1);
ps.gen_i    = sparse(ps.gen(:,1),1,(1:ng)',max_bus_no,1);

%% check branch data
% fix status
status = ps.branch(:,C.br.status);
if any( status~=0 & status~=1 )
    warning('PowerSystems:BranchStatus','Strange branch status variables found');
end
ps.branch(status>1,C.br.status)=1;
ps.branch(status<1,C.br.status)=0;

% fix any transformer taps listed as zero
tap = ps.branch(:,C.br.tap);
ps.branch(tap==0,C.br.tap) = 1;

% populate branch ids if needed and update branch_i
needs_id = ps.branch(:,C.br.id) <= 0;
if any( needs_id )
    max_id = max(max(ps.branch(:,C.br.id)),0);
    ps.branch(needs_id,C.br.id) = max_id + (1:sum(needs_id))';
end
max_br_no   = max(ps.branch(:,C.br.id));
ps.branch_i = sparse(ps.branch(:,C.br.id),1,(1:m)',max_br_no,1);

%% check shunt data
% check shunt status
if isfield(ps,'shunt') && ~isempty(ps.shunt)
    status = ps.shunt(:,C.sh.status);
    if any( status<0 | status>1 )
        error('strange status/factor value found in shunt matrix');
    end
end

% add in shunt id#s
needs_id = ps.shunt(:,C.sh.id) <= 0;
if any( needs_id )
    max_id = max(max(ps.shunt(:,C.sh.id)),0);
    ps.shunt(needs_id,C.sh.id) = max_id + (1:sum(needs_id))';
end
ps.shunt_i = sparse(ps.shunt(:,1),1,(1:size(ps.shunt,1))',max_bus_no,1);

%% addcolumns function
function M = addcolumns(M,ncols)

[n,m] = size(M);

M = [M zeros(n,ncols-m)];
