function ps = subsetps(ps,bus_subset)
% usage: ps = subsetps(ps,bus_subset)
% subset a power system given the nodes specified in 
% "subset."  The input "subset" can be either a binary/logical vector
% with one entry per bus, or a list of bus numbers

C   = psconstants;

if islogical(bus_subset)
    bus_nos = ps.bus(bus_subset,1);
else
    bus_nos = bus_subset;
    bus_subset = ismember(ps.bus(:,1),bus_nos);    
end
% subset the buses
ps.bus = ps.bus(bus_subset,:);
n = size(ps.bus,1);

% figure out which branches to keep
F = ps.branch(:,1);
T = ps.branch(:,2);
br_keep = ismember(F,bus_nos) & ismember(T,bus_nos);
br_keep_nos = ps.branch(br_keep,C.br.id);
% subset the branches
ps.branch = ps.branch(br_keep,:);
m = size(ps.branch,1);

% subset the generators (including dynamics)
ge_keep = ismember(ps.gen(:,1),bus_nos);
ge_keep_nos = ps.gen(ge_keep,C.ge.id);
ps.gen    = ps.gen(ge_keep,:);
if isfield(ps,'mac'), ps.mac = ps.mac(ge_keep,:); end
if isfield(ps,'exc'), ps.exc = ps.exc(ge_keep,:); end
if isfield(ps,'gov'), ps.gov = ps.gov(ge_keep,:); end

% subset the shunts
if isfield(ps,'shunt')
    sh_keep = ismember(ps.shunt(:,1),bus_nos);
    sh_keep_nos = ps.shunt(sh_keep,C.sh.bus);
    ps.shunt = ps.shunt(sh_keep,:);
end

% subset branch relays
if isfield(ps,'relay')
    re_keep_br = ismember(ps.relay(:,C.re.branch_loc),br_keep_nos);    
    re_keep_bu = ismember(ps.relay(:,C.re.bus_loc),bus_nos);
    re_keep_ge = ismember(ps.relay(:,C.re.gen_loc),ge_keep_nos);
    re_keep_sh = ismember(ps.relay(:,C.re.shunt_loc),sh_keep_nos);
    re_keep = re_keep_br | re_keep_bu | re_keep_ge | re_keep_sh;
    ps.relay = ps.relay(re_keep,:);
end
nge = size(ps.gen);
nsh = size(ps.shunt);

% rebuild the bus and branch indices
max_bus_no  = max(ps.bus(:,1));
ps.bus_i    = sparse(ps.bus(:,1),1,(1:n)',max_bus_no,1);
if ~isempty(ps.branch)
    max_br_no   = max(ps.branch(:,C.br.id));
    ps.branch_i = sparse(ps.branch(:,C.br.id),1,(1:m)',max_br_no,1);
end
if ~isempty(ps.gen)
    ps.gen_i = sparse(ps.gen(:,1),1,(1:nge)',max_bus_no,1);
end
if ~isempty(ps.shunt)
    ps.shunt_i = sparse(ps.shunt(:,1),1,(1:nsh)',max_bus_no,1);
end


