function ps = update_load_freq_source(ps)
% usage: ps = update_load_freq_source(ps)
% this function assign a generator freqency signal to each bus. The
% frequency source is one of the nearest generators with largest inertia M

C  = psconstants;
ps.neighbors = get_neighbors(ps.bus(:,1),ps.branch(:,1:2)); % get neighbors of each bus
ps.Ebus = electrical_distance(ps.Ybus);                     % electrical distance matrix
is_load = ps.bus_i(ps.shunt(:,1));                          % find the load bus ids
n = size(is_load,1);
Near_Large_Gen2Load = nan(n,1);

if size(ps.bus(:,1))<30
    n_gen_limit = 1;
else
    n_gen_limit = 2;  % size of result pool for each load bus.(ie. the nearest 2 generators)
end

% get the nearest generator id for each load bus
for i = 1:n
    j = is_load(i);
    Near_Large_Gen2Load(i,1) = get_nearest_gen(ps,j,n_gen_limit); 
end
ps.shunt(:,C.sh.near_gen) = Near_Large_Gen2Load;

end