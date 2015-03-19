clc
clear all

% test update_load_freq_source.m
C  = psconstants;
ps = updateps(case39_ps);
% ps = replicate_case(ps,2);
% ps = unify_generators(ps);  

% data = load('case3120sp_ps.mat');   
% ps = data.ps;
% ps = unify_generators(ps);     % converts multiple generators connected to a bus into a single generator
opt = psoptions;
[ps.Ybus,ps.Yf,ps.Yt] = getYbus(ps,false);
ps.Ebus = electrical_distance(ps.Ybus);                     % electrical distance matrix
ps.neighbors = get_neighbors(ps.bus(:,1),ps.branch(:,1:2)); % get neighbors of each bus
is_load = ps.bus_i(ps.shunt(:,1));                          % find the load bus ids
n = size(is_load,1);
Near_Large_Gen2Load = nan(n,1);

if size(ps.bus(:,1))<30
    n_gen_limit = 1;
else
    n_gen_limit = 2;  % size of result pool for each load bus.(ie. the nearest 3 generators)
end

% get the nearest generator id for each load bus
for i = 1:n
    j = is_load(i);
    Near_Large_Gen2Load(i,1) = get_nearest_gen(ps,j,n_gen_limit); 
end
ps.shunt(:,C.sh.near_gen) = Near_Large_Gen2Load;

Load2Gen = [ps.shunt(:,1),ps.shunt(:,C.sh.near_gen)];
disp('|Load#|NearGen#|');
disp(Load2Gen);
