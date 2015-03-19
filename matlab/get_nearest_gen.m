function Near_LargeM_Gen = get_nearest_gen(ps,source_bus_id,n_gen_limit)
% usage: Near_LargeM_Gen = get_nearest_gen(ps,source_bus_id,n_gen_limit)
% This function gets the relatively nearest generator based on electrical
% distance, and picks up one with largest inertia M
% Algorithm: Dijkstra algorithm in social network search


C  = psconstants;
neighbors = ps.neighbors;
num_node = length(neighbors);
mac_buses   = ps.gen(:,C.gen.bus);
elec_dist   = ps.Ebus;

dist = Inf(1, num_node);        % shortest distance to source up to now
prev = zeros(1, num_node);      % previous node in optimal path from source
visited = false(1, num_node);   % assume no one has been visited

dist(source_bus_id) = 0;
nearestGen = zeros(1, num_node); % you could possibly find num_node generators if every bus has one
genInd = 0;

for ii = 1:num_node
    temp_dist = dist;
    temp_dist(visited) = Inf;
    [min_dist, index] = min(temp_dist);
    visited(index) = true;
    
    if ismember(ps.bus(index,C.bu.id),mac_buses) % whether ind is a gen bus
        genInd = genInd + 1;
        nearestGen(1,genInd) = ps.bus(index,C.bu.id);
    end
    
    if genInd >= n_gen_limit          % find n_gen_limit gens that are close to source
        break;
    end
    
    if dist(index) == Inf             % index is in different network
        break;
    end
    
    neighbors_j = neighbors{index};
    for jj = 1 : length(neighbors_j)
        n = neighbors_j(jj);
        neighborDist = elec_dist(n,index); % electrical distance between n and ind
        newDist = min_dist + neighborDist;
        if newDist < dist(n)
            dist(n) = newDist;
            prev(n) = index;
        end
    end
end

% from the 3 or n_gen_limit of gernerators, pick up one with largest
% inertia M
all_gens = nan(n_gen_limit,2);
all_gens(:,1) = nearestGen(1,1:n_gen_limit)';
all_gens_loc = ismember(mac_buses,all_gens(:,1));

if isempty(all_gens_loc)
    all_gens(:,2) = NaN;
else
    all_gens(:,2) = ps.mac(all_gens_loc,C.ma.M);     % C.ge.Pg if there's no machine data
end
[~,near_large_gen_index] = max(all_gens(:,2));
Near_LargeM_Gen = all_gens(near_large_gen_index,1);

end

