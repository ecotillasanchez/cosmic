function ps = unify_generators(ps)
% usage: ps = unify_generators(ps)
% converts multiple generators connected to a bus into a single generator
% eduardo.cotilla-sanchez

% find buses with generators
C  = psconstants; 
ps = updateps(ps);

gen_buses       = unique(ps.gen(:,1));
remove_set      = [];

for i = 1:size(gen_buses,1)
    is_merge = ps.gen(:,1) == gen_buses(i);
    % look for multiple generators and merge them
    if sum(is_merge) > 1
        % mark the generators that will be kept/removed
        merge_set   = find(is_merge);
        keep        = merge_set(1);
        remove_set  = union(remove_set,merge_set(2:end));

        % calculate/update equivalent generator
        ps.gen(keep,C.gen.Pg)   = sum(ps.gen(is_merge,C.gen.Pg));
        ps.gen(keep,C.gen.Qg)   = sum(ps.gen(is_merge,C.gen.Qg));
        ps.gen(keep,C.gen.Qmax) = sum(ps.gen(is_merge,C.gen.Qmax));
        ps.gen(keep,C.gen.Qmin) = sum(ps.gen(is_merge,C.gen.Qmin));
        ps.gen(keep,C.gen.Pmax) = sum(ps.gen(is_merge,C.gen.Pmax));
        ps.gen(keep,C.gen.Pmin) = sum(ps.gen(is_merge,C.gen.Pmin));
        
        % assume linear cost
        ps.gen(keep,C.gen.cost) = sum(ps.gen(is_merge,C.gen.cost));
        ps.gencost = [];
    end
end

% remove the generators
ps.gen(remove_set,:) = [];
ps = updateps(ps);