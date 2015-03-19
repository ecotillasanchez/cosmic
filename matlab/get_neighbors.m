function neighbors = get_neighbors(nodes,links)
% build the "neighbors" structure, which is used by "connections2clusters"
% usage: neighbors = get_neighbors(nodes,links)

% renumber to get sequential numbering
n = length(nodes);
e2i = sparse(nodes,1,(1:n)',max(nodes),1);
nodes = e2i(nodes);
links(:,1) = e2i(links(:,1));
links(:,2) = e2i(links(:,2));

% figure out who the neighbors are
neighbors = cell(n,1);
for i = 1:n
    nodeset1 = links(links(:,1)==nodes(i),2);
    nodeset2 = links(links(:,2)==nodes(i),1);
    neighbors{i} = union(nodeset1,nodeset2);
end
