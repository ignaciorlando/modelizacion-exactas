function [ feasible ] = is_feasible(MIN_TIME, adjacency, costs, begin_end)

IN = 1; DEST = -1;

graph = sparse(length(costs),length(costs));

[row_i,~] = find(adjacency == IN);
[row_e,~] = find(adjacency == DEST);

for i = 1:length(costs)
    graph(row_i(i),row_e(i)) = costs(i);
end

orig = find(begin_end == IN);
dest = find(begin_end == DEST);
[dist,~,~] = graphshortestpath(graph, orig, dest);

feasible = dist <= MIN_TIME;
end

