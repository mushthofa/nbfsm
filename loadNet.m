function net = loadNet(filename)

edges = readtable(filename, 'ReadVariableNames', 1, 'FileType', 'Text');

nodes = []; 
node1 = edges.node1;
node2 = edges.node2;

nodes = union(nodes, node1); 
nodes = union(nodes, node2);

net.size = length(nodes);
net.nodes = nodes;

name2node = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i=1:length(nodes)
    name2node(nodes{i}) = i;
end

net.mat = single(zeros(net.size, net.size));
idxrow = cell2mat(values(name2node, node1));
idxcol = cell2mat(values(name2node, node2));
idx = sub2ind(size(net.mat), idxrow, idxcol);
net.mat(idx) = 1;

net.name2node = name2node;

