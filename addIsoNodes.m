function net = addIsoNodes(net, nnodes)
% Adding dummy isolated nodes to a network

    nsize = length(nnodes);

    zz = zeros(net.size, nsize);
    net.mat = [net.mat zz; zz' eye(nsize)];
    
    net.nodes = [net.nodes; nnodes];
    for i=1:nsize
        net.name2node(nnodes{i}) = net.size + i;
    end
    net.size = net.size + nsize;
end
