function [NODE, RBE0] = sortRBE0Array(NODE,RBE0,INFO)

nrbe0 = INFO.nrbe0;

if nrbe0
    
    for n=1:ngrid
        
        NODE.Aero.Coord(n).data = [];
        NODE.Aero.Index(n).data = [];
        
    end
    
    
    for n = 1:nrbe0
        
        i = find(RBE0.Master(n) == NODE.ID);
        if isempty(i)
            
            error('Unable to find master node %d in node database for aero rigid element %d.', ...
                RBE0.Master(n), RBE0.ID(n));
            
        end
        
        RBE0.Master(n) = i;
        [~, index] = intersect(NODE.ID, RBE0.Node(n).data);
        
        ne = length(RBE0.Node(n).data);
        
        if length(index) ~= ne
            
            j = setdiff([1:ne], index);
            error('Unable to find slave aero node %d for aero rigid set %d.', ...
                RBE0.Node(n).data(j(1)), RBE0.ID(n));
            
        end
        
        RBE0.Node(n).data = index;
        % add to master node, slave aero points
        coord = NODE.Coord(RBE0.Node(n).data, 1:3) - repmat(NODE.Coord(RBE0.Master(n), 1:3), ne, 1);
        NODE.Aero.Coord(RBE0.Master(n)).data = [NODE.Aero.Coord(RBE0.Master(n)).data, coord'];
        NODE.Aero.Index(RBE0.Master(n)).data = [NODE.Aero.Index(RBE0.Master(n)).data, RBE0.Node(n).data];
        
    end
    
end
end