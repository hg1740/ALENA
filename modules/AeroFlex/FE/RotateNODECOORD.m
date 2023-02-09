function [NODE] = RotateNODECOORD(NODE,COORD)

% define GRID in GLOBAL reference frame
index = find (NODE.CS);

for n=1:length(index)
    
    node = index(n);
    i = find(COORD.ID == NODE.CS(node));
    if isempty(i)
        
        error('Unable to find coordinate system %d for node %d definition.', ...
            NODE.CS(node), NODE.ID(node));
        
    end
    NODE.Coord(node, 1:3) = COORD.Origin(i,1:3) + (COORD.R(:,:,i) * NODE.Coord(node, 1:3)')';
    NODE.CS(node) = i; % overwrite id with index for rapid accessing
    
end

index = find(NODE.CD);
for n=1:length(index)
    
    node = index(n);
    i = find(COORD.ID == NODE.CD(node));
    
    if isempty(i)
        
        error('Unable to find coordinate system %d for node %d definition.', ...
            NODE.CD(node), NODE.ID(node));
        
    end
    
    NODE.CD(node) = i; % overwrite with index for rapid accessing
    
end

end