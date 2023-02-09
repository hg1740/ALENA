function [CAERO] = sortINTERPdata(CAERO,NODE,INFO)

% set interpolation database

nset = INFO.nset;

if nset
    
    index = unique(CAERO.Set.ID);
    
    if length(index) ~= nset
        n = setdiff(CAERO.Set.ID, index);
        error('Duplicated ID for structural SET1/AELIST %d.', n(1));
        
    end
    
end

for n=1:nset
    
    [~, index_n, index] = intersect(NODE.ID, CAERO.Set.Node(n).data,'legacy');
    
    ne = length(CAERO.Set.Node(n).data);
    
    if length(index_n) ~= ne
        
        j = setdiff([1:ne], index);
        error('Unable to find node %d for structural set %d.', ...
            CAERO.Set.Node(n).data(j(1)), CAERO.Set.ID(n));
        
    end
    
    CAERO.Set.Node(n).data = index_n';
    
end

end