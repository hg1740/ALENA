function [RBAR] = sortRBARArray(RBAR,NODE,INFO)

nrbar = INFO.nrbar;

if nrbar > 0
    fprintf('\n\tSorting Rigid Bar elements...');
    
    RBARID = sort(RBAR.ID);
    if ~isempty(find(diff(RBARID)==0,1))
        error('\n\t: duplicated RBar.');
    end
    
    for i = 1 : nrbar
        if RBAR.GA(i) == RBAR.GB(i)
            error('\n\tRigid bar grids must be different %d. ',RBAR.ID(i));
        end
        RBarDOF2 = RBAR.DOFC(i).list;
        
        RBAR.DOF(i,1).data = NODE.DOF(NODE.ID == RBAR.GA(i),RBarDOF2);
        RBAR.DOF(i,2).data = NODE.DOF(NODE.ID == RBAR.GB(i),RBarDOF2);
        master = find(NODE.ID == RBAR.GA(i));
        slave = find(NODE.ID == RBAR.GB(i));
        RBAR.Node(i,1) = master;
        RBAR.Node(i,2) = slave;
         
        r = NODE.Coord(slave,:) - NODE.Coord(master,:);
        
        RBAR.R(i).data = crossm(r);
        RBAR.rel(i).data = r;
    end
    
    fprintf('done.');
end
end
