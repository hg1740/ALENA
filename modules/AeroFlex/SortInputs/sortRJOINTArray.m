function [RJOINT] = sortRJOINTArray(RJOINT,NODE,INFO)

nrjoint = INFO.nrjoint;

if nrjoint > 0
    fprintf('\n\tSorting Rigid Joints...');
    
    RJOINTID = sort(RJOINT.ID);
    if ~isempty(find(diff(RJOINTID)==0,1))
        error('\n\t: duplicated RJoint.');
    end
    
    for i = 1 : nrjoint
        if RJOINT.GA(i) == RJOINT.GB(i)
            error('\n\tRigid joint grids must be different %d. ',RJOINT.ID(i));
        end
        RJointDOF2 = RJOINT.DOFC(i).list;
        
        RJOINT.DOF(i,1).data = NODE.DOF(NODE.ID == RJOINT.GA(i),RJointDOF2);
        RJOINT.DOF(i,2).data = NODE.DOF(NODE.ID == RJOINT.GB(i),RJointDOF2);
        RJOINT.Node(i,1) = find(NODE.ID == RJOINT.GA(i));
        RJOINT.Node(i,2) = find(NODE.ID == RJOINT.GB(i));
        
        x1 = RJOINT.Orient(i,:);
        x2 = [0,1,0];
        x3 = cross(x1,x2);
        x2 = cross(x3,x1);
        x1 = x1 ./ norm(x1);
        x2 = x2 ./ norm(x2);
        x3 = x3 ./ norm(x3);
        R = [x1', x2', x3'];
        RJOINT.R(:,:,1,i) = R;
        RJOINT.R(:,:,2,i) = R;
    end
    
    fprintf('done.');
end
end
