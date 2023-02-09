function [CBUSH] = sortCBUSHArray(CBUSH,PBUSH,NODE,INFO)

ncbush = INFO.ncbush;

if ncbush>0
    fprintf('\n\tSorting Cbush database...');
    CBUSHID = sort(CBUSH.ID);
    if ~isempty(find(diff(CBUSHID)==0,1))
        error('\n\t: duplicated CBush.');
    end
    
    CBUSH.R = zeros(3, 3, 2, ncbush);
    
    for i = 1 : ncbush
        if CBUSH.Node(i,1) == CBUSH.Node(i,2)
            error('\n\tCBush grids must be different %d. ',CBUSH.ID(i));
        end
        
        CbushDOF2 = [1,2,3,4,5,6];

        CBUSH.DOF(i,1).data = NODE.DOF(NODE.ID == CBUSH.Node(i,1),CbushDOF2);
        CBUSH.DOF(i,2).data = NODE.DOF(NODE.ID == CBUSH.Node(i,2),CbushDOF2);
        idx = find(PBUSH.ID == CBUSH.PID(i));
        CBUSH.K(i).data = PBUSH.K(idx,:);
        CBUSH.Node(i,1) = find(NODE.ID == CBUSH.Node(i,1));
        CBUSH.Node(i,2) = find(NODE.ID == CBUSH.Node(i,2));
        
        % Store the Rotation matrix associated the hinge orientation
        x1 = CBUSH.Orient(i,:);
        x2 = [0,1,0];
        x3 = cross(x1,x2);
        x2 = cross(x3,x1);
        x1 = x1 ./ norm(x1);
        x2 = x2 ./ norm(x2);
        x3 = x3 ./ norm(x3);
        R = [x1', x2', x3'];
        CBUSH.R(:,:,1,i) = R;
        CBUSH.R(:,:,2,i) = R;
    end
    
    fprintf('done.');
end
end
