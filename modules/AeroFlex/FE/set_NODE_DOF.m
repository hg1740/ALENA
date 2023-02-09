function [NODE,INFO] = set_NODE_DOF(PARAM,NODE,SPC,INFO)

% Initialise the DOF structure array
NODE.DOF = int32(zeros(INFO.ngrid, 6));

% If SPCs exist take their DOF into account
if PARAM.SPC
    
    % Check SPC IDs match
    spcindex = find(SPC.ID == PARAM.SPC);
    if isempty(spcindex)
        error('Unable to find SPC set %d as required.', PARAM.SPC);
    else
        
        for i = spcindex
            NODE.DOF(SPC.Nodes(i).list, SPC.DOF(i).list) = int32(1);
        end
        
        set = (1:6);
        ndof = 0;
        
        % Count DOF
        for n=1:INFO.ngrid
            
            if NODE.Index(n)
                
                index = find(NODE.DOF(n, :));
                index = setdiff(set, index);
                NODE.DOF(n, :) = int32(zeros(1,6));
                NODE.DOF(n, index) = int32((ndof+1 : ndof+length(index)));
                ndof = ndof + length(index);
                
            end
        end
    end 
else
    
    ndof = 0;
    
    % Count DOF
    for n = 1:INFO.ngrid
        
        if NODE.Index(n)
            
            NODE.DOF(n, :) = int32(zeros(1,6));
            NODE.DOF(n, :) = int32((ndof+1 : ndof+6));
            ndof = ndof + 6;
            
        end
    end
    
end 

INFO.ndof = ndof;

end