function [NODE] = sortNODEArray(NODE,INFO)

% sort NODE database
ngrid = INFO.ngrid;
if ngrid >0
    
    %fprintf('\n\tSorting Node database...');
    
    [NODE.ID, index] = sort(NODE.ID);
    [labels, i] = unique(NODE.ID);
    % check for duplicated labels
    if (length(labels) ~= ngrid)
        
        n = [1 : ngrid];
        dof = NODE.ID(setdiff(n, i));
        
        for k=1:length(dof)
            fprintf('\n\t### Warning: duplicated labels for grid: %d.', NODE.ID(dof(k)));
        end
        error('Grid entries have duplicated labels.');
    end
    
    NODE.CS = NODE.CS(index);
    NODE.Coord = NODE.Coord(index,1:3);
    NODE.CD = NODE.CD(index);
    %NODE.PS = NODE.PS(index);
    
    %fprintf('done.');
    
end

end