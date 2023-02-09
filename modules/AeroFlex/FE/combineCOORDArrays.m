function [COORD,INFO] = combineCOORDArrays(COORD,COORD1,COORD2,INFO)

% create a unique COORD database
ncord1 = INFO.ncord1;
ncord2 = INFO.ncord2;

ntot = ncord1 + ncord2;

if ntot >0
    
    fprintf('\n\tSorting Coordinate System database...');
    
    COORD.ID = [COORD1.ID, COORD2.ID];
    COORD.Origin = zeros(ntot, 3);
    COORD.R = zeros(3, 3, ntot);
    
    COORD.Origin = [COORD1.Origin;
        COORD2.Origin];
    
    for n=1:ncord1
        
        COORD.R(:,:,n) = COORD1.R(:,:,n);
        
    end
    i = 0;
    
    for n=ncord1+1:ntot
        i = i + 1;
        COORD.R(:,:,n) = COORD2.R(:,:,i);
        
    end
    
    [COORD.ID, index] = sort(COORD.ID);
    [labels, i] = unique(COORD.ID);
    
    if (length(labels) ~= ntot)
        
        n = [1 : ntot];
        
        error('COORD1 and COORD2 entries have duplicated labels: %d', COORD.ID(setdiff(n, i)));
        
    end
    
    COORD.Origin = COORD.Origin(index,1:3);
    COORD.R = COORD.R(:,:,index);
    
    fprintf('done.');
    
end

INFO.ncord = ntot;
end
