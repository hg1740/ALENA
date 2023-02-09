function [COORD2] = sortCOORD2Array(COORD2,INFO)

CORDIT = 10;

ncord2 = INFO.ncord2;

if ncord2>0
    
    [COORD2.ID, index] = sort(COORD2.ID);
    [labels, i] = unique(COORD2.ID);
    
    if (length(labels) ~= ncord2)
        n = [1 : ncord2];
        error(['CORD2R entries have duplicated labels: ', num2str(COORD2.ID(setdiff(n, i)))]);
    end
    
    COORD2.Nodes = COORD2.Nodes(:,:,index);
    COORD2.RID = COORD2.RID(index);
    
    for n=1:ncord2
        if COORD2.RID(n)~=0
            index = find(COORD2.ID == COORD2.RID(n));
            if isempty(index)
                error(['Unable to find CORD2R ', num2str(COORD2.RID(n)), ' declared in CORD2R ', num2str(COORD2.ID(n))]);
            else
                COORD2.RID(n) = index;
            end
        end
    end
    
    % add two fields
    COORD2.R = zeros(3, 3, ncord2);
    COORD2.Origin = zeros(ncord2,3);
    TAG = zeros(ncord2,1);
    
    % loop on master COORD
    for n=1:ncord2
        if COORD2.RID(n) == 0
            TAG(n) = 1;
            [COORD2.Origin(n,1:3), COORD2.R(:,:,n)] = cord2r_data(COORD2, n);
        end
    end
    
    % loop on 1st and 2nd type dependence
    for k=1:CORDIT
        for n=1:ncord2
            m = COORD2.RID(n);
            if m ~= 0
                if TAG(m) == 1 && TAG(n) == 0
                    TAG(n) = 1;
                    [COORD2.Origin(n,1:3), COORD2.R(:,:,n)] = cord2r_data(COORD2, n);
                    COORD2.Origin(n,1:3) = COORD2.Origin(m,1:3) + (COORD2.R(:,:,m) * COORD2.Origin(n,1:3)')';
                    COORD2.R(:,:,n) = COORD2.R(:,:,m) * COORD2.R(:,:,n);
                end
            end
        end
        if sum(TAG) == ncord2
            break;
        end
    end
    
    if k == CORDIT && sum(TAG) < ncord2
        error('Unable to solve CORD2R dependency within 10 iterations. Please change CORD2R dependency.');
    end
end

end

function [O, R] = cord2r_data(COORD2, n)

O = COORD2.Nodes(1,1:3,n);
x3 = COORD2.Nodes(2,1:3,n) - COORD2.Nodes(1,1:3,n);
x1 = COORD2.Nodes(3,1:3,n) - COORD2.Nodes(1,1:3,n);
x2 = cross(x3, x1);
if det([x1; x2; x3]) == 0
    error('COORD2 %d definition has collinear points.', COORD2.ID(n));
end
x3 = x3 ./ norm(x3);
x2 = x2 ./ norm(x2);
x1 = cross(x2, x3);
x1 = x1 ./ norm(x1);
R = [x1', x2', x3'];

end