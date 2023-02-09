function [CAERO] = sortCAEROArray(CAERO,PARAM,NODE,COORD,INFO)

ncaero = INFO.ncaero;

if ncaero > 0
    
    %fprintf(PARAM.FID, '\n\tSetting Aero database...');
    
    if PARAM.GRDPNT ~= 0
        CAERO.geo.ref_point = NODE.Coord(PARAM.GRDPNT, 1:3);
    else
        CAERO.geo.ref_point = zeros(1,3);
    end
    
    for n=1:ncaero
        
        if CAERO.CP(n) ~= 0
            
            j = find(COORD.ID == CAERO.CP(n));
            if isempty(j)
                error('Unable to find reference frame %d for CAERO entry %d.', CAERO.CP(n), CAERO.ID(n));
            else
                CAERO.CP(n) = j;
            end
            
            v = COORD.R(:,:,j) * [CAERO.geo.startx(n); CAERO.geo.starty(n); CAERO.geo.startz(n)];
            
            CAERO.geo.startx(n) = v(1) + COORD.Origin(j,1);
            CAERO.geo.starty(n) = v(2) + COORD.Origin(j,2);
            CAERO.geo.startz(n) = v(3) + COORD.Origin(j,3);
            
        end
        
    end
    
    CAERO.geo.nwing = ncaero;
    CAERO.geo.symetric = zeros(ncaero,1)';
    CAERO.geo.flap_vector = zeros(ncaero,max(CAERO.geo.nelem));
    
    fprintf(PARAM.FID,'done.');
    
end
end