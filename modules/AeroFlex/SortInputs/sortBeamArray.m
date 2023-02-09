function [BEAM,NODE,INFO] = sortBeamArray(BEAM,PBEAM,NODE,MAT,COORD,INFO)

nbeam = INFO.nbeam;
ngrid = INFO.ngrid;
nrbe0 = INFO.nrbe0;

if nbeam >0
        
    [labels, i] = unique(BEAM.ID);
    
    if (length(labels) ~= nbeam)  
        n = [1 : nbeam];
        error('Bar entries have duplicated labels: %d', BEAM.ID(setdiff(n, i)));
    end
    
    BEAM.Colloc = zeros(2, 3, nbeam); % collocation point coords
    BEAM.R = zeros(3, 3, 5, nbeam);   % rotation material matrix
    BEAM.D = zeros(6, 6, 2, nbeam);   % stiffness matrix
    BEAM.M = zeros(6, 6, 3, nbeam);   % lumped body mass matrix associated to beams nodes
    
    % sort all database and use local indeces
    for n=1:nbeam
        
        % substitute NODE ID with NODE index for rapid accessing
        
        i = find(NODE.ID == BEAM.Conn(n,1));
        if isempty(i)
            error('Unable to determine node %d position for CBAR %d.',BEAM.Conn(n,1), BEAM.ID(n));
        end
        
        j = find(NODE.ID == BEAM.Conn(n,2));
        if isempty(j)
            error('Unable to determine node %d position for CBAR %d.',BEAM.Conn(n,2), BEAM.ID(n));
        end
        
        % add the third node
        NODE.ID(ngrid + n) = int32(NODE.ID(ngrid) + n);
        NODE.CS(ngrid + n) = int32(0);
        NODE.CD(ngrid + n) = int32(0);
        index = [i, ngrid + n, j];
        BEAM.Conn(n,:) = index;
        % Create a node at the mid-point of the beam
        NODE.Coord(ngrid + n,:) = mean([NODE.Coord(i,:) ; NODE.Coord(j,:)] ,1);
        % Create a node at the mid-point of the curved beam
        % Could add this as a 3d vector of curvature
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%         curvature = 1/100;
%         beam_length = norm(NODE.Coord(i,:) - NODE.Coord(j,:));
%         Midpoint = mean([NODE.Coord(i,:) ; NODE.Coord(j,:)] ,1);
%         NormalBeam = [-NODE.Coord(j,2)+NODE.Coord(i,2),NODE.Coord(j,1)-NODE.Coord(i,1),0.5*(NODE.Coord(j,3)+NODE.Coord(i,3))];
%         NODE.Coord(ngrid + n,:) = Midpoint + (1/beam_length)*(1/curvature - 0.5*sqrt(4/(curvature^2)-beam_length^2))*NormalBeam;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % substitute PID ID with PBAR index for rapid accessing
        
        i = find(PBEAM.ID == BEAM.PID(n));
        
        if length(i) ~= 1
            error('Unable to determine property card for CBEAM %d.', BEAM.ID(n));
        end
        
        BEAM.PID(n) = i;
        
        % overwrite offset if g0 is integer
        
        if BEAM.beamg0(n)
            
            i = find(NODE.ID == int32(BEAM.Orient(n,1)));
            if isempty(i)
                
                error('Unale to find node G0 %d in node database for bar %d.', ...
                    int32(BEAM.Orient(n,1)), BEAM.ID(n));
                
            end
            
            BEAM.Orient(n,:) = NODE.Coord(i, :) - NODE.Coord(BEAM.Conn(n,1), :);
            
        end
        
        % convert offset and orientation in BASIC reference frame
        
        % orient
        
        switch BEAM.OffsetT(n)
            
            case {1,3,5,7}
                
                if NODE.CD(BEAM.Conn(n,1)) ~= 0
                    
                    BEAM.Orient(n,:) = (COORD.R(:, :, NODE.CD(BEAM.Conn(n,1))) * BEAM.Orient(n,:)')';
                    
                end
        end
        
        % offset
        
        switch BEAM.OffsetT(n)
            
            case {1,2} % GGG , BGG
                
                % offset1
                if NODE.CD(BEAM.Conn(n,1)) ~= 0
                    
                    BEAM.Offset(n,1:3) = (COORD.R(:, :, NODE.CD(BEAM.Conn(n,1))) * BEAM.Offset(n,1:3)')';
                    
                end
                % offset3
                if NODE.CD(BEAM.Conn(n,3)) ~= 0
                    
                    BEAM.Offset(n,7:9) = (COORD.R(:, :, NODE.CD(BEAM.Conn(n,3))) * BEAM.Offset(n,7:9)')';
                    
                end
                
            case {3,4} % GGO , BGO
                
                % offset1
                if NODE.CD(BEAM.Conn(n,1)) ~= 0
                    
                    BEAM.Offset(n,1:3) = (COORD.R(:, :, NODE.CD(BEAM.Conn(n,1))) * BEAM.Offset(n,1:3)')';
                    
                end
                
                R = get_basic_coord(BEAM.Orient(n,:), NODE.Coord(BEAM.Conn(n,3),:), NODE.Coord(BBEAMAR.Conn(n,1),:));
                
                % offset3
                BEAM.Offset(n,7:9) = (R * BEAM.Offset(n,7:9)')';
                
                
            case {5,6} %GOG , BOG
                
                % offset3
                if NODE.CD(BEAM.Conn(n,3)) ~= 0
                    
                    BEAM.Offset(n,7:9) = (COORD.R(:, :, NODE.CD(BEAM.Conn(n,3))) * BEAM.Offset(n,7:9)')';
                    
                end
                
                R = get_basic_coord(BEAM.Orient(n,:), NODE.Coord(BEAM.Conn(n,3),:), NODE.Coord(BEAM.Conn(n,1),:));
                
                % offset1
                BEAM.Offset(n,1:3) = (R * BEAM.Offset(n,1:3)')';
                
                
            case {7,8} %GOO , BOO
                
                R = get_basic_coord(BEAM.Orient(n,:), NODE.Coord(BEAM.Conn(n,3),:), NODE.Coord(BEAM.Conn(n,1),:));
                
                % offset1
                BEAM.Offset(n,1:3) = (R * BEAM.Offset(n,1:3)')';
                % offset3
                BEAM.Offset(n,7:9) = (R * BEAM.Offset(n,7:9)')';
                
                
        end %switch
        
        % offset for the midline node
        BEAM.Offset(n, 4:6) = get_midline_offset(NODE.Coord(BEAM.Conn(n,1),:), ...
            NODE.Coord(BEAM.Conn(n,2),:), NODE.Coord(BEAM.Conn(n,3),:), BEAM.Offset(n,1:3), BEAM.Offset(n,7:9));
        %%%%%%%%%%%%%%%%%%%%%%%%%
        BEAM.Offset(n, 4:6) = 0.5*(BEAM.Offset(n,1:3) + BEAM.Offset(n,7:9));
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % store global collocation point coordinates
        BEAM.Colloc(:, :, n) = interp_colloc_pos(NODE.Coord(BEAM.Conn(n,1),:) + BEAM.Offset(n, 1:3), ...
            NODE.Coord(BEAM.Conn(n,2),:) + BEAM.Offset(n, 4:6), ...
            NODE.Coord(BEAM.Conn(n,3),:) + BEAM.Offset(n, 7:9));
        
        % build nodes and collocation point reference frames
        
        x1 = NODE.Coord(BEAM.Conn(n,3),:) + BEAM.Offset(n, 7:9) - (NODE.Coord(BEAM.Conn(n,1),:) + BEAM.Offset(n, 1:3));
        x2 = BEAM.Orient(n,:);
        x3 = cross(x1, x2);
        if norm(x3) == 0
            
            error('Wrong Orientation direction definition for Bar %d. Please check', BAR.ID(n));
            
        end
        x2 = cross(x3, x1);
        
        x1 = x1 ./ norm(x1);
        x2 = x2 ./ norm(x2);
        x3 = x3 ./ norm(x3);
        R= [x1', x2', x3'];
        % when assembled all the reference frame belonging to the BAR are the same
        for i=1:5 % 3 reference for 3 nodes + 2 reference for two collocation points
            
            BEAM.R(:,:,i, n) =  R;
            
        end
        
        %WARNING: offset are expressend in the basic coordinates system which
        %         correspond before deformation to each node reference frame
        % get local offset vector
        %BAR.Offset(n, 1:3) = (R' * BAR.Offset(n, 1:3)')';
        %BAR.Offset(n, 4:6) = (R' * BAR.Offset(n, 4:6)')';
        %BAR.Offset(n, 7:9) = (R' * BAR.Offset(n, 7:9)')';
        
        % assembly stiffnes matrix
        %
    end % bar loop
    % set stiffness and mass database
    BEAM = set_cbeam_database(nbeam, BEAM, PBEAM, MAT, NODE, COORD);
    %
    if nrbe0
        for n=ngrid+1:(ngrid + nbeam)
            
            NODE.Aero.Coord(n).data = [];
            NODE.Aero.Index(n).data = [];
            
        end
    end
    
    ngrid = ngrid + nbeam; % update node counter
    
end

INFO.ngrid = ngrid;

end
