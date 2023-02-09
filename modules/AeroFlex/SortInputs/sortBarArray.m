function [BAR,NODE,INFO] = sortBarArray(BAR,PBAR,NODE,MAT,COORD,INFO)

nbar = INFO.nbar;
ngrid = INFO.ngrid;
nrbe0 = INFO.nrbe0;


if nbar >0

    [labels, i] = unique(BAR.ID);
    
    if (length(labels) ~= nbar)
        
        n = [1 : nbar];
        error('Bar entries have duplicated labels: %d', BAR.ID(setdiff(n, i)));
        
    end
    
    BAR.Colloc = zeros(2, 3, nbar); % collocation point coords
    BAR.R = zeros(3, 3, 5, nbar);   % rotation material matrix
    BAR.D = zeros(6, 6, 2, nbar);   % stiffness matrix
    BAR.M = zeros(6, 6, 3, nbar);   % lumped body mass matrix associated to beams nodes
    % sort all database and use local indeces
    for n=1:nbar
        
        % substitute NODE ID with NODE index for rapid accessing
        
        i = find(NODE.ID == BAR.Conn(n,1));
        if isempty(i)
            error('Unable to determine node %d position for CBAR %d.',BAR.Conn(n,1), BAR.ID(n));
        end
        
        j = find(NODE.ID == BAR.Conn(n,2));
        if isempty(j)
            error('Unable to determine node %d position for CBAR %d.',BAR.Conn(n,2), BAR.ID(n));
        end
        
        % add the third node
        NODE.ID(ngrid + n) = int32(NODE.ID(ngrid) + n);
        NODE.CS(ngrid + n) = int32(0);
        NODE.CD(ngrid + n) = int32(0);
        index = [i, ngrid + n, j];
        BAR.Conn(n,:) = index;
        NODE.Coord(ngrid + n,:) = mean([NODE.Coord(i,:) ; NODE.Coord(j,:)] ,1);
        % substitute PID ID with PBAR index for rapid accessing
        
        i = find(PBAR.ID == BAR.PID(n));
        
        if length(i) ~= 1
            
            error('Unable to determine property card for CBAR %d.', BAR.ID(n));
            
        end
        
        BAR.PID(n) = i;
        
        % overwrite offset if g0 is integer
        
        if BAR.barg0(n)
            
            i = find(NODE.ID == int32(BAR.Orient(n,1)));
            if isempty(i)
                
                error('Unale to find node G0 %d in node database for bar %d.', ...
                    int32(BAR.Orient(n,1)), BAR.ID(n));
                
            end
            
            BAR.Orient(n,:) = NODE.Coord(i, :) - NODE.Coord(BAR.Conn(n,1), :);
            
        end
        
        % convert offset and orientation in BASIC reference frame
        
        % orient
        
        switch BAR.OffsetT(n)
            
            case {1,3,5,7}
                
                if NODE.CD(BAR.Conn(n,1)) ~= 0
                    
                    BAR.Orient(n,:) = (COORD.R(:, :, NODE.CD(BAR.Conn(n,1))) * BAR.Orient(n,:)')';
                    
                end
        end
        
        % offset
        
        switch BAR.OffsetT(n)
            
            case {1,2} % GGG , BGG
                
                % offset1
                if NODE.CD(BAR.Conn(n,1)) ~= 0
                    
                    BAR.Offset(n,1:3) = (COORD.R(:, :, NODE.CD(BAR.Conn(n,1))) * BAR.Offset(n,1:3)')';
                    
                end
                % offset3
                if NODE.CD(BAR.Conn(n,3)) ~= 0
                    
                    BAR.Offset(n,7:9) = (COORD.R(:, :, NODE.CD(BAR.Conn(n,3))) * BAR.Offset(n,7:9)')';
                    
                end
                
            case {3,4} % GGO , BGO
                
                % offset1
                if NODE.CD(BAR.Conn(n,1)) ~= 0
                    
                    BAR.Offset(n,1:3) = (COORD.R(:, :, NODE.CD(BAR.Conn(n,1))) * BAR.Offset(n,1:3)')';
                    
                end
                
                R = get_basic_coord(BAR.Orient(n,:), NODE.Coord(BAR.Conn(n,3),:), NODE.Coord(BAR.Conn(n,1),:));
                
                % offset3
                BAR.Offset(n,7:9) = (R * BAR.Offset(n,7:9)')';
                
                
            case {5,6} %GOG , BOG
                
                % offset3
                if NODE.CD(BAR.Conn(n,3)) ~= 0
                    
                    BAR.Offset(n,7:9) = (COORD.R(:, :, NODE.CD(BAR.Conn(n,3))) * BAR.Offset(n,7:9)')';
                    
                end
                
                R = get_basic_coord(BAR.Orient(n,:), NODE.Coord(BAR.Conn(n,3),:), NODE.Coord(BAR.Conn(n,1),:));
                
                % offset1
                BAR.Offset(n,1:3) = (R * BAR.Offset(n,1:3)')';
                
                
            case {7,8} %GOO , BOO
                
                R = get_basic_coord(BAR.Orient(n,:), NODE.Coord(BAR.Conn(n,3),:), NODE.Coord(BAR.Conn(n,1),:));
                
                % offset1
                BAR.Offset(n,1:3) = (R * BAR.Offset(n,1:3)')';
                % offset3
                BAR.Offset(n,7:9) = (R * BAR.Offset(n,7:9)')';
                
                
        end %switch
        
        % offset for the midline node
        BAR.Offset(n, 4:6) = get_midline_offset(NODE.Coord(BAR.Conn(n,1),:), ...
            NODE.Coord(BAR.Conn(n,2),:), NODE.Coord(BAR.Conn(n,3),:), BAR.Offset(n,1:3), BAR.Offset(n,7:9));
        
        % store global collocation point coordinates
        BAR.Colloc(:, :, n) = interp_colloc_pos(NODE.Coord(BAR.Conn(n,1),:) + BAR.Offset(n, 1:3), ...
            NODE.Coord(BAR.Conn(n,2),:) + BAR.Offset(n, 4:6), ...
            NODE.Coord(BAR.Conn(n,3),:) + BAR.Offset(n, 7:9));
        
        % build nodes and collocation point reference frames
        
        x1 = NODE.Coord(BAR.Conn(n,3),:) + BAR.Offset(n, 7:9) - (NODE.Coord(BAR.Conn(n,1),:) + BAR.Offset(n, 1:3));
        x2 = BAR.Orient(n,:);
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
            
            BAR.R(:,:,i, n) =  R;
            
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
    BAR = set_cbar_database(nbar, BAR, PBAR, MAT, NODE, COORD);
    %
    if nrbe0
        for n=ngrid+1:(ngrid + nbar)
            
            NODE.Aero.Coord(n).data = [];
            NODE.Aero.Index(n).data = [];
            
        end
    end
    
    ngrid = ngrid + nbar; % update node counter
    
end

INFO.ngrid = ngrid;

end
