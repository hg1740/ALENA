% The following function combines the conm1 and conm2 entries as well as 
% any conms associated to the same node 

function [CONM,INFO] = sortCONMArray(CONM1,CONM2,NODE,COORD,PARAM,INFO)

 ncom1 = INFO.ncom1;
 ncom2 = INFO.ncom2;
 

 CONM = [];
    mtot = ncom1 + ncom2;
    
    if mtot > 0
        
        if ncom1>0
            CONM1.M = PARAM.WTMASS .* CONM1.M;
        end
        
        if ncom2 > 0
            CONM2.M = PARAM.WTMASS .* CONM2.M;
        end
        
        CONM.Offset = zeros(mtot, 3);
                
        % CONM1 struct
        i = 0;
        for n = 1 : ncom1
            
            i = i + 1;
            % set NODE index for rapid accessing
            j = find(CONM1.Node(i) == NODE.ID);
            
            if isempty(j)
                
                error('Unable to find node %d in Node database for CONM1 %d.', CONM1.Node(i), CONM1.ID(i));
                
            end
            
            CONM.Node(n) = j; % set index for rapid accessing
            
            CONM.M(:,:,n) = zeros(6,6);
            CONM.M(:,:,n) = CONM1.M(:,:,i);
            
            if CONM1.CID(i) > 0
                
                j = find(CONM1.CID(i) == COORD.ID);
                
                if isempty(j)
                    error('Unable to find reference frame %d in COORD database for CONM1 %d.', CONM1.CID(i), CONM1.ID(i));
                end
                
                zcg = CONM.M(5,1,n) / CONM.M(1,1,n);
                ycg = -CONM.M(6,1,n) / CONM.M(1,1,n);
                xcg = CONM.M(6,2,n) / CONM.M(1,1,n);
                v = [xcg; ycg; zcg];
                
                v = COORD.R(:,:,j) * v; % get mass cg in global reference frame
                
                CONM.Offset(n,:) = v'; % store offset
                
                S = zeros(3,3);
                S(2,1) =  CONM.M(1,1,n) * v(3);
                S(3,1) = -CONM.M(1,1,n) * v(2);
                S(1,2) = -S(2,1);
                S(3,2) =  CONM.M(1,1,n) * v(1);
                S(1,3) = -S(3,1);
                S(2,3) = -S(3,2);
                
                CONM.M(4:6, 1:3, n) = S;
                CONM.M(1:3, 4:6, n) = S';
                % mass matrix in basic reference frame in the node
                RJ = COORD.R(:,:,j) * sqrt(CONM1.M(4:6,4:6,i));
                CONM.M(4:6, 4:6, n) = RJ * RJ' - ...
                    CONM1.M(1,1,i) * crossm(CONM.Offset(n,:)) * crossm(CONM.Offset(n,:));
                
            end
            
            if CONM1.CID(n) == 0
                
                zcg = CONM.M(5,1,n) / CONM.M(1,1,n);
                ycg = -CONM.M(6,1,n) / CONM.M(1,1,n);
                xcg = CONM.M(6,2,n) / CONM.M(1,1,n);
                v = [xcg; ycg; zcg];
                
                CONM.Offset(n,:) = v';
                
                S = zeros(3,3);
                S(2,1) =  CONM.M(1,1,n) * v(3);
                S(3,1) = -CONM.M(1,1,n) * v(2);
                S(1,2) = -S(2,1);
                S(3,2) =  CONM.M(1,1,n) * v(1);
                S(1,3) = -S(3,1);
                S(2,3) = -S(3,2);
                
                CONM.M(4:6, 1:3, n) = S;
                CONM.M(1:3, 4:6, n) = S';
                CONM.M(4:6, 4:6, n) = CONM1.M(4:6,4:6,i) - CONM1.M(1,1,i) * crossm(CONM.Offset(n,:)) * ...
                    crossm(CONM.Offset(n,:));
                
            end
            
        end
        
        % CONM2 struct
        i = 0;
        for n = (ncom1+1) : mtot
            
            i = i+1;
            j = find(CONM2.Node(i) == NODE.ID);
            
            if isempty(j)
                
                error('Unable to find node %d in Node database for CONM2 %d.', CONM2.Node(i), CONM2.ID(i));
            end
            
            CONM.Node(n) = j; % set NODE index for rapid accessing
            
            if CONM2.CID(n) == 0
                
                CONM.Offset(n,:) = CONM2.Offset(n, :);
                CONM.M(:,:,n) = zeros(6,6);
                CONM.M(1:3,1:3,n) = CONM2.M(1:3,1:3,i);
                % set static moment
                S = zeros(3,3);
                S(2,1) =  CONM.M(1,1,n) * CONM2.Offset(i, 3);
                S(3,1) = -CONM.M(1,1,n) * CONM2.Offset(i, 2);
                S(1,2) = -S(2,1);
                S(3,2) =  CONM.M(1,1,n) * CONM2.Offset(i, 1);
                S(1,3) = -S(3,1);
                S(2,3) = -S(3,2);
                
                CONM.M(4:6, 1:3, n) = S;
                CONM.M(1:3, 4:6, n) = S';
                % get Inertia in a reference frame with origin in the NODE and parallel
                % to the basic reference frame
                CONM.M(4:6, 4:6, n) = CONM2.M(4:6,4:6,i) - CONM2.M(1,1,i) .* ...
                    (crossm(CONM2.Offset(i, 1:3)) * crossm(CONM2.Offset(i, 1:3)));
                
            end
            
            if CONM2.CID(i) > 0
                
                j = find(CONM2.CID(i) == COORD.ID);
                
                if isempty(j)
                    
                    error('Unable to find reference frame %d in COORD database for CONM2 %d.', ...
                        CONM2.CID(i), CONM2.ID(i));
                    
                end
                % set offset in the global reference frame
                CONM2.Offset(i, 1:3) = (COORD.R(:,:,j) * CONM2.Offset(i, 1:3)')';
                CONM.Offset(n,:) = CONM2.Offset(i, :);
                CONM.M(:,:,n) = zeros(6,6);
                CONM.M(1:3,1:3,n) = CONM2.M(1:3,1:3,i);
                % set static moment
                S = zeros(3,3);
                S(2,1) =  CONM.M(1,1,n) * CONM2.Offset(i, 3);
                S(3,1) = -CONM.M(1,1,n) * CONM2.Offset(i, 2);
                S(1,2) = -S(2,1);
                S(3,2) =  CONM.M(1,1,n) * CONM2.Offset(i, 1);
                S(1,3) = -S(3,1);
                S(2,3) = -S(3,2);
                
                CONM.M(4:6, 1:3, n) = S;
                CONM.M(1:3, 4:6, n) = S';
                % get Inertia in a reference frame with origin in the NODE and parallel to the basic reference frame
                RJ = COORD.R(:,:,j) * sqrt(CONM2.M(4:6,4:6,i));
                CONM.M(4:6, 4:6, n) = RJ * RJ' - ...
                    CONM2.M(1,1,i) .* (crossm(CONM2.Offset(i, 1:3)) * crossm(CONM2.Offset(i, 1:3)));
                
            end
            
            if CONM2.CID(i) < 0
                
                % set offset in the global reference frame as the difference between the CG and the NODE
                CONM2.Offset(i, 1:3) = CONM2.Offset(i, 1:3) - NODE.Coord(CONM.Node(n),1:3);
                CONM.Offset(n,:) = CONM2.Offset(i, 1:3);
                CONM.M(:,:,n) = zeros(6,6);
                CONM.M(1:3,1:3,n) = CONM2.M(1:3,1:3,i);
                % set static moment
                S = zeros(3,3);
                S(2,1) =  CONM.M(1,1,n) * CONM2.Offset(i, 3);
                S(3,1) = -CONM.M(1,1,n) * CONM2.Offset(i, 2);
                S(1,2) = -S(2,1);
                S(3,2) =  CONM.M(1,1,n) * CONM2.Offset(i, 1);
                S(1,3) = -S(3,1);
                S(2,3) = -S(3,2);
                
                CONM.M(4:6, 1:3, n) = S;
                CONM.M(1:3, 4:6, n) = S';
                % get Inertia in a reference frame with origin in the NODE
                % and parallel to the basic reference frame
                CONM.M(4:6, 4:6, n) = CONM2.M(4:6,4:6,i) - CONM2.M(1,1,i) .* ...
                    (crossm(CONM2.Offset(i, 1:3)) * crossm(CONM2.Offset(i, 1:3)));
            end
            % constrain matrix symmetry
            %    CONM.M(:, :, n) = (CONM.M(:, :, n) + CONM.M(:, :, n)') ./ 2.0;
            
        end
        
    end
    INFO.nconm = mtot;
 end