%% VLM mesh deformation script for the SPLINE2 card
% The following file takes the undeformed (NODE), deformed (UPD_NODEPOS)
% node coordinates, node rotations (UPD_NODER) and the AERO object to
% deform the mesh, recalculate the panel unit normals and if requested,
% rotation of the unit normal to account for twist. The file has been
% modified from the original NeoCASS file.

function lattice_defo = update_vlm_mesh(NODE, UPD_NODEPOS, UPD_NODER, AERO, beam_model)

lattice_defo = AERO.lattice_vlm;

try
    AERO.Interp.Ic;
catch
    error('Collocation points interface matrix not available.');
end

try
    AERO.Interp.In;
catch
    error('Panel nodes interface matrix not available.');
end

try
    AERO.Interp.Iv;
catch
    error('Vortex points interface matrix not available.');
end

try
    AERO.Interp.Imv;
catch
    error('Midpoint vortex points interface matrix not available.');
end

% Initialise some parameters
ngrid = size(NODE.Coord, 1);
np = size(AERO.Interp.Ic, 1);
COLLOC_C = zeros(np, 3);
N = zeros(np, 3);
NODE_C = zeros(np, 5, 3);
VORTEX_C = zeros(np, 8, 3);

% Find the delta coordinates for the mesh (relative to the undeformed mesh)
[COLLOC_D, NODE_D, VORTEX_D] = deform_vlm_mesh(NODE, UPD_NODEPOS, UPD_NODER, AERO.Interp.Ic, AERO.Interp.In, AERO.Interp.Iv);

% new colloc position
COLLOC_C = AERO.lattice_vlm.COLLOC + COLLOC_D;
% new node position
NODE_C = AERO.lattice_vlm.XYZ + NODE_D;
% new vortex position
VORTEX_C = AERO.lattice_vlm.VORTEX + VORTEX_D;

% VORTEX_C(1,[1,2,3,4],1) = AERO.lattice_vlm.VORTEX(1,[1,2,3,4],1);
% VORTEX_C(1,[1,2,3,4],2) = AERO.lattice_vlm.VORTEX(1,[1,2,3,4],2);
% VORTEX_C(1,[1,2,3,4],3) = AERO.lattice_vlm.VORTEX(1,[1,2,3,4],3);
% VORTEX_C(51,[1,2,3,4],1) = AERO.lattice_vlm.VORTEX(51,[1,2,3,4],1);
% VORTEX_C(51,[1,2,3,4],2) = AERO.lattice_vlm.VORTEX(51,[1,2,3,4],2);
% VORTEX_C(51,[1,2,3,4],3) = AERO.lattice_vlm.VORTEX(51,[1,2,3,4],3);

% Determine the new flat plate normal
N = zeros(size(COLLOC_C));
for nID = 1 : length(AERO.ID)
    nPT = find(AERO.lattice_vlm.DOF(nID,:,1));
    for nn = nPT
        N(AERO.lattice_vlm.DOF(nID,nn,1):AERO.lattice_vlm.DOF(nID,nn,2),:) = get_defo_normal(COLLOC_C(AERO.lattice_vlm.DOF(nID,nn,1):AERO.lattice_vlm.DOF(nID,nn,2),:),...
            VORTEX_C(AERO.lattice_vlm.DOF(nID,nn,1):AERO.lattice_vlm.DOF(nID,nn,2),:,:),...
            AERO.lattice_vlm.DN(AERO.lattice_vlm.DOF(nID,nn,1):AERO.lattice_vlm.DOF(nID,nn,2),:),AERO.geo.b(nID));
        % There exists a get_defo_normal_diag script as well!
    end
end

lattice_defo.NOrig = N;

% The following routine has been added to rotate the lattice normal vector
% to account for twist without physically rotating the mesh.
%beam_model.camber_corr = 0;
if beam_model.twist_corr == 1 || beam_model.camber_corr == 1
    % Find the master grid nodes within the SET object for the mesh
    % assuming that all the beam nodes are grouped into one set
    MasterNode = cell(numel(AERO.Set.Node),1);
    for i = 1:numel(AERO.Set.Node)
        for j = AERO.Set.Node(i).data' % This calls all the grid nodes in that set
            if ~isempty(NODE.Aero.Index(j).data)
                MasterNode{i,1} = [MasterNode{i,1} j]; % This list the masternodes for that set
            end
        end
    end
    
    %AERO.Interp.Set ---> Use to isolate sets
    
    % Isolate the leading collocation points for interpolation of the wing
    % local system
    %leading_colloc = [];
    %SI = unique(AERO.Interp.Set);
    yloc_axis = zeros(sum(AERO.geo.ny),3);
    for i  = 1:numel(AERO.geo.ny)
        
        nPT = AERO.lattice_vlm.DOF(i,:,1);
        
        for j = 1:AERO.geo.ny(i)
            
            n1(1:3) = VORTEX_C(nPT-1+j,3,:);
            n2(1:3) = VORTEX_C(nPT-1+j,4,:);
            n3(1:3) = VORTEX_C(nPT-1+j,5,:);
            n4(1:3) = VORTEX_C(nPT-1+j,6,:);
            
            mid_x = 0.5*(n1-n2 + n4-n3);
            
            mid_x = mid_x./norm(mid_x);
            
            yloc_axis(i,:) = cross(mid_x,N(nPT-1+j,:));
            yloc_axis(i,:) = yloc_axis(i,:)./norm(yloc_axis(i,:));
        end
        
        % Take us to the end of the trailing edge panel
    end
    
%     % Count the number of leading edge panels in that set   
%     leading_count = num2cell(zeros(size(SI)));
%     for i = 1:sum(AERO.geo.ny)
%         leading_count{AERO.Interp.Set(i)} = leading_count{AERO.Interp.Set(i)} + 1;
%         for j = 1:AERO.geo.nx(i) + AERO.geo.fnx(i)
%             panel_count = panel_count + 1;
%         end
%         leading_colloc{AERO.Interp.Set(i)}(leading_count{AERO.Interp.Set(i)},:) = COLLOC_C(panel_count,:);
%         
%         Dihedral(AERO.Interp.Set(i)) = AERO.geo.dihed(i);
%     end
%     
%     % Interpolate the grid associated coordinate system to the mesh collocation
%     % points
%     yloc_axis = [];
%     for i = SI
%         
%         if Dihedral(i) == 0.5*pi
%             eyx = interp1(NODE.Coord(MasterNode{i,1},3),squeeze(UPD_NODER(1,3,MasterNode{i,1})),leading_colloc{i}(:,3));
%             eyy = interp1(NODE.Coord(MasterNode{i,1},3),squeeze(UPD_NODER(2,3,MasterNode{i,1})),leading_colloc{i}(:,3));
%             eyz = interp1(NODE.Coord(MasterNode{i,1},3),squeeze(UPD_NODER(3,3,MasterNode{i,1})),leading_colloc{i}(:,3)); 
%         else
%             
%             eyx = interp1(NODE.Coord(MasterNode{i,1},2),squeeze(UPD_NODER(1,2,MasterNode{i,1})),leading_colloc{i}(:,2));
%             eyy = interp1(NODE.Coord(MasterNode{i,1},2),squeeze(UPD_NODER(2,2,MasterNode{i,1})),leading_colloc{i}(:,2));
%             eyz = interp1(NODE.Coord(MasterNode{i,1},2),squeeze(UPD_NODER(3,2,MasterNode{i,1})),leading_colloc{i}(:,2));
%         end
%         % Store as vector and normalise for guarantee of unit
%         yloc_axis = [yloc_axis; [eyx, eyy, eyz]./repmat(sqrt(sum([eyx, eyy, eyz].^2,2)),1,3)];
%         % TODO check CPU time of repmat
%         
%         % Routine to check the interpolation
%         InterpCheck = 0;
%         if InterpCheck == 1
%             figure;
%             plot(NODE.Coord(MasterNode{i,1},2),squeeze(UPD_NODER(1,2,MasterNode{i,1})),'k-');
%             hold on;
%             plot(leading_colloc{i}(:,2),eyx,'r-');
%             figure;
%             plot(NODE.Coord(MasterNode{i,1},2),squeeze(UPD_NODER(2,2,MasterNode{i,1})),'k-');
%             hold on;
%             plot(leading_colloc{i}(:,2),eyy,'r-');
%             figure;
%             plot(NODE.Coord(MasterNode{i,1},2),squeeze(UPD_NODER(3,2,MasterNode{i,1})),'k-');
%             hold on;
%             plot(leading_colloc{i}(:,2),eyz,'r-');
%         end
%     end
    
    % Find the mean twist of a panel
    twist = mean([squeeze(AERO.geo.TW(:,:,1)),squeeze(AERO.geo.TW(:,:,2))],2);
    %twist = AERO.geo.TWIST;
    count = 0;
    
    % Rotation matrix for twist about the yloc_axis vector
    for i = 1:numel(AERO.geo.ny)
        
        for j = 1:AERO.geo.ny(i)
            % Need to check and alter for full wing
            alpha = twist(i);
            
            % Call the vector for that panel
            X = yloc_axis(i,1);
            Y = yloc_axis(i,2);
            Z = yloc_axis(i,3);
            
            % Generate twist rotation matrix
            Twist_mat = RotVecHinge(X,Y,Z,alpha);
            
            % I should be saving the camber angle here so that I can
            % recalculate my rotation matrix to account for the camber
            
            for k = 1:AERO.geo.nx(i) + AERO.geo.fnx(i)
                count = count + 1;
                
                % Rotate the panel unit normal
                N(count,:) = (Twist_mat*N(count,:)')';
                
                % Generate camber rotation matrix
                if beam_model.camber_corr == 1
                    camber = AERO.lattice_vlm.Camber(count);
                    Camber_mat = RotVecHinge(X,Y,Z,-camber);
                    N(count,:) = (Camber_mat* N(count,:)')';
                else
                    Camber_mat = eye(3,3);
                end
                
                lattice_defo.TRM(:,:,count) = Twist_mat;
                lattice_defo.CRM(:,:,count) = Camber_mat;
            end
        end
    end
end

lattice_defo.N = N;

% Update the lattice
lattice_defo.COLLOC = COLLOC_C;
lattice_defo.VORTEX = VORTEX_C;
lattice_defo.XYZ = NODE_C;

% Determine the new wake shape
lattice_defo = defo_wakesetup(lattice_defo, AERO.state, AERO.ref);

end
%***********************************************************************************************************************
