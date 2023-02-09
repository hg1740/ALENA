function [Matrices,Aero,MassProps] = InitialiseRobbie(aircraft_data ,extra_flag) %,hp)

if nargin < 2 || isempty(extra_flag)
    extra_flag = 0;
end

%% Inputs

%Number of points describing the aerofoil upper/lower surface for rendering
nProfilePoints = 10;

%% Set up part data

%Grab unique part names
partNames  = {aircraft_data.PartId.Part};
allPartIDs = unique(partNames, 'stable');
nPart      = numel(allPartIDs);

%Find the wing
WingIdx = find(ismember(allPartIDs, 'PortWing'));

%Grab the node and aero ID numbers for each part
[node_IDs, aero_IDs] = grabNodeAndAeroIDs(nPart, allPartIDs, partNames, aircraft_data);

fuse_nodes    = [];

flap_length_pc_chord = [0.25;0.25;0.25;0.25;0.25];

%% Set up temp variables to allow easier indexing
nodeCoords = aircraft_data.Node.Coord;
nodeIDs    = aircraft_data.Node.ID;
aeroIDs    = aircraft_data.Aero.ID;

%% Function defintion(s)

    function [node_IDs, aero_IDs] = grabNodeAndAeroIDs(nPart, allPartIDs, partNames, aircraft_data)
        %grabNodeAndAeroIDs Returns a cell array of ID numbers for all NODE
        %and CAERO1 data in the aircraft
        %
        %   - TODO : This should map across to the FEModel class quite well.
        
        node_IDs = cell(1, nPart);
        aero_IDs = cell(1, nPart);
        
        for iPart = 1 : nPart
            
            %Grab associated parts
            %   - In Dario's structure a 'part' is just a definition of ID numbers
            %   for a set of GRID, CBAR or CAERO1 cards that are associated with
            %   the ALENA objects by using the same name, i.e. 'StbdWing'.
            part_idx     = ismember(partNames, allPartIDs{iPart});
            part         = aircraft_data.PartId(part_idx);
            partType     = {part.Type};
            
            %Structural data
            part_nodeidx = ismember(partType, 'GRID');
            part_node    = part(part_nodeidx);
            node_IDs{iPart} = part_node.data;
            
            %Aero data
            part_aeroidx = ismember(partType, 'CAERO');
            part_aero    = part(part_aeroidx);
            if isempty(part_aero)
                aero_IDs{iPart} = 0;
            else
                aero_IDs{iPart} = part_aero.data;
            end
            
        end
        
    end

    function [beamLength, feBeamLength] = getBeamPositions(nElem, partNodes, nodeIDs, nodeCoords)
        %getBeamPositions Calculates the straight line distance for each beam
        %element as well as the total length of the geometry beam.
        
        %   - Calculate the straight line distance between nodes
        s_coord   = zeros(3, nElem);
        for iElem = 1 : nElem
            idxNodeB = nodeIDs == partNodes(iElem + 1);
            idxNodeA = nodeIDs == partNodes(iElem);
            s_coord(:, iElem) = nodeCoords(idxNodeB, :) - nodeCoords(idxNodeA, :);
        end
        
        %Overall length of this beam and the beam element coords
        %   - TODO : Get this from the beam model!
        hyp = sqrt(sum(s_coord.^2, 1));
        beamLength   = sum(hyp);
        feBeamLength = cumsum([0 hyp]);
        
    end

    function interpCoords = s2xyz(nodeCoords, nodeS, sSample)
        %interpCoords Interpolates the s-coordinates (nodeS) of
        %'nodeCoords' at the sample point 'sSample' to return the 3 x N
        %coordinates in the (x,y,z) domain.
  
       interpCoords = i_interp1(nodeS, nodeCoords, sSample);
        
    end

    function hG = plotCoordSys(hAx, O, RMatrix, clr, nam)
        %plotCoordSys Plots the rotation matrix 'rMatrix' at the (x,y,z)
        %position 'r0' in the axes object 'hAx'.
        
        %Set up origin data
        O = O';
        Ox   = O(1, :); Oy = O(2, :); Oz = O(3, :);
        nan_ = nan(size(Ox));
        
        %Extract direction vectors
        eX = squeeze(RMatrix(:, 1, :)) + O;
        eY = squeeze(RMatrix(:, 2, :)) + O;
        eZ = squeeze(RMatrix(:, 3, :)) + O;        
        
        %Set up vector data - pad with nan
        OXx = [Ox ; eX(1, :) ; nan_];
        OXy = [Oy ; eX(2, :) ; nan_];
        OXz = [Oz ; eX(3, :) ; nan_];
        OYx = [Ox ; eY(1, :) ; nan_];
        OYy = [Oy ; eY(2, :) ; nan_];
        OYz = [Oz ; eY(3, :) ; nan_];
        OZx = [Ox ; eZ(1, :) ; nan_];
        OZy = [Oy ; eZ(2, :) ; nan_];
        OZz = [Oz ; eZ(3, :) ; nan_];
        
        %Column format => vectorised plotting        
        xData = [OXx(:) ; OYx(:) ; OZx(:)];
        yData = [OXy(:) ; OYy(:) ; OZy(:)];
        zData = [OXz(:) ; OYz(:) ; OZz(:)];
        
        hG = line('Parent', hAx  , ...
            'XData'      , xData , ...
            'YData'      , yData , ...
            'ZData'      , zData , ...
            'Marker'     , 'none', ...
            'LineStyle'  , '-'   , ...
            'LineWidth'  , 1.5   , ...
            'Color'      , clr   , ...
            'DisplayName', nam);    
        
        %Add text
        text(eX(1, :), eX(2, :), eX(3, :), 'X');
        text(eY(1, :), eY(2, :), eY(3, :), 'Y');
        text(eZ(1, :), eZ(2, :), eZ(3, :), 'Z');
        
    end
    
%% Preallocate
start_node = cellfun(@(x) x(1), node_IDs);
n_node     = cellfun(@numel, node_IDs);
n_elem     = n_node - 1;
n_points   = zeros(1, nPart);

if (extra_flag == 1) && (n_elem > 1)
    fprintf(['\nThis method is not supported for aircraft with more '  , ...
        'than one member. Only using the first defined member for the ', ...
        'analysis...\n'])
    n_elem = n_elem(1);
end

%Preallocate
%   - Generic
cell_1_nP       = cell(1, nPart);
cell_nP_1       = cell(nPart, 1);
cell_zrs_1_nE   = arrayfun(@(nE) zeros(1, nE)   , n_elem, 'Unif', false);
cell_zrs_3_nE   = arrayfun(@(nE) zeros(3, nE)   , n_elem, 'Unif', false);
cell_zrs_3_3_nE = arrayfun(@(nE) zeros(3, 3, nE), n_elem, 'Unif', false);
cell_zrs_6_6_nE = arrayfun(@(nE) zeros(6, 6, nE), n_elem, 'Unif', false);
%   - Geometry
l    = zeros(1, nPart); %Length of the geometry beam
s    = cell_1_nP;       %Straight line distance along the FE beams
type = cell_1_nP;
local_dihedral = cell_1_nP;
local_twist    = cell_1_nP;
local_sweep    = cell_1_nP;
%   - Structural & aero mesh
s2            = cell_1_nP;
s_mp          = cell_1_nP;
s_aero        = cell_1_nP;
s_aero_mp     = cell_1_nP;
s_cust        = cell_1_nP;
s_out         = cell_1_nP;
s_out2        = cell_1_nP;
s_node_ind    = cell_1_nP;
s_node_ind2   = cell_1_nP;
s_mp_ind      = cell_1_nP;
s_aero_ind    = cell_1_nP;
s_aero_mp_ind = cell_1_nP;

p  = cell_1_nP;
p0 = cell_1_nP;
%   - Rotation matrices
CaB0   = cell_1_nP;
CaBi   = cell_zrs_3_3_nE;
CBB    = cell_zrs_3_3_nE;
CBA0   = cell_zrs_3_3_nE;
CBB_tw = cell_zrs_3_3_nE;

K      = cell_zrs_6_6_nE;
CFFbar = cell_zrs_6_6_nE;
C      = cell_zrs_6_6_nE;
%   - Mass matrix
M_pt_store = cell_zrs_6_6_nE;
%   - Aero terms
aerotwist = cell_zrs_1_nE;
chord     = cell_zrs_1_nE;
dCLdd     = cell_zrs_1_nE;
dCMdd     = cell_zrs_1_nE;
le        = cell_zrs_3_nE;
cp        = cell_zrs_3_nE;
te        = cell_zrs_3_nE;
tq        = cell_zrs_3_nE;
%   - Aero settings
delta_flap = arrayfun(@(nE) zeros(nE, 1), n_elem, 'Unif', false);
tau        = ones(1, nPart);    
lift_dist  = arrayfun(@(nN) ones(1, nN) , n_node, 'Unif', false); 
%   - Rendering
nP =  2 * (nProfilePoints + 1);
RenderPoints = arrayfun(@(nE) struct( ...
    'x_upper', zeros(nE, nP), ...
    'x_lower', zeros(nE, nP), ...
    'y_upper', zeros(nE, nP), ...
    'y_lower', zeros(nE, nP), ...
    'z_upper', zeros(nE, nP), ...
    'z_lower', zeros(nE, nP)), n_elem, 'Unif', false);
%   - Initial conditions
x        = arrayfun(@(nE) zeros(nE * 6, 1), n_elem, 'Unif', false); %Initial values for strains & curvatures in each element
%   - Forces
f        = arrayfun(@(nN) zeros(nN * 6, 1), n_node, 'Unif', false); %Follower force (f_b??) (at grid points - hence nElem + 1)
f_a      = arrayfun(@(nN) zeros(nN * 6, 1), n_node, 'Unif', false); %Non-follower force in 'aircraft frame' (at grid points)
f_mp_pnt = arrayfun(@(nE) zeros(nE * 6, 1), n_elem, 'Unif', false); %Point force applied to middle of element

%% Define the geometry, stiffness matrix, mass/inertia matrix and 
%% aerodynamic mesh for each part in turn

% %Generate a method object for each part
% MethodObj = arrayfun(@(~) awi.methods.IntrinsicStrainFE, 1 : nPart);
% 
% for jj = 1 : nPart
%     
%     
%     
%     
%     
%     
%     
% end


%TODO - Convert from cell notation to structure arrays.

for jj = 1 : nPart
    
    partName = lower(allPartIDs{jj});
    
    fprintf(['Setting up element properties for Member ',num2str(jj),' (''', partName, ''')...\n'])
    
    %Grab node & element data
    partNodes = node_IDs{jj};
    partAero  = aero_IDs{jj};
    nNode = n_node(jj);
    nElem = n_elem(jj);
     
    %Coordinates of this part 
    partCoords = nodeCoords(any(nodeIDs == partNodes, 2), :);
    
    % Beam Lengths
    [l(jj), s{jj}] = getBeamPositions(nElem, partNodes, nodeIDs, nodeCoords);
    
    % Surface Type
    type{jj}   = 'LS';
    
    % Dihedral, Twist, Sweep --> Unused??
    local_dihedral{jj} =  zeros(1, nNode);
    local_twist{jj}    =  zeros(1, nNode);
    local_sweep{jj}    =  zeros(1, nNode);
    
    %Get strip geometry
    indAero   = arrayfun(@(aID) find(aeroIDs == aID), partAero, 'Unif', false);
    indAero(cellfun(@isempty, indAero)) = {false};
    indAero = horzcat(indAero{:});
    chord_vec = aircraft_data.Aero.geo.c(indAero);
    %   - (x,y,z) position of LE 
    xAero     = aircraft_data.Aero.geo.startx(indAero);
    yAero     = aircraft_data.Aero.geo.starty(indAero);
    zAero     = aircraft_data.Aero.geo.startz(indAero);
    tw_vec    = aircraft_data.Aero.geo.TW(indAero, :, 1);
    
    %% Beam s-coordinates
    
    %Structural points
    %   - Offset
    s2{jj}     = sort([s{jj}, s{jj}(2:end-1)-1e-5]);
    %   - Mid-point
    s_mp{jj}   = (s{jj}(1 : end - 1) + s{jj}(2 : end)) / 2;%/4

    %Aero points
    dsAero     = l(jj) / (nElem * 2);   %2 strips per element?
    s_aero{jj} = 0 : dsAero : l(jj);%*2*4/2/4
    %   - Midpoint
    s_aero_mp{jj} = (s_aero{jj}(1 : end - 1) + s_aero{jj}(2 : end)) / 2;%(s{1}(1:end-1)+s{1}(2:end))/2;
    s_cust{jj} = [];%14.3;
    
    %Unique points
    s_out{jj}    = unique([s{jj}, s_mp{jj}, s_cust{jj}, s_aero{jj}, s_aero_mp{jj}]);
    n_points(jj) = numel(s_out{jj});
    
    %   - Offset
    s_out2{jj} = sort([s_out{jj}, s{jj}(2 : end - 1) - 1e-5]);
    
    %Indexing...not sure why at the moment
    s_node_ind{jj}    = find(any((s{jj}' - s_out{jj}) == 0, 1));
    s_node_ind2{jj}   = find(any((s2{jj}' - s_out2{jj}) == 0, 1));
    s_mp_ind{jj}      = find(any((s_mp{jj}' - s_out{jj}) == 0, 1));
    s_aero_ind{jj}    = find(any((s_aero{jj}' - s_out{jj}) == 0, 1));
    s_aero_mp_ind{jj} = find(any((s_aero_mp{jj}' - s_out{jj}) == 0, 1));
    
    p{jj}              = zeros(nNode*3, 1);
    p{jj}(1 : 3 : end) = s{jj};
    
    %Offset distance from the ref node (here this is the wing root node)
    %The ref node has to be at the start of a beam due to boundary
    %conditions.
    p0{jj} = nodeCoords(nodeIDs == start_node(jj), :)' - nodeCoords(nodeIDs == node_IDs{WingIdx}(1),:)';  %Reference node??
        
    %% Element set up (structural, mass/inertia & aero)
    
    %Calculate rotation matrix at the first element of this part    
    %   - eX is the vector between nodes (along the beam)
    ex = nodeCoords(nodeIDs == partNodes(2), :)' - nodeCoords(nodeIDs == partNodes(1), :)';
    %   - Set initial eY as the orientation from the Nastran-esque beam 
    idx = and(aircraft_data.Beam.Conn(:, 1) == partNodes(1), aircraft_data.Beam.Conn(:, 2) == partNodes(2));
    ey_tmp = aircraft_data.Beam.Orient(idx, :)'; %Global CS? Must be...
    ez = cross(ex,ey_tmp);
    ey = cross(ez,ex);
    CaB0{jj} = [ex / norm(ex), ey / norm(ey), ez / norm(ez)];

    for ii = 1 : nElem
        
        %Indexing
        idxA = nodeIDs == partNodes(ii);
        idxB = nodeIDs == partNodes(ii + 1);
        idxBeam     = and( ...
            aircraft_data.Beam.Conn(:, 1) == partNodes(ii), ...
            aircraft_data.Beam.Conn(:, 2) == partNodes(ii + 1));
        idxBeamProp = ...
            aircraft_data.PBeam.ID == aircraft_data.Beam.PID(idxBeam);
        
        %Node coordinates
        coordsA = nodeCoords(idxA, :);
        coordsB = nodeCoords(idxB, :);
        
        %% Structural Setup        
                
        %Calculate rotation matrix for each element
        %   - TODO --> Can be vectorised!!
        ex_beam = coordsB' - coordsA';
        ey_tmp  = aircraft_data.Beam.Orient(idxBeam,:)';
        ez      = cross(ex_beam, ey_tmp);
        ey      = cross(ez, ex_beam);
        CaBi{jj}(:,:,ii) = [ex_beam/norm(ex_beam), ey/norm(ey), ez/norm(ez)];
        
        %Transformation from one element to the next
        if ii == 1
            CBB{jj}(:,:,ii) = (CaB0{jj}'           * CaBi{jj}(:,:,ii))';
        else
            CBB{jj}(:,:,ii) = (CaBi{jj}(:,:,ii-1)' * CaBi{jj}(:,:,ii))';
        end
        
        %Set up stiffness/flexibility matrices...
        
        %Set up the beam sectional & material properties
        %   - Material
        if isfield(aircraft_data.PBeam, 'Mat')
            E   = aircraft_data.Mat.E(aircraft_data.PBeam.Mat(idxBeamProp));
            G   = aircraft_data.Mat.G(aircraft_data.PBeam.Mat(idxBeamProp));
            rho = aircraft_data.Mat.Rho(aircraft_data.PBeam.Mat(idxBeamProp));
        else
            idxMat = ismember(aircraft_data.Mat.ID, aircraft_data.PBeam.MID(ii));
            E     = aircraft_data.Mat.E(idxMat);
            G     = aircraft_data.Mat.G(idxMat);
            rho   = aircraft_data.Mat.Rho(idxMat);
        end        
        %   - Sectional
        A        = mean(aircraft_data.PBeam.A(idxBeamProp).data);
        I1       = mean(aircraft_data.PBeam.I(idxBeamProp).data(2, :));
        I2       = mean(aircraft_data.PBeam.I(idxBeamProp).data(1, :));
        I12      = mean(aircraft_data.PBeam.I(idxBeamProp).data(3, :));
        J_       = mean(aircraft_data.PBeam.J(idxBeamProp).data);
        offset_a = aircraft_data.Beam.Offset(idxBeamProp, 1 : 3);
        offset_b = aircraft_data.Beam.Offset(idxBeamProp, 7 : 9);
        R_shear  = (offset_a' + offset_b') / 2; %average of offsets at end B & end A
        
        %Here we determine the rotation matrix between the shear axis and
        %the element axis
        %   - N.B. As any rotation matrix (R) must be orthogonal it has the
        %   property transpose(R) = inv(R).        
        ex       = (offset_b + coordsB - coordsA - offset_a)';
        ey_tmp   = aircraft_data.Beam.Orient(idxBeam, :)';
        ez       = cross(ex, ey_tmp);
        ey       = cross(ez, ex);
        %   - Rotation matrix from element axis to shear axis
        %   R_shear = CBBbar * R_element
        CBBbar   = CaBi{jj}(:,:,ii)' * [ex/norm(ex), ey/norm(ey), ez/norm(ez)];
        R_tilde  = [eye(3) -skew(CaBi{jj}(:,:,ii)' * R_shear) ; zeros(3) eye(3)]; %R_tilde accounts for axial load induced by shear centre offset when beam element bends (cross product of shear centre moment arm)
        dsprime  = norm(ex);
        ds       = s{jj}(ii + 1) - s{jj}(ii);
        length_fact = ds/dsprime;
        
        %Set up stiffness matrix
        K{jj}(:,:,ii)  = eye(6);
        K{jj}(1,1,ii)  = E*A;
        K{jj}(2,2,ii)  = G*A; %shear terms (Timoshenko?)
        K{jj}(3,3,ii)  = G*A; %shear terms
        K{jj}(4,4,ii)  = G*J_;
        K{jj}(5,5,ii)  = E*I1;
        K{jj}(6,6,ii)  = E*I2;
        K{jj}(5,6,ii)  = E*I12; %TODO - When working with the HARTEN model we need to multiply this by 0.
        K{jj}(6,5,ii)  = E*I12; %TODO - When working with the HARTEN model we need to multiply this by 0.
        %Pre and post multiply stiffness matrix by rotation matrices in
        %order to evaluate loads about the beam axis and not the shear axis
        %(rotate and shift loads)
        K{jj}(:,:,ii)  = R_tilde'*blkdiag(CBBbar,CBBbar)*K{jj}(:,:,ii)*blkdiag(CBBbar',CBBbar')*R_tilde*length_fact; 
        
        %Rotation matrix to transform loads from beam axis to shear axis 
        CFFbar{jj}(:,:,ii) = inv(R_tilde'*blkdiag(CBBbar,CBBbar));
        
        %Compliance matrix
        C{jj}(:,:,ii)  = inv(K{jj}(:,:,ii));
        
        %% Set up mass & inertia properties...
        %   - Need to transform the point mass definition into a smeared
        %   mass and inertia. 
        %   - TODO : Lots of work here to understand!!!
        
        %Something to do with the radius of gyration??
        h = sqrt(sqrt(I1/I2)*A);
        b = A/h;
                       
        M_pt       = zeros(6);
        
        % Search for lumped masses associated with the end element
        CONM2_list = find(aircraft_data.Conm2.Node==node_IDs{jj}(ii+1));
        for mm = 1:length(CONM2_list)
            M_pt = M_pt + blkdiag(CaBi{jj}(:,:,ii)',CaBi{jj}(:,:,ii)')*[eye(3) zeros(3);skew(aircraft_data.Conm2.Offset(CONM2_list(mm),:))  eye(3)]*aircraft_data.Conm2.M(:,:,CONM2_list(mm))*[eye(3) -skew(aircraft_data.Conm2.Offset(CONM2_list(mm),:));zeros(3)  eye(3)]*blkdiag(CaBi{jj}(:,:,ii),CaBi{jj}(:,:,ii));
        end
        
        M_pt_store{jj}(:,:,ii+1) = zeros(6);%M_pt;
        if ii == 1
            M_pt_store{jj}(:,:,1) = zeros(6);
        end
        
        m_pt0     =    M_pt(1,1);
        if m_pt0 == 0
            M_pt0  = zeros(6);
            cg_pt0 = [0;0;0];
        else
            cg_pt     = -[-M_pt(2,6);M_pt(1,6);-M_pt(1,5)]/m_pt0;
            M_pt0     =    M_pt - [zeros(3) -skew(cg_pt)*m_pt0; skew(cg_pt)*m_pt0 -skew(cg_pt)*skew(cg_pt)*m_pt0];
            cg_pt0    =    cg_pt;
        end
        
        J_pt0 = M_pt(4:6,4:6);
        
        m{jj}(ii)      = rho*A + m_pt0/ds; %Smear point mass over the length of the element
        i1             = rho*A/12*(b^2 + h^2) ;
        i2             = rho*A*(b^2 + ds^2)/12;%1e-5;%
        i3             = rho*A*(h^2 + ds^2)/12;%1e-5;%
        J{jj}(:,:,ii)  = diag([i1,i2,i3]) + J_pt0/ds;
        cg{jj}(:,:,ii) = [0;0;0] + cg_pt0;
        
        %% Aero Setup
        
        %Position of element mid-points
        rMid_ = (nodeCoords(idxB, :) + nodeCoords(idxA, :)) / 2;
        
        % TODO - Remove this from the code fudge needed to account for
        % strip theory
        tw_vec = smooth(tw_vec, 20);
        
        switch partName
            case 'vtp'
                %VTP is different because the interpolation is performed
                %w.r.t the z coordinates. Could we just use s-coordinates
                %instead? Then it would be generalisable...
                
                %Interpolate LE position at element mid-points (global sys)
                xLE = i_interp1(zAero, xAero, rMid_(3), 'linear', 'extrap');
                yLE = i_interp1(zAero, yAero, rMid_(3), 'linear', 'extrap');
                %Interpolate data to the midpoint of the nodes...
                aerotwist{jj}(ii) = i_interp1(zAero, tw_vec   , rMid_(3), 'linear', 'extrap');
                chord{jj}(ii)     = i_interp1(zAero, chord_vec, rMid_(3), 'linear', 'extrap');
                rLE = [xLE, yLE, rMid_(3)];
                
                CBB_tw{jj}(:,:,ii) = eye(3);
                
            otherwise
                       
                %Interpolate LE position at element mid-points (global sys)
                xLE = i_interp1(yAero, xAero, rMid_(2), 'linear', 'extrap');
                zLE = i_interp1(yAero, zAero, rMid_(2), 'linear', 'extrap');
                %Interpolate twist & chord data to the nodes midpoints
                aerotwist{jj}(ii) = i_interp1(yAero, tw_vec,    rMid_(2), 'linear', 'extrap');
                chord{jj}(ii)     = i_interp1(yAero, chord_vec, rMid_(2), 'linear', 'extrap');
                
                rLE = [xLE, rMid_(2), zLE];
                
                %TODO - Robbie & Dario to fix the twist implementation.
                %Possibly an issue with sign convention between the port &
                %starboard wing.
                CBB_tw{jj}(:,:,ii) =  eye(3);%CaBi{jj}(:,:,ii)'*RotVecHinge(ex_beam(1),ex_beam(2),ex_beam(3),(pi/180)*aerotwist{jj}(ii))*CaBi{jj}(:,:,ii);

        end
        
        %   - Offset between node position and LE position in the
        %   global CS at the element mid-point.
        rLE = rLE - rMid_;
        
        %Offset of LE from mid-point node in beam CS (vector from
        %beam node to LE)
        le{jj}(:, ii) = CaBi{jj}(:,:,ii)' * rLE';
                
        %Construct aero rotation matrix (local CS)
        ez      = cross(le{jj}(:,ii), [1;0;0]);
        ex      = cross(le{jj}(:,ii), ez);
        % Remove the x contribution - * * EXTANT CODE * * (Not sure why
        % this was ever applicable but keep in for now)
        %ex_beam(1) = 0;
        %ex_beam = ex_beam/norm(ex_beam);
        CB_aero = [ex/norm(ex), le{jj}(:,ii)/norm(le{jj}(:,ii)), ez/norm(ez)];
        
        %Aerodynamic reference frame (A)
        %   - Rotate by twist angle
        CBA0{jj}(:,:,ii) = CBB_tw{jj}(:,:,ii) * CB_aero;
        %   - Calculate the position of the Centre of Pressure, TE
        %   and 3/4 chord point
        r_cp = [0.0 ; -0.25 * chord{jj}(ii); 0.0]; 
        r_te = [0.0 ; -1.00 * chord{jj}(ii); 0.0]; 
        r_tq = [0.0 ; -0.75 * chord{jj}(ii); 0.0];
        cp{jj}(:,ii)     = CBA0{jj}(:,:,ii) * r_cp + le{jj}(:, ii);
        te{jj}(:,ii)     = CBA0{jj}(:,:,ii) * r_te + le{jj}(:, ii);
        tq{jj}(:,ii)     = CBA0{jj}(:,:,ii) * r_tq + le{jj}(:, ii);
        
        %% Member flaps
        %   - Theoretical Cl/Cm for trailing edge flaps
        flap_length      = flap_length_pc_chord(jj)*chord{jj}(ii);
        k                = (chord{jj}(ii)-flap_length)/chord{jj}(ii);
        theta_k          =  acos(1-2*k);
        dCLdd{jj}(:,ii)  = (2.00*(pi-theta_k)   + 2.00*sin(theta_k));
        dCMdd{jj}(:,ii)  = (0.25*sin(2*theta_k) - 0.50*sin(theta_k));
        
        xs{1}(ii)        = 0.0;
        
    end      
    
    %% Plotting
    %   - Nodes
    %   - Elements
    %   - Coordinate systems
    hF   = figure('Name', sprintf('Memeber %i (''%s'')', jj, partName));
    hAx  = axes('Parent', hF, 'NextPlot', 'add');
    xlabel(hAx, 'X [m]'); ylabel(hAx, 'Y [m]'); zlabel(hAx, 'Z [m]');
    hLeg = legend(hAx, 'Location', 'NorthEastOutside');
    %   - Node coordinates in global frame
    hG(1) = plot3(hAx, partCoords(:, 1), partCoords(:, 2), partCoords(:, 3), ...
        'LineStyle'      , 'none', ...
        'MarkerFaceColor', 'g'   , ...
        'MarkerEdgeColor', 'k'   , ...
        'Marker'         , 'o'   , ...
        'DisplayName'    , 'Structural Nodes');
    %   - Beam elements
    hG(end + 1) = plot3(hAx, partCoords(:, 1), partCoords(:, 2), partCoords(:, 3), ...
        'LineStyle'  , '--'  , ...
        'Color'      , 'k'   , ...
        'LineWidth'  , 1.5   , ...
        'Marker'     , 'none', ...
        'DisplayName', 'Structural Elements');
    %   - Planform
    xyzLE = [xAero, yAero, zAero];
    xyzTE = xyzLE + [chord_vec, zeros(numel(xAero), 1), zeros(numel(xAero), 1)];
    xyzPlanform = [xyzLE ; flipud(xyzTE)];
    hG(end + 1) = plot3(hAx, xyzPlanform(:, 1), xyzPlanform(:, 2), xyzPlanform(:, 3), ...
        'r--', 'DisplayName', 'Aero panel inboard (x,y,z)');
    %   - Mid coordinates (where the states are evaluated)
    mpCoords    = s2xyz(partCoords, s{jj}, s_mp{jj});
    hG(end + 1) = plot3(hAx(1), mpCoords(:, 1)     , mpCoords(:, 2)     , mpCoords(:, 3)     , 'rx', ...
        'LineStyle'  , 'none', ...
        'DisplayName', 'Structural mid-points');
    %   - Edge of aero strips
    aeroCoords  = s2xyz(partCoords, s{jj}, s_aero{jj});
    hG(end + 1) = plot3(hAx(1), aeroCoords(:, 1)   , aeroCoords(:, 2)   , aeroCoords(:, 3)   , 'b^', ...
        'LineStyle'  , 'none', ...
        'DisplayName', 'Edge of Aero Strips');
    %   - Centre of aero strips
    aero_mpCoords = s2xyz(partCoords, s{jj}, s_aero_mp{jj});
    hG(end + 1)   = plot3(hAx(1), aero_mpCoords(:, 1), aero_mpCoords(:, 2), aero_mpCoords(:, 3), 'k*', ...
        'LineStyle'  , 'none', ...
        'DisplayName', 'Aero mid-points');    
    %   - Beam coordinate system at first element (CaB0)
    hG(end + 1) = plotCoordSys(hAx, partCoords(1, :), CaB0{jj}, 'k', 'Beam CS');
    %   - Element coordinate system (CaBi)
    hG(end + 1) = plotCoordSys(hAx, partCoords(2 : end, :), CaBi{jj}, 'b', 'Element CS');
    %   - Transformation matrix from one element to the next (CBB)
    %   - Tranformation matrix from element CS to shear CS (CBBbar)
    %   - CBB_tw
    %   - Aerodynamic coordinate system (CBA0) TODO Replace with matrices
    rCP = mpCoords + [-0.25 .* chord{jj}' , zeros(nElem, 1), zeros(nElem, 1)];
    rCP_ = zeros(3, nElem);
    for iE = 1 : nElem
        %Transform from aero CS to beam element CS
        rCP_B = CBA0{jj}(:, :, iE)' * cp{jj}(:, iE);
        %Transform from beam element CS to global CS
        rCP_G = CaBi{jj}(:, :, iE)' * rCP_B;
        rCP_(:, iE) = rCP_G;
    end
    hG(end + 1) = plot3(hAx, rCP(:, 1), rCP(:, 2), rCP(:, 3), ...
        'LineStyle', 'none', 'Marker', 's', 'MarkerFaceColor', 'm', ...
        'MarkerEdgeColor', 'k', 'DisplayName', 'Aero CP');
    
    %     axis(hAx, 'equal'); 
    
    %% Profile information
    
    %Define the aerofoil profile
    iaf.designation       = '0012';
    iaf.n                 = nProfilePoints;
    iaf.HalfCosineSpacing = 1;
    iaf.wantFile          = 0;
    iaf.is_finiteTE       = 0;
    
    af                    = naca4gen(iaf);
    
    %Shortcuts
    nP  = nProfilePoints + 1;
    onz = ones(1, nP);
    
    %Interpolate the panel coordinates to the structural nodes    
    idx = any(nodeIDs == partNodes, 2); 
    uCoords = zeros(nElem + 1, nProfilePoints + 1, 3);
    lCoords = zeros(nElem + 1, nProfilePoints + 1, 3);
    switch partName
        case 'vtp'
            %VTP is different because the interpolation is performed
            %w.r.t the z coordinates. Could we just use s-coordinates
            %instead? Then it would be generalisable...
            renderchord =  i_interp1(zAero, chord_vec, nodeCoords(idx, 3), 'linear', 'extrap');
            zProfile = nodeCoords(idx, 3);
            xProfile = i_interp1(zAero, xAero, zProfile, 'linear', 'extrap');
            yProfile = i_interp1(zAero, yAero, zProfile, 'linear', 'extrap');
            rLE     = [xProfile, yProfile, zProfile];
            uCoords(:, :, 3) = nodeCoords(idx, 3) * onz;
            lCoords(:, :, 3) = nodeCoords(idx, 3) * onz;
            uCoords(:, :, 2) = (renderchord + rLE(1)) * (1 - af.xU');% - norm(le0);
            lCoords(:, :, 2) = (renderchord + rLE(1)) * (1 - af.xL');% - norm(le0);
            uCoords(:, :, 1) = renderchord * af.zU(end : -1 : 1)';% + le_a(3);
            lCoords(:, :, 1) = renderchord * af.zL(end : -1 : 1)';% + le_a(3);
        otherwise
            renderchord =  i_interp1(yAero, chord_vec, nodeCoords(idx, 2), 'linear', 'extrap');
            yProfile = nodeCoords(idx, 2);
            xProfile = i_interp1(yAero, xAero    , yProfile, 'linear', 'extrap');
            zProfile = i_interp1(yAero, zAero    , yProfile, 'linear', 'extrap');
            rLE        = [xProfile, yProfile, zProfile];
            uCoords(:, :, 1) = nodeCoords(idx, 2) * onz;
            lCoords(:, :, 1) = nodeCoords(idx, 2) * onz;
            uCoords(:, :, 2) = (renderchord + rLE(:, 1)) * (1 - af.xU');% - norm(le0);
            lCoords(:, :, 2) = (renderchord + rLE(:, 1)) * (1 - af.xL');% - norm(le0);
            uCoords(:, :, 3) = (renderchord + rLE(:, 3)) * af.zU(end:-1:1)';
            lCoords(:, :, 3) = (renderchord + rLE(:, 3)) * af.zL(end:-1:1)';
    end   
    
    %Calculate coordinates at the mid point of the structural nodes
    idxA = any(nodeIDs == partNodes(1 : end - 1), 2);
    idxB = any(nodeIDs == partNodes(2 : end), 2);
    xLE = repmat((nodeCoords(idxA , 1) + nodeCoords(idxB, 1)) / 2, [1, 2 * nP]);
    yLE = repmat((nodeCoords(idxA , 2) + nodeCoords(idxB, 2)) / 2, [1, 2 * nP]);
    zLE = repmat((nodeCoords(idxA , 3) + nodeCoords(idxB, 3)) / 2, [1, 2 * nP]);    
    midCoords = cat(3, xLE, yLE, zLE);
    
    %Create array of upper and lower surfaces to plot in a single patch
    uCoords = cat(2, uCoords(2 : end, :, :), fliplr(uCoords(1 : end - 1, :, :)));
    lCoords = cat(2, lCoords(2 : end, :, :), fliplr(lCoords(1 : end - 1, :, :)));
    
    %Shift to location of structural nodes (centre of element) and stack
    uProfileCoords = permute(uCoords - midCoords, [3, 2, 1]);
    lProfileCoords = permute(lCoords - midCoords, [3, 2, 1]);
        
    %Transform by...????? (AoA, dihedral & sweep??)
    upper_coords = arrayfun(@(i) CaBi{jj}(:,:,i)' * uProfileCoords(:, :, i), 1 : nElem, 'Unif', false);
    upper_coords = cat(3, upper_coords{:});  
    upper_coords = permute(upper_coords, [3, 2, 1]);
    lower_coords = arrayfun(@(i) CaBi{jj}(:,:,i)' * lProfileCoords(:, :, i), 1 : nElem, 'Unif', false);
    lower_coords = cat(3, lower_coords{:});
    lower_coords = permute(lower_coords, [3, 2, 1]);
    
    % * * * I may have broken this bit of code...
    %   --> Plotting from this and from 'readDario' are different!
    hF  = figure('Name', partName);
    hAx = axes('Parent', hF, 'NextPlot', 'add');
    patch(hAx, upper_coords(:, :, 1)', upper_coords(:, :, 2)', upper_coords(:, :, 3)', 'white');
    patch(hAx, lower_coords(:, :, 1)', lower_coords(:, :, 2)', lower_coords(:, :, 3)', 'white');
    plot3(hAx, xAero, yAero, zAero, 'r-');
        
    RenderPoints{jj}.x_upper = upper_coords(:, :, 1);
    RenderPoints{jj}.x_lower = lower_coords(:, :, 1);
    RenderPoints{jj}.y_upper = upper_coords(:, :, 2);% - norm(le0);
    RenderPoints{jj}.y_lower = lower_coords(:, :, 2);% - norm(le0);
    RenderPoints{jj}.z_upper = upper_coords(:, :, 3);
    RenderPoints{jj}.z_lower = lower_coords(:, :, 3);
    
%     %% Members Aero Settings
%     delta_flap{jj} = zeros(nElem,1);
%     tau(jj)        = 1.0;
%     lift_dist{jj}  = ones(size(s{jj}));

    %% Forces and initial conditions
%     x{jj}     = zeros(nElem*6,1);
%     
%     f               = zeros((nElem+1)*6,1); %Follower force (f_b??) (at grid points - hence nElem + 1)
%     f_a             = zeros((nElem+1)*6,1); %Non-follower force in 'aircraft frame' (at grid points)
%     f_mp_pnt        = zeros( nElem   *6,1); %Point force applied to middle of element
%     Matrices.f{jj}        =  f;
%     Matrices.f_a{jj}      =  f_a;
%     Matrices.f_mp_pnt{jj} =  f_mp_pnt;
    
    %Can ignore these...(different formulation)
    %x2{1} = zeros(nElem*12,1);
    f2    = zeros(nElem*12,1);
    f_a2  = zeros(nElem*12,1);
    Matrices.f2{jj}   =    f2;
    Matrices.f_a2{jj} =  f_a2;
    
end

function yq = i_interp1(x, y, xq, method, str)
    
    if nargin < 4
        method = 'linear';
    end
    if nargin < 5
        str = 'extrap';
    end
    
    if numel(x) > 1
        yq = interp1(x, y, xq, method, str);
    else
        yq = repmat(y, size(xq));
    end
    
        
end

start_node = unique(start_node);

%% RB Data

% If there exists any nodes which are rigidly connected to any start nodes
% then lump them up onto the fuselage mass matrix.
%   - Fuselage can have a single point mass at the origin on the 'a' frame.

M_rb_fus = zeros(6);
for ii = 1:length(start_node)
    RB_CONM2_list    = find(aircraft_data.Conm2.Node==start_node(ii));
    M_cg_rb_fus      = aircraft_data.Node.Coord(aircraft_data.Node.ID==start_node(ii),:)-aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{WingIdx}(1),:);
    for mm = 1:length(RB_CONM2_list)
        M_rb_fus         = M_rb_fus + [eye(3) zeros(3);skew(M_cg_rb_fus)  eye(3)]*[eye(3) zeros(3);skew(aircraft_data.Conm2.Offset(RB_CONM2_list(mm),:))  eye(3)]*aircraft_data.Conm2.M(:,:,RB_CONM2_list(mm))*[eye(3) -skew(aircraft_data.Conm2.Offset(RB_CONM2_list(mm),:)); zeros(3) eye(3)]*[eye(3) -skew(M_cg_rb_fus); zeros(3) eye(3)];
    end
end
AllRB_CONM2 = [];
% Iterate through the start nodes
if isfield(aircraft_data, 'RBE2')
    for ii = 1:length(start_node)
        % Find out if any of those start nodes are member of an Rbe2
        slave_nodes   = [aircraft_data.RBE2.IDS.data];
        Rebe2_list    = find(slave_nodes==start_node(ii));
        % Find out if there are any lumped masses connected to the master node
        RB_CONM2_list = find(aircraft_data.Conm2.Node==aircraft_data.RBE2.IDM(Rebe2_list));
        % if any RB_CONM2_list are contained in the the all list remove them
        [~,~,coin] = intersect(AllRB_CONM2,RB_CONM2_list);
        % Remove any nodes in the the AllRB_CONM2 from RB_CONM2_list
        RB_CONM2_list(coin) = [];
        % Store the conm2s so that they do not get repeated
        AllRB_CONM2  = [AllRB_CONM2;RB_CONM2_list];
        % Calculate the offset of the rigidly connected node to the start node
        M_cg_rb_fus   = aircraft_data.Node.Coord(aircraft_data.Node.ID==aircraft_data.RBE2.IDM(Rebe2_list),:)-aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{WingIdx}(1),:);
        % Iterate through the conm2s and get the offsets
        for mm = 1:length(RB_CONM2_list)
            M_rb_fus         = M_rb_fus + [eye(3) zeros(3);skew(M_cg_rb_fus)  eye(3)]*[eye(3) zeros(3);skew(aircraft_data.Conm2.Offset(RB_CONM2_list(mm),:))  eye(3)]*aircraft_data.Conm2.M(:,:,RB_CONM2_list(mm))*[eye(3) -skew(aircraft_data.Conm2.Offset(RB_CONM2_list(mm),:)); zeros(3) eye(3)]*[eye(3) -skew(M_cg_rb_fus); zeros(3) eye(3)];
        end
    end
end

M_cg_full        = -[-M_rb_fus(2,6);M_rb_fus(1,6);-M_rb_fus(1,5)];%M_cg_rb_fus*M_rb_fus(1);
Mf_tot           =    M_rb_fus(1);

Matrices.M_cg_rb =  M_cg_full/Mf_tot(1);    %yes
Matrices.M_rb    =  M_rb_fus;               %yes
%
% theta_interp = [0:5:180];
% fuse_interp  = 0.5*(1-cos(theta_interp*pi/180));
% r_interp     = i_interp1(aircraft_data.BAero.geo.fs{2},aircraft_data.BAero.geo.Rs{2},fuse_interp);
% [Z,Y,X]      = cylinder(r_interp);
% X            = repmat(fuse_interp',1,size(X,2));
% X            = X*aircraft_data.BAero.geo.L(2);
%
% Z_shift      = i_interp1([0 5 35 44],[-0.1 0 0 1],X(:,1));
% Z            = Z + repmat(Z_shift,1,size(Z,2));
%
% for ii = 1:size(X,1)
%     coords_rot = [X(ii,:);Y(ii,:);Z(ii,:)];
%     X_out(ii,:) = coords_rot(1,:);
%     Y_out(ii,:) = coords_rot(2,:);
%     Z_out(ii,:) = coords_rot(3,:);
% end
% surf(X_out,Y_out,Z_out,'FaceColor','w')
% RenderPoints_RB.x = X_out - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{1}(1),1);
% RenderPoints_RB.y = Y_out - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{1}(1),2);
% RenderPoints_RB.z = Z_out - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{1}(1),3);

%% Engines
% theta_interp = [0:30:180];
% eng_interp   = 0.5*(1-cos(theta_interp*pi/180));
% r_interp     = i_interp1(aircraft_data.BAero.geo.fs{3},aircraft_data.BAero.geo.Rs{3},eng_interp);
% [Z,Y,X]      = cylinder(r_interp);
% X            = repmat(eng_interp',1,size(X,2));
% X            = X*aircraft_data.BAero.geo.L(3);
% clear X_out Y_out Z_out
% for ii = 1:size(X,1)
%     coords_rot = CaBi{1}(:,:,10)'*[X(ii,:) + aircraft_data.BAero.geo.ref_point(3,1) - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{1}(10),1);...
%                                    Y(ii,:) + aircraft_data.BAero.geo.ref_point(3,2) - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{1}(10),2);...
%                                    Z(ii,:) + aircraft_data.BAero.geo.ref_point(3,3) - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{1}(10),3)];
%     X_out(ii,:) = coords_rot(1,:);
%     Y_out(ii,:) = coords_rot(2,:);
%     Z_out(ii,:) = coords_rot(3,:);
% end
% surf(X+aircraft_data.BAero.geo.ref_point(3,1),Y+aircraft_data.BAero.geo.ref_point(3,2),Z+aircraft_data.BAero.geo.ref_point(3,3),'FaceColor','w')
% RenderPoints_RB.Extra{1}.x = X_out;
% RenderPoints_RB.Extra{1}.y = Y_out;
% RenderPoints_RB.Extra{1}.z = Z_out;
% r_interp     = i_interp1(aircraft_data.BAero.geo.fs{1},aircraft_data.BAero.geo.Rs{1},eng_interp);
% [Z,Y,X]      = cylinder(r_interp);
% X            = repmat(eng_interp',1,size(X,2));
% X            = X*aircraft_data.BAero.geo.L(1);
% clear X_out Y_out Z_out
% for ii = 1:size(X,1)
%     coords_rot = CaBi{2}(:,:,10)'*[X(ii,:) + aircraft_data.BAero.geo.ref_point(1,1) - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{2}(10),1);...
%                                    Y(ii,:) + aircraft_data.BAero.geo.ref_point(1,2) - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{2}(10),2);...
%                                    Z(ii,:) + aircraft_data.BAero.geo.ref_point(1,3) - aircraft_data.Node.Coord(aircraft_data.Node.ID==node_IDs{2}(10),3)];
%     X_out(ii,:) = coords_rot(1,:);
%     Y_out(ii,:) = coords_rot(2,:);
%     Z_out(ii,:) = coords_rot(3,:);
% end
% surf(X+aircraft_data.BAero.geo.ref_point(1,1),Y+aircraft_data.BAero.geo.ref_point(1,2),Z+aircraft_data.BAero.geo.ref_point(1,3),'FaceColor','w')
% RenderPoints_RB.Extra{2}.x = X_out;
% RenderPoints_RB.Extra{2}.y = Y_out;
% RenderPoints_RB.Extra{2}.z = Z_out;
% RenderPoints_RB.ExtraInd   = [1 10;
%                               2 10];

%% Aero Distributions
% load('\\ads\filestore\Engineering\Research\Projects\AWI\Models\Straight_Tapered\NoTwist\Cl_alpha.mat')
% load('\\ads\filestore\Engineering\Research\Projects\AWI\Models\Scitech2018Straight_Tapered\Cl_alphaAR10.mat')
% PanelArea = PanelArea';
%
% wing1CL          = [Cl_alpha( 1:22);(Cl_alpha(23:2:33 ).*PanelArea(23:2:33 )+Cl_alpha(24:2:34 ).*PanelArea(24:2:34 ))./(PanelArea(23:2:33 )+PanelArea(24:2:34 ));Cl_alpha(35:38)];
% wing2CL          = [Cl_alpha(39:60);(Cl_alpha(61:2:71 ).*PanelArea(61:2:71 )+Cl_alpha(62:2:72 ).*PanelArea(62:2:72 ))./(PanelArea(61:2:71 )+PanelArea(62:2:72 ));Cl_alpha(73:76)];
% wing3CL          = [                (Cl_alpha(77:2:83 ).*PanelArea(77:2:83 )+Cl_alpha(78:2:84 ).*PanelArea(78:2:84 ))./(PanelArea(77:2:83 )+PanelArea(78:2:84 ));Cl_alpha(85)   ];
% wing4CL          = [                (Cl_alpha(86:2:92 ).*PanelArea(86:2:92 )+Cl_alpha(87:2:93 ).*PanelArea(87:2:93 ))./(PanelArea(86:2:92 )+PanelArea(87:2:93 ));Cl_alpha(94)   ];
% wing5CL          = [                (Cl_alpha(95:2:101).*PanelArea(95:2:101)+Cl_alpha(96:2:102).*PanelArea(96:2:102))./(PanelArea(95:2:101)+PanelArea(96:2:102));Cl_alpha(103)  ];
%
% rootCL1          = i_interp1(s_mp{1},wing1CL/2/pi,s{1}(1),'linear','extrap');
% rootCL2          = i_interp1(s_mp{2},wing2CL/2/pi,s{2}(1),'linear','extrap');
% rootCL3          = i_interp1(s_mp{3},wing3CL/2/pi,s{3}(1),'linear','extrap');
% rootCL4          = i_interp1(s_mp{4},wing4CL/2/pi,s{4}(1),'linear','extrap');
% rootCL5          = i_interp1(s_mp{5},wing5CL/2/pi,s{5}(1),'linear','extrap');
%
% n_wing           = length(wing1CL);
% Awing            = zeros(n_wing+2,n_wing+1);
% Bwing            = zeros(n_wing+1,n_wing);
% Cwing            = zeros(n_wing+1, 1);
%
% Awing(1:n_wing+1,1:n_wing+1) =                                eye(n_wing+1);
% Awing(2:n_wing+2,1:n_wing+1) = Awing(2:n_wing+2,1:n_wing+1) + eye(n_wing+1);
% Awing                        = Awing(1:n_wing+1,1:n_wing+1)                ;
% Bwing(2:n_wing+1,1:n_wing)   =                              2*eye(n_wing)  ;
% Cwing(1)                     =                              1              ;
%
% lift_dist{1}     = Awing\(Bwing*wing1CL/2/pi + Cwing*rootCL1);
% lift_dist{2}     = Awing\(Bwing*wing2CL/2/pi + Cwing*rootCL2);
%
% n_htp            = length(wing3CL);
% Ahtp             = zeros(n_htp+2,n_htp+1);
% Bhtp             = zeros(n_htp+1,n_htp);
% Chtp             = zeros(n_htp+1,1);
%
% Ahtp(1:n_htp+1,1:n_htp+1) =                             eye(n_htp+1);
% Ahtp(2:n_htp+2,1:n_htp+1) = Ahtp(2:n_htp+2,1:n_htp+1) + eye(n_htp+1);
% Ahtp                      = Ahtp(1:n_htp+1,1:n_htp+1)               ;
% Bhtp(2:n_htp+1,1:n_htp)   =                           2*eye(n_htp)  ;
% Chtp(1)                   =                           1             ;
%
% lift_dist{3}     = Ahtp\(Bhtp*wing3CL/2/pi + Chtp*rootCL3);
% lift_dist{4}     = Ahtp\(Bhtp*wing4CL/2/pi + Chtp*rootCL4);
% lift_dist{5}     = Ahtp\(Bhtp*wing5CL/2/pi + Chtp*rootCL5);

%% Preallocate system matrices

Dtot2   = cell_nP_1;
Btot2   = cell_nP_1;
e1tot2  = cell_nP_1;
Mtot2   = cell_nP_1;
Ftot2   = cell_nP_1;

R3_bc   = cell_nP_1;
R4_bc   = cell_nP_1;

loadInt2 = cell_nP_1;

Cfull         = cell_nP_1;
CFFbar_full   = cell_nP_1;

Ab_sep  = cell_nP_1;
Ab2_sep = cell_nP_1;
Aa_sep  = cell_nP_1;
Aa2_sep = cell_nP_1;

Ab1var_sep = cell_nP_1;
Ab2var_sep = cell_nP_1;
Ab3var_sep = cell_nP_1;
Ab4var_sep = cell_nP_1;
Aa1var_sep = cell_nP_1;
Aa2var_sep = cell_nP_1;
Aa3var_sep = cell_nP_1;
Aa4var_sep = cell_nP_1;
Aa5var_sep = cell_nP_1;
Aa6var_sep = cell_nP_1;

%Preallocate the zero matrices for each part
cell_6nN_6nE   = arrayfun(@(i)  zeros(n_node(i) * 6, n_elem(i) * 6)  , 1 : nPart, 'Unif', false);   %[nNode * 6 , nElem * 6]
cell_6nN_6nN   = arrayfun(@(nN) zeros(6 * nN, 6 * nN)                , n_node   , 'Unif', false);   %[nNode * 6 , nNode * 6]
cell_6nN_6nP   = arrayfun(@(i)  zeros(n_node(i) * 6, n_points(i) * 6), 1 : nPart, 'Unif', false);   %[nNode * 6 , nPoint * 6]
cell_6nE_1     = arrayfun(@(nE) zeros(nE * 6, 1)                     , n_elem   , 'Unif', false);   %[nElem * 6 , 1]
cell_12nN_12nN = arrayfun(@(nN) zeros(nN * 12, nN * 12)              , n_node   , 'Unif', false);   %[nNode * 12, nNode * 12]
cell_12nE_6    = arrayfun(@(nE) zeros(nE * 12, 6)                    , n_elem   , 'Unif', false);   %[nElem * 12, 6]

%Some finite element matrices...
Dtot        = cell_6nN_6nE;
Btot        = cell_6nN_6nN;
Btot_mp     = cell_6nN_6nE;
Btot_mp_pnt = cell_6nN_6nE;
B_aero      = cell_6nN_6nP;
e1tot       = cell_6nE_1;
Gvtot       = cell_12nN_12nN;
Gftot       = cell_12nN_12nN;
G_bc        = cell_12nE_6;

R1_bc = arrayfun(@(nE) [zeros(nE * 6, 6), eye(nE * 6)], n_elem, 'Unif', false);
R2_bc = arrayfun(@(nE) [eye(nE * 6), zeros(nE * 6, 6)], n_elem, 'Unif', false);

vel_in = arrayfun(@(nN) zeros(nN * 6, 6), n_node, 'Unif', false);  

Mtot = arrayfun(@(i) zeros(n_node(i) * 6, n_elem(i) * 6), 1 : nPart, 'Unif', false);
Ftot = arrayfun(@(i) zeros(n_node(i) * 6, n_elem(i) * 6), 1 : nPart, 'Unif', false);

loadInt     = arrayfun(@(nN) zeros(6, nN * 6) , n_node  , 'Unif', false);
loadIntAero = arrayfun(@(nP) zeros(6, nP * 6) , n_points, 'Unif', false);
PntInt      = arrayfun(@(nE) zeros(6, nE * 12), n_elem  , 'Unif', false);

%% * * * Possibly remove...

%Ignore this - Different formulation
e1tilde = skew([1;0;0]);
e2tilde = skew([0;1;0]);
e3tilde = skew([0;0;1]);
E1tilde = [ zeros(3), -e1tilde,   zeros(3),-e2tilde,  zeros(3), -e3tilde,   zeros(3),  zeros(3), zeros(3),  zeros(3), zeros(3),  zeros(3);
    -e1tilde,   zeros(3), -e2tilde,  zeros(3),-e3tilde,   zeros(3),  zeros(3), -e1tilde,  zeros(3), -e2tilde,  zeros(3), -e3tilde];
E2tilde = [ zeros(3),  zeros(3), zeros(3),  zeros(3), zeros(3),  zeros(3), -e1tilde,   zeros(3), -e2tilde,  zeros(3),-e3tilde,   zeros(3);
    -e1tilde,   zeros(3), -e2tilde,  zeros(3),-e3tilde,   zeros(3),  zeros(3), -e1tilde,   zeros(3),-e2tilde,  zeros(3), -e3tilde];
E3tilde = [ zeros(3), -e1tilde,  zeros(3), -e2tilde,  zeros(3), -e3tilde,  -e1tilde,   zeros(3), -e2tilde,  zeros(3),-e3tilde,   zeros(3);
    zeros(3),  zeros(3), zeros(3),  zeros(3), zeros(3),  zeros(3),  zeros(3), -e1tilde,  zeros(3), -e2tilde,  zeros(3), -e3tilde];
P_36    = sparse(36,36);
for ii = 1:6
    ind6_6       = [1:6] + (ii-1)*6;
    for kk = 1:6
        ind6_2 = [1:6] + (kk-1)*6;
        P_36(ind6_6(kk),ind6_2(ii)) = 1;
    end
end
%^^^ Can get rid of this.

% if ~exist('p0','var') || isempty(p0{jj})
%     p0     = cell(1,nPart);
% end
if ~exist('eps0','var') %TODO - Try to avoid using 'exist' where possible.
    eps0   = cell_1_nP;
end
if ~exist('kappa0','var')
    kappa0 = cell_1_nP;
end
if ~exist('CBB','var')
    CBB = cell_1_nP;
end
if ~exist('Parent','var')
    Parent = cell_1_nP;
end
if ~exist('CFFbar','var')
    CFFbar = cell(1,nPart);
end
if ~exist('M_pt_store','var')
    M_pt_store = cell_1_nP;
end

Children = cell_1_nP;

%% Set up FE matrices for each part

for jj = 1 : nPart
    
    nElem    = n_elem(jj);    
    partName = lower(allPartIDs{jj});    
    fprintf(['Setting up finite element matrices for Member ', ...
        num2str(jj),' (''', partName, ''')...\n'])
    
    %Some finite element matrices...
    Dtot2{jj}   = zeros(nElem*12  ,nElem*12);
    Btot2{jj}   = zeros(nElem*12  ,nElem*12);
    e1tot2{jj}  = zeros(nElem*12  ,1);
    
    %Is this applying a boundary condition??
    G_bc{jj}(1 : 6,:) = eye(6);
    
    %Logical flag
    if strcmp(type{jj},'LS')
        isLS{jj} = 1;
    elseif strcmp(type{jj},'FUS')
        isLS{jj} = 0;
    else
        isLS{jj} = 0;
    end
        
    %Get rid of this...
    if ~exist('p0','var')     || isempty(p0{jj})
        p0{jj} = [0;0;0];
    end
    if ~exist('eps0','var')   || isempty(eps0{jj})
        eps0{jj}   = zeros(3,1,nElem);
    end
    if ~exist('kappa0','var') || isempty(kappa0{jj})
        kappa0{jj} = zeros(3,1,nElem);
    end
    if ~exist('CBB','var')    || isempty(CBB{jj})
        CBB{jj} = repmat(eye(3),1,1,nElem);
    end
    if ~exist('Parent','var') || isempty(Parent{jj})
        Parent{jj} = 0;
    end
    if ~exist('CFFbar','var') || isempty(CFFbar{jj})
        CFFbar{jj} = repmat(eye(6),1,1,nElem);
    end
    if isempty(M_pt_store{jj})
        M_pt_store{jj} = zeros(6,6,nElem+1);
    end
    
    %Is this ever actually used?
    if Parent{jj}(1)
        Children{Parent{jj}(1)} = [Children{Parent{jj}(1)};jj];
        parent_ind               = Parent{jj}(1);
        ind_tmp                  = [1:6] + (Parent{jj}(2)-2)*6;
        R3_bc{jj}                = sparse(zeros(n_elem(parent_ind)*6,(nElem+1)*6));
        R3_bc{jj}(ind_tmp,1:6)   = blkdiag(CaB0{Parent{jj}(1)}'*CaB0{jj},CaB0{Parent{jj}(1)}'*CaB0{jj});
        R4_bc{jj}                = sparse(zeros(nElem*6,(n_elem(parent_ind)+1)*6));
        R4_bc{jj}(1:6,end-5:end) = blkdiag(CaB0{Parent{jj}(1)}'*CaB0{jj},CaB0{Parent{jj}(1)}'*CaB0{jj})';
    end
        
    %??????
    vel_in{jj}(1:6,1:6) = [CaB0{jj}',zeros(3);zeros(3),CaB0{jj}']*[eye(3),-skew(p0{jj});zeros(3),eye(3)];
    
    %Mtot{jj}    = zeros((nElem+1)*6, nElem*6 );
    Mtot2{jj}   = zeros( nElem   *12,nElem*12);
    %Ftot{jj}    = zeros((nElem+1)*6, nElem*6 );
    Ftot2{jj}   = zeros( nElem   *12,nElem*12);
    
    %loadInt{jj}     = zeros(6,(nElem+1)*6);
    loadInt2{jj}    = zeros(6, nElem   *12);
    %loadIntAero{jj} = zeros(6,length(s_out{jj})*6);
    %PntInt{jj}      = zeros(6, nElem   *12);    
    
    %What is all this???
    B1 = zeros(nElem*3+6,nElem*3+3);
    B1(1:end-3,:) = eye(nElem*3+3);
    B1(4:end,:)   = B1(4:end,:) + eye(nElem*3+3);
    B1 = B1(1:end-3,:);
    
    B2 = zeros(nElem*3+3,nElem*3);
    B2(4:end,:)   = 2*eye(nElem*3);
    
    x_p_recov1{jj} = inv(B1)*B2;
    x_p_recov2{jj} = zeros(nElem*3+3,3);
    x_p_recov2{jj}(1:3,:) = eye(3);
    x_p_recov2{jj} = inv(B1)*x_p_recov2{jj};
    
    B1 = zeros(nElem*6+12,nElem*6+6);
    B1(1:end-6,:) = eye(nElem*6+6);
    B1(7:end,:)   = B1(7:end,:) + eye(nElem*6+6);
    B1 = B1(1:end-6,:);
    
    B2 = zeros(nElem*6+6,nElem*6);
    B2(7:end,:)   = 2*eye(nElem*6);
    
    x_v_recov1_tmp = inv(B1)*B2;
    x_v_recov2_tmp = zeros(nElem*6+6,6);
    x_v_recov2_tmp(1:6,:) = eye(6);
    x_v_recov2_tmp = inv(B1)*x_v_recov2_tmp;
    
    x_v_recov1{jj} = [];
    x_v_recov2{jj} = [];
    for ii = 1 : n_node(jj)
        ind = [1 : 6] + (ii - 1) * 6;
        orig_v         = zeros(6, nElem*6);
        orig_v(:,ind)  = eye(6);
        if ii == nElem+1
            x_v_recov1{jj} = [x_v_recov1{jj};x_v_recov1_tmp(ind,:)];
            x_v_recov2{jj} = [x_v_recov2{jj}();x_v_recov2_tmp(ind,:)];
        else
            x_v_recov1{jj} = [x_v_recov1{jj};x_v_recov1_tmp(ind,:);orig_v];
            x_v_recov2{jj} = [x_v_recov2{jj}();x_v_recov2_tmp(ind,:);zeros(6)];
        end
    end
    
    %Does this ever get used?
    if extra_flag
        Q1_quad{jj}  = sparse(zeros((nElem+1)*6,  nElem*6    ));
        Q2_quad{jj}  = sparse(zeros((nElem+1)*6,  nElem*6    ));
        Q3_quad{jj}  = sparse(zeros((nElem+1)*6,  nElem*6    ));
        Q4_quad{jj}  = sparse(zeros((nElem+1)*6,  nElem*6    ));
        Q5_quad{jj}  = sparse(zeros((nElem+1)*6,  nElem*6    ));
        Q6_quad{jj}  = sparse(zeros((nElem+1)*6,  nElem*6    ));
        Q7_quad{jj}  = sparse(zeros((nElem+1)*6,  nElem*6    ));
        Q8_quad{jj}  = sparse(zeros((nElem+1)*6,  nElem*6    ));
        Q9_quad{jj}  = sparse(zeros((nElem+1)*6,  nElem^2*36 ));
        Q10_quad{jj} = sparse(zeros((nElem+1)*6,  nElem^2*36 ));
        Q11_quad{jj} = sparse(zeros((nElem+1)*6,  nElem^2*36 ));
        Q12_quad{jj} = sparse(zeros((nElem+1)*6,  nElem^2*36 ));
        Q13_quad{jj} = sparse(zeros((nElem+1)*6,  nElem^2*36 ));
        Q14_quad{jj} = sparse(zeros((nElem+1)*6,  nElem^2*36 ));
        Q15_quad{jj} = sparse(zeros((nElem+1)*6,  nElem^2*36 ));
        Q16_quad{jj} = sparse(zeros((nElem+1)*6,  nElem^2*36 ));
        P_quad       = sparse(       nElem^2*144, nElem^2*144);
        %         P_quad_test  = sparse(       nElem^2*144, nElem^2*144);
    end
    
    %% Start populating FE matrices
    %   - For a detailed description of this process the user is directed
    %   to Section II.D.2 of the following paper: "Aeroelastic Modelling of
    %   Highly Flexible Wings", C.Howcroft et. al, AIAA SciTech Forum, 2016
    %   
    ind5_start  = 1;
    
    %Indexing
    ind6_6   = arrayfun(@(i) (1 : 6) + ((i - 1) * 6)   , 1 : nElem, 'Unif', false);
    ind6_6   = vertcat(ind6_6{:});
    ind12_12 = arrayfun(@(i) (1 : 12) + ((i - 1) * 12), 1 : nElem, 'Unif', false);
    ind12_12 = vertcat(ind12_12{:});
    ind12_6  = arrayfun(@(i) (1 : 12) + ((i - 1) * 6) , 1 : nElem, 'Unif', false);
    ind12_6  = vertcat(ind12_6{:});
    
    for ii = 1 : nElem
        
        %System matrices without shape function terms
        %   - mass matrix
        massMatrix = [ ...
            m{jj}(ii) * eye(3)              , -m{jj}(ii) * skew(cg{jj}(:,:,ii));
            m{jj}(ii) * skew(cg{jj}(:,:,ii)),  eye(3)    * J{jj}(:,:,ii)      ];        
        %   - stiffness matrix 
        complianceMatrix  = C{jj}(:,:,ii);
        
        %Equation 2 of Robbie's paper - Extension is along the local x!
        e1 = zeros(6,1); e1(1) = 1;
        %   - eps0 = pre strain, kappa0 = pre-curvature
        e1 = e1 + [eps0{jj}(:,:,ii) ; kappa0{jj}(:,:,ii)];
        e1 = K{jj}(:,:,ii) * e1;
        
        %Incremental beam length
        ds = s{jj}(ii + 1) - s{jj}(ii);
        
        %Some kind of rotation of the geometry...
        %   - CBB is the rotation matric from one element to the next.
        CBB_tot = blkdiag(CBB{jj}(:,:,ii), CBB{jj}(:,:,ii), eye(6));
        %   - Linear shape functions for the FE formulation
        %   - See eqn 18
        A1 = [-eye(6)*ds,  zeros(6)]/(-ds * CBB_tot);
        A2 = [ eye(6)   , -eye(6)  ]/(-ds * CBB_tot);
        
        Ab  = A1'*ds + A2'*ds^2/2; %Coefficient of element stiffness matrix
        Ab2 = eye(6)*ds;
        
        %% Shape functions for Galerkin FE element 
        Ab1var  = (A1'*A1*ds     + A1'*A2*ds^2/2 + ...
            A2'*A1*ds^2/2 + A2'*A2*ds^3/3);
        
        Ab2var  = (A1'*A1*ds^2/2 + A1'*A2*ds^3/3 + ...
            A2'*A1*ds^3/3 + A2'*A2*ds^4/4);
        
        Ab3var  = (A1*ds     + A2*ds^2/2);
        
        Ab4var  = (A1*ds^2/2 + A2*ds^3/3); %Should this be over 6?
        
        %Building the state matrices for nonlinear corriolis terms, etc.
        Aa1var = complianceMatrix * A1;  %Stiffness
        Aa2var = complianceMatrix * A2;
        Aa3var = massMatrix  * A1;  %Mass
        Aa4var = massMatrix  * A2;
        Aa5var = A1;                 %Some system matrix...
        Aa6var = A2;
        %% System matrices
        
        %
        Bdst = A1'*A1*ds     + A1'*A2*ds^2/2 + ...
            A2'*A1*ds^2/2 + A2'*A2*ds^3/3;
        
        %Minus A2 because of integration by parts
        D   = -A2'*eye(6)*ds;
        %         D2  = -(A2'*A1*ds + A2'*A2*ds^2/2);
        D2  = (A1'*A2*ds + A2'*A2*ds^2/2);
        
        elemMassMatrix  = Ab  * massMatrix;     %Elemental mass matrix        
        elemStiffMatrix = Ab  * complianceMatrix;    %Elemental stiffness matrix
        
        %Stuff below is for alternative formulation (I think...)
        Ms2 = ...
            A1' * massMatrix * A1 * ds     + A1' * massMatrix * A2 * ds^2/2 + ...
            A2' * massMatrix * A1 * ds^2/2 + A2' * massMatrix * A2 * ds^3/3;
        Fs2 = ...
            A1' * complianceMatrix * A1 * ds     + A1' * complianceMatrix * A2 * ds^2/2 + ...
            A2' * complianceMatrix * A1 * ds^2/2 + A2' * complianceMatrix * A2 * ds^3/3;
        
        C_dihedral_v = [cos(local_dihedral{jj}(ii)) 0 -sin(local_dihedral{jj}(ii)) ; 0 1 0 ; sin(local_dihedral{jj}(ii)) 0 cos(local_dihedral{jj}(ii))];
        C_sweep_v    = [cos(local_sweep{jj}(ii)) -sin(local_sweep{jj}(ii)) 0 ; sin(local_sweep{jj}(ii)) cos(local_sweep{jj}(ii)) 0 ; 0 0 1];
        C_twist_v    = [1 0 0 ; 0 cos(local_twist{jj}(ii)) sin(local_twist{jj}(ii)) ; 0 -sin(local_twist{jj}(ii)) cos(local_twist{jj}(ii))];
        
        C_rot_v      = C_dihedral_v*C_sweep_v*C_twist_v;
        
        C_dihedral_f = [cos(local_dihedral{jj}(ii+1)) 0 -sin(local_dihedral{jj}(ii+1)) ; 0 1 0 ; sin(local_dihedral{jj}(ii+1)) 0 cos(local_dihedral{jj}(ii+1))];
        C_sweep_f    = [cos(local_sweep{jj}(ii+1)) -sin(local_sweep{jj}(ii+1)) 0 ; sin(local_sweep{jj}(ii+1)) cos(local_sweep{jj}(ii+1)) 0 ; 0 0 1];
        C_twist_f    = [1 0 0 ; 0 cos(local_twist{jj}(ii+1)) sin(local_twist{jj}(ii+1)) ; 0 -sin(local_twist{jj}(ii+1)) cos(local_twist{jj}(ii+1))];
        
        C_rot_f      = C_dihedral_f*C_sweep_f*C_twist_f;
        
        C_rot{jj}(:,:,ii) = C_dihedral_v*C_sweep_v*C_twist_v';
        if ii == nElem
            C_rot{jj}(:,:,ii+1) = C_dihedral_f*C_sweep_f*C_twist_f';
        end
        
        %         G   = [zeros(6),zeros(6);eye(6),-eye(6)];
        Gv   = [zeros(6),zeros(6);eye(6),-blkdiag(C_rot_v,C_rot_v)];
        Gf   = [zeros(6),zeros(6);eye(6),-blkdiag(C_rot_f,C_rot_f)];

        %What is this all about?
        nodes_in_elem = find(s_out{jj}>=s{jj}(ii)&s_out{jj}<=s{jj}(ii+1));
        nNodesInElem  = numel(nodes_in_elem);
        B_aero_e      = zeros(12, nNodesInElem * 6);
        loadIntAero_e = zeros( 6, nNodesInElem * 6);
        for kk = 1 : nNodesInElem - 1
            ind4 = [1:12] + (kk-1)*6;
            ds_elem    = s_out{jj}(nodes_in_elem(kk+1))-s_out{jj}(nodes_in_elem(kk));
            B_aero_e_n = A1'*A1*(ds_elem)     + A1'*A2*(ds_elem)^2/2 + ...
                A2'*A1*(ds_elem)^2/2 + A2'*A2*(ds_elem)^3/3;
            
            B_aero_e(:,ind4)      = B_aero_e(:,ind4) + B_aero_e_n;
            loadIntAero_e(:,ind4) = loadIntAero_e(:,ind4) + Ab';
        end
        ind5       = ind5_start:ind5_start+nNodesInElem*6-1;
        ind5_start = ind5_start+nNodesInElem*6-6;
        
        %% Alt cross method (using kron)
        MkronI                     = kron(massMatrix,                        eye(6));
        CkronI                     = kron(C{jj}(:,:,ii),            eye(6));
        Ikronef                    = kron(eye(6),                   e1    );
        MptkronI                   = kron(M_pt_store{jj}(:,:,ii+1), eye(6));
        
        e1ds                       = [1;0;0]*ds/2;
        T_pt{jj}(:,:,ii)           = [eye(3) -skew(e1ds);zeros(3) eye(3)];
        T_ptkronT_pt               = kron(T_pt{jj}(:,:,ii),T_pt{jj}(:,:,ii));
        
        C_pre1{jj}(:,:,ii)         = Ab*E1tilde*MkronI;
        C_pre2{jj}(:,:,ii)         = C_pre1{jj}(:,:,ii)*P_36;
        A_pre1{jj}(:,:,ii)         = Ab*E2tilde*CkronI;
        A_pre2{jj}(:,:,ii)         = A_pre1{jj}(:,:,ii)*P_36;
        A_pre3{jj}(:,:,ii)         = A_pre2{jj}(:,:,ii)*Ikronef;
        E_pre1{jj}(:,:,ii)         = Ab*E3tilde*CkronI;
        E_pre2{jj}(:,:,ii)         = E_pre1{jj}(:,:,ii)*P_36;
        E_pre3{jj}(:,:,ii)         = E_pre2{jj}(:,:,ii)*Ikronef;
        E2tildeIkronM{jj}(:,:,ii)  = E2tilde*kron(eye(6),massMatrix);
        
        C_pt_pre1{jj}(:,:,ii)      = E1tilde*MptkronI*T_ptkronT_pt;
        C_pt_pre2{jj}(:,:,ii)      = C_pt_pre1{jj}(:,:,ii)*P_36;
        
        %% Combine Everything
        Ab_sep{jj}(:, :, ii)  = Ab;
        Ab2_sep{jj}(:, :, ii) = Ab2;
        Aa_sep{jj}(:, :, ii)  = complianceMatrix;
        Aa2_sep{jj}(:, :, ii) = massMatrix;
        
        %Shape functions for nonlinear terms (state * state)
        Ab1var_sep{jj}(:,:,ii)  = Ab1var;
        Ab2var_sep{jj}(:,:,ii)  = Ab2var;
        Ab3var_sep{jj}(:,:,ii)  = Ab3var;
        Ab4var_sep{jj}(:,:,ii)  = Ab4var;
        Aa1var_sep{jj}(:,:,ii)  = Aa1var;
        Aa2var_sep{jj}(:,:,ii)  = Aa2var;
        Aa3var_sep{jj}(:,:,ii)  = Aa3var;
        Aa4var_sep{jj}(:,:,ii)  = Aa4var;
        Aa5var_sep{jj}(:,:,ii)  = Aa5var;
        Aa6var_sep{jj}(:,:,ii)  = Aa6var;
        
        %Replace with lb/ub method (cumsum)
        ind1 = ind12_12(ii, :);
        ind2 = ind12_6(ii, :);
        ind3 = ind6_6(ii, :);

        Dtot{jj}(ind2,ind3)        = Dtot{jj}(ind2,ind3)     + D; %Differentiation matrix for differentiating the states
        Dtot2{jj}(ind1,ind1)       = Dtot2{jj}(ind1,ind1)    + D2;
        Btot{jj}(ind2,ind2)        = Btot{jj}(ind2,ind2)     + Bdst;%Applies a force  
        Btot2{jj}(ind1,ind1)       = Btot2{jj}(ind1,ind1)    + Bdst;%
        Btot_mp{jj}(ind2,ind3)     = Btot_mp{jj}(ind2,ind3)  + Ab;% Applies a distributed force to the midpoints
        Btot_mp_pnt{jj}(ind2,ind3) = Btot_mp_pnt{jj}(ind2,ind3)  + [eye(6);eye(6)]/2;%Point force applied to the midpoints
        B_aero{jj}(ind2,ind5)      = B_aero{jj}(ind2,ind5)   + B_aero_e;%Possibly not used anymore
        e1tot{jj}(ind3,1)          = e1; %Vector of precurvature/strains
        e1tot2{jj}(ind1,1)         = [e1;e1];
        Mtot{jj}(ind2,ind3)        = Mtot{jj}(ind2,ind3)      + elemMassMatrix; %Mass Matrix
        Mtot2{jj}(ind1,ind1)       = Mtot2{jj}(ind1,ind1)     + Ms2;% + M_pnt{jj}(ind1,ind1);
        Ftot{jj}(ind2,ind3)        = Ftot{jj}(ind2,ind3)      + elemStiffMatrix; %Compliance matrix???? (Including shape functions)
        Ftot2{jj}(ind1,ind1)       = Ftot2{jj}(ind1,ind1)     + Fs2;
        Gvtot{jj}(ind1+6,ind1+6)   = Gvtot{jj}(ind1+6,ind1+6) + Gv; %For the other formulation
        Gftot{jj}(ind1+6,ind1+6)   = Gftot{jj}(ind1+6,ind1+6) + Gf;
        
        %If point masses are included adjust the mass matrix accordingly
        Mtot{jj}(ind3+6,ind3)      = Mtot{jj}(ind3+6,ind3)     + M_pt_store{jj}(:,:,ii+1)*T_pt{jj}(:,:,ii);
        
        loadInt{jj}(1:6,ind2)      = loadInt{jj}(1:6,ind2)     + Ab'; %Loads integration?? During post-processing
        loadInt2{jj}(1:6,ind1)     = loadInt2{jj}(1:6,ind1)    + Ab';
        loadIntAero{jj}(1:6,ind5)  = loadIntAero{jj}(1:6,ind5) + loadIntAero_e;
        PntInt{jj}(1:6,ind1)       = [eye(6),eye(6)];
        
        if extra_flag
            ind36 = [1:36] + (ii-1)*36;
            Q1_quad{jj}( ind2,ind3 ) = Q1_quad{jj}( ind2,ind3 ) + elemMassMatrix;
            Q2_quad{jj}( ind2,ind3 ) = Q2_quad{jj}( ind2,ind3 ) + D*1e-6;
            Q4_quad{jj}( ind2,ind3 ) = Q4_quad{jj}( ind2,ind3 ) + elemStiffMatrix;
            Q6_quad{jj}( ind2,ind3 ) = Q6_quad{jj}( ind2,ind3 ) - D;
            Q7_quad{jj}( ind2,ind3 ) = Q7_quad{jj}( ind2,ind3 ) - D;
            Q9_quad{jj}( ind2,ind36) = Q9_quad{jj}( ind2,ind36) + Ab*E1tilde*kron(massMatrix,            eye(6));
            Q12_quad{jj}(ind2,ind36) = Q12_quad{jj}(ind2,ind36) + Ab*E2tilde*kron(C{jj}(:,:,ii),eye(6));
            Q14_quad{jj}(ind2,ind36) = Q14_quad{jj}(ind2,ind36) + Ab*E3tilde*kron(C{jj}(:,:,ii),eye(6));
        end
    end
    if extra_flag
        for ii = 1:(12*nElem)
            ind192       = [1:12*nElem] + (ii-1)*12*nElem;
            for kk = 1:(12*nElem)
                ind192_2 = [1:12*nElem] + (kk-1)*12*nElem;
                P_quad(ind192(kk),ind192_2(ii)) = 1;
            end
        end
    end
    
    %Block diagonal of compliance matrices (no shape functions)
    %   - Multiply the force & moment vector by 'Cfull' to obtain the
    %     curvatures & strains
    Cdiag     = mat2cell(C{jj}, 6, 6, ones(1, nElem));
    Cfull{jj} = blkdiag(Cdiag{:});
    
    %Block diagonal of shear stiffness matrices (no shape functions)
    CFFdiag         = mat2cell(CFFbar{jj}, 6, 6, ones(1, nElem));
    CFFbar_full{jj} = blkdiag(CFFdiag{:});
    
end

if ~exist('CBA0','var')
    CBA0 = [];
end
if ~exist('le','var')
    le = [];
end
if ~exist('te','var')
    te = [];
end
if ~exist('cp','var')
    cp = [];
end
if ~exist('RenderPoints','var')
    RenderPoints = [];
end
if ~exist('RenderPoints_RB','var')
    RenderPoints_RB = [];
end
if ~exist('dCLdd','var')
    dCLdd = [];
end
if ~exist('dCMdd','var')
    dCMdd = [];
end
if ~exist('delta_flap','var')
    delta_flap = [];
end

%Store in 'Matrices' structure
%   - TODO: Define as sparse...

Matrices.grav_switch     = 1;   %yes
    
Matrices.f        =  f;         %yes
Matrices.f_a      =  f_a;       %yes
Matrices.f_mp_pnt =  f_mp_pnt;  %yes

Matrices.Ab  = Ab_sep;  %yes
Matrices.Ab2 = Ab2_sep; 
Matrices.Aa  = Aa_sep;  %yes
Matrices.Aa2 = Aa2_sep; %yes

Matrices.Ab1var = Ab1var_sep;
Matrices.Ab2var = Ab2var_sep;
Matrices.Ab3var = Ab3var_sep;
Matrices.Ab4var = Ab4var_sep;
Matrices.Aa1var = Aa1var_sep;
Matrices.Aa2var = Aa2var_sep;
Matrices.Aa3var = Aa3var_sep;
Matrices.Aa4var = Aa4var_sep;
Matrices.Aa5var = Aa5var_sep;
Matrices.Aa6var = Aa6var_sep;

Matrices.C_pre1    = C_pre1;
Matrices.C_pre2    = C_pre2;
Matrices.A_pre1    = A_pre1;
Matrices.A_pre2    = A_pre2;
Matrices.A_pre3    = A_pre3;
Matrices.E_pre1    = E_pre1;
Matrices.E_pre2    = E_pre2;
Matrices.E_pre3    = E_pre3;
Matrices.C_pt_pre1 = C_pt_pre1;
Matrices.C_pt_pre2 = C_pt_pre2;
Matrices.E1tilde = E1tilde;
% Matrices.E2tilde = E2tilde;
Matrices.E2tildeIkronM  = E2tildeIkronM;

Matrices.Dtot    = Dtot;    %yes
Matrices.Dtot2   = Dtot2;
Matrices.Btot    = Btot;
Matrices.Btot2   = Btot2;
Matrices.Btot_mp = Btot_mp; %yes
Matrices.Btot_mp_pnt = Btot_mp_pnt;
Matrices.B_aero  = B_aero;

Matrices.Gvtot    = Gvtot;
Matrices.Gftot    = Gftot;
Matrices.G_bc     = G_bc;

Matrices.R1_bc    = R1_bc;  %yes
Matrices.R2_bc    = R2_bc;
Matrices.R3_bc    = R3_bc;  %yes
Matrices.R4_bc    = R4_bc;

Matrices.vel_in   = vel_in;

Matrices.CFFbar   = CFFbar_full;

Matrices.M   = Mtot;
Matrices.M2  = Mtot2;%M_pnt;%
Matrices.T1  = Ftot;
Matrices.T12 = Ftot2;

Matrices.e1tot  = e1tot;    %yes
Matrices.e1tot2 = e1tot2;

Matrices.eps0   = eps0;     %yes
Matrices.kappa0 = kappa0;   %yes    

Matrices.M_pt = M_pt_store; %yes
Matrices.T_pt = T_pt;

Matrices.m  = m;
Matrices.J  = J;
Matrices.cg = cg;

Matrices.C = C;     %yes

Matrices.C_rot = C_rot;

Matrices.loadInt     = loadInt;
Matrices.loadInt2    = loadInt2;
Matrices.loadIntAero = loadIntAero;
Matrices.pntInt      = PntInt;

Matrices.Cfull = Cfull;     %yes

Matrices.n_elem = n_elem;   %yes
Matrices.n_node = n_node;   %yes
Matrices.l      = l;        %yes

Matrices.type   = type;

Matrices.s             = s;             %yes
Matrices.s2            = s2;
Matrices.s_node_ind    = s_node_ind;    %yes
Matrices.s_node_ind2   = s_node_ind2;
Matrices.s_aero        = s_aero;        %yes
Matrices.s_aero_mp     = s_aero_mp;     %yes
Matrices.s_out         = s_out;         %yes
Matrices.s_out2        = s_out2;
Matrices.s_mp_ind      = s_mp_ind;      %yes
Matrices.s_aero_ind    = s_aero_ind;    %yes    
Matrices.s_aero_mp_ind = s_aero_mp_ind; %yes

Matrices.p0            = p0;    %yes
Matrices.CaB0          = CaB0;  %yes
Matrices.CBA0          = CBA0;  %yes
Matrices.CBB           = CBB;   %yes

Matrices.chord         = chord;
Matrices.le            = le;
Matrices.te            = te;
Matrices.cp            = cp;
% Matrices.xs            = xs;

Matrices.RenderPoints     = RenderPoints;
Matrices.RenderPoints_RB  = RenderPoints_RB;

Matrices.x_p_recov1    = x_p_recov1;
Matrices.x_p_recov2    = x_p_recov2;
Matrices.x_v_recov1    = x_v_recov1;
Matrices.x_v_recov2    = x_v_recov2;

Matrices.Parent        = Parent;    %yes
Matrices.Children      = Children;  %yes

if extra_flag
    reorder_mat1            = sparse(36*n_elem^2,144*n_elem^2);
    reorder_mat2            = sparse(36*n_elem^2,144*n_elem^2);
    reorder_mat3            = sparse(36*n_elem^2,144*n_elem^2);
    reorder_mat4            = sparse(36*n_elem^2,144*n_elem^2);
    
    counter = 0;
    for ii = 1:n_elem
        ind6_6 = [1:6] + (ii-1)*6;
        for jj = 1:6
            VO_tmp1           = zeros(6*n_elem,1);
            FO_tmp1           = zeros(6*n_elem,1);
            VO_tmp1(ind6_6(jj)) = 1;
            for kk = 1:6
                counter = counter + 1;
                VO_tmp2                 = zeros(6*n_elem,1);
                FO_tmp2                 = zeros(6*n_elem,1);
                VO_tmp2(ind6_6(kk))       = 1;
                reorder_mat1(counter,:) = kron([VO_tmp1;FO_tmp1],[VO_tmp2;FO_tmp2]);
            end
        end
    end
    counter = 0;
    for ii = 1:n_elem
        ind6_6 = [1:6] + (ii-1)*6;
        for jj = 1:6
            VO_tmp1           = zeros(6*n_elem,1);
            FO_tmp1           = zeros(6*n_elem,1);
            FO_tmp1(ind6_6(jj)) = 1;
            for kk = 1:6
                counter = counter + 1;
                VO_tmp2                 = zeros(6*n_elem,1);
                FO_tmp2                 = zeros(6*n_elem,1);
                VO_tmp2(ind6_6(kk))       = 1;
                reorder_mat2(counter,:) = kron([VO_tmp1;FO_tmp1],[VO_tmp2;FO_tmp2]);
            end
        end
    end
    counter = 0;
    for ii = 1:n_elem
        ind6_6 = [1:6] + (ii-1)*6;
        for jj = 1:6
            VO_tmp1           = zeros(6*n_elem,1);
            FO_tmp1           = zeros(6*n_elem,1);
            FO_tmp1(ind6_6(jj)) = 1;
            for kk = 1:6
                counter = counter + 1;
                VO_tmp2                 = zeros(6*n_elem,1);
                FO_tmp2                 = zeros(6*n_elem,1);
                FO_tmp2(ind6_6(kk))       = 1;
                reorder_mat4(counter,:) = kron([VO_tmp1;FO_tmp1],[VO_tmp2;FO_tmp2]);
            end
        end
    end
    reorder_mat             = [reorder_mat1;reorder_mat2;reorder_mat3;reorder_mat4];
    reorder_mat_red         = reorder_mat;
    reorder_mat_red(sum(reorder_mat,2)==0,:) = [];
    
    Q1_tot                  = [R1_bc{1}*Q1_quad{1}, R1_bc{1}*Q2_quad{1} ;...
        R2_bc{1}*Q3_quad{1}, R2_bc{1}*Q4_quad{1} ];
    Q2_tot                  = [R1_bc{1}*Q5_quad{1}, R1_bc{1}*Q6_quad{1} ;...
        R2_bc{1}*Q7_quad{1}, R2_bc{1}*Q8_quad{1} ];
    Q3_tot                  = [R1_bc{1}*Q9_quad{1}, R1_bc{1}*Q10_quad{1},R1_bc{1}*Q11_quad{1},R1_bc{1}*Q12_quad{1};...
        R2_bc{1}*Q13_quad{1},R2_bc{1}*Q14_quad{1},R2_bc{1}*Q15_quad{1},R2_bc{1}*Q16_quad{1}];
    Matrices.A_quad         = Q2_tot;%Q1_tot\Q2_tot;
    Matrices.B_quad         = Q3_tot*reorder_mat;%Q1_tot\Q3_tot*reorder_mat;
    Matrices.C_quad         = Q1_tot;
    Matrices.I_quad         = speye(size(Q1_tot));
    Matrices.P_quad         = P_quad;
    Matrices.B_x_P_quad_p_I = Matrices.B_quad*(P_quad + speye(size(P_quad)));
    Matrices.e1_quad        = sparse([zeros(size(e1tot{1}));e1tot{1}]);
    Matrices.A_quad         = Matrices.A_quad + Matrices.B_quad*kron(Matrices.e1_quad,Matrices.I_quad);
    Matrices.D_quad         = sparse(n_elem*12,(n_elem+1)*6);
    Matrices.D_quad(1:(n_elem+1)*6,1:(n_elem+1)*6) = speye((n_elem+1)*6);
end

%%
Aero.chord = chord;
Aero.xs    = xs;
Aero.tau   = tau;

Aero.lift_dist = lift_dist;

Aero.dCLdd = dCLdd;
Aero.dCMdd = dCMdd;
Aero.delta_flap = delta_flap;

Matrices.isLS  = isLS;

%% Initialise labels
loadLabels  = {'\epsilon_x','\epsilon_y','\epsilon_z','\kappa_x (m^{-1})','\kappa_y (m^{-1})','\kappa_z (m^{-1})'};
forceLabels = {'F_x (N)','F_y (N)','F_z (N)','M_x (Nm)','M_y (Nm)','M_z (Nm)'};

%% Calculate Mass Properties
MassProps = calculateCoM(x,Matrices,1);

%% Change Mass
% deltaMass              = 61441.5892 - MassProps.m;
% Matrices.M_rb(1:3,1:3) = Matrices.M_rb(1:3,1:3) + eye(3)*deltaMass;
% MassProps = calculateCoM(x,Matrices,1);

%% Clear all variables from the Workspace, but for x, Matrices, Aero and MassProps
% clearvars -except x Matrices MassProps Aero ff flex_vec

end
