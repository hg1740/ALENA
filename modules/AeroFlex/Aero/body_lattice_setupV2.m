%%BODY_LATTICE_SETUPV2
% 
% 
%   Author: Dario Calderon 
function BAERO = body_lattice_setupV2(BAERO, COORD)

nbaero = length(BAERO.ID);

% Define the number of sections along the circumference
NP = 30;
THETA = [0:2*pi/NP:2*pi];

nt = length(THETA);

for i=1:nbaero
    
    Nsec = length(BAERO.geo.fs{i})-1;
    
    Ref = BAERO.geo.ref_point(i,:);
    
    % input geometry data
    Xgeo = BAERO.geo.fs{i}* BAERO.geo.L(i);
    R    = BAERO.geo.Rs{i};
    Rmat = eye(3);
    
    Corners = zeros(nt*(Nsec+1),3);
    for k = 1:nt
        Corners((k-1)*(Nsec+1) + 1:k*(Nsec+1),1) = Ref(1) + Xgeo;
        Corners((k-1)*(Nsec+1) + 1:k*(Nsec+1),2) = Ref(2) + R * sin(THETA(k));
        Corners((k-1)*(Nsec+1) + 1:k*(Nsec+1),3) = Ref(3) + R * cos(THETA(k));
    end
    
    count = 0;
    CornerIdx = zeros(Nsec*(nt-1),4);
    for idx_l = 1:Nsec
        for idx_r = 1:nt-1
            count = count+1;
            CornerIdx(count,:) = [idx_r + (idx_r-1)*(Nsec), idx_r + idx_r*(Nsec)+1, ...
                idx_r + idx_r*(Nsec) + 2 ,idx_r + (idx_r-1)*(Nsec) + 1] + (idx_l-1);
        end
    end
    
    NODE_GLOBAL = zeros(nt*(Nsec+1),3);
    for k = 1:size(Corners,1)
        NODE_GLOBAL(k,:) = Ref + Corners(k,:);
    end
    
    npanels = size(CornerIdx,1);
    PanelCentroid = zeros(npanels, 3);
    for k=1:npanels
        PanelCentroid(k,:) = sum(NODE_GLOBAL(CornerIdx(k,1:4),:),1)./4;
    end
    
    % Store the vertex coordinates and connecting panel nodes for patch
    BAERO.lattice.Elem.Node{i} = Corners;
    BAERO.lattice.Elem.Conn{i} = CornerIdx;
    BAERO.Rmat(:,:,i) = Rmat;
    BAERO.lattice.Elem.Midpoint_loc{i} = PanelCentroid;
end


