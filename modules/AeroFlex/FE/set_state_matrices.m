% The following script has been written in order to reduce the
% computational cost of the dynamic simulations by taking a vectorised
% approach to the solver


function [beam_model] = set_state_matrices(beam_model)

% Sectional Stiffness matrix
D = [];
% Interpolation matrix
NShape = [];
% Interpolation derivative matrix
NShapeD = [];
% Node Rotation Matrix (Collocation)
NR0 = [];
% Bar Rotation Matrix (Evaluation)
BARR0 = [];
% Bar Offsets Vecotr
Offsets = [];
% Strain Matrix
Strain = [];
% Node Coordinate points
NodeCoord = [];
% State to angle of attack conversion matrix
Cangle = [];
% Determine initial tangent vector in bar material reference frame
PO = [];

% Number of bar elements
nbar = beam_model.Info.nbeam;

Bar = beam_model.Beam;
Node = beam_model.Node;

rangle = [];rdof = [];rdofs1 = [];rdofs2 = [];rdofrot = [];rdoftr = [];
cangle = [];cdof = [];cdofs1 = [];cdofs2 = [];cdofrot = [];cdoftr = [];
vangle = [];vdof = [];vdofs1 = [];vdofs2 = [];vdofrot = [];vdoftr = [];

ndof  = double(beam_model.Info.ndof2);
ngrid = beam_model.Info.ngrid;
if beam_model.SolParam.Aeroelastic == 1
    npanels = size(beam_model.Aero.lattice_vlm.COLLOC,1);
else
    npanels = 0;
end

% Beam Interpolation matrices:

% COLLOC 1
eta = -1/sqrt(3);
NIp = [(0.5 * (2 * eta-1)) (-2 * eta) (0.5 * (2 * eta+1))]; % shape functions derivative evaluated in COLLOC1
NI = [(0.5 * eta * (eta-1)) (1-eta^2) (0.5 * eta * (eta+1))];
NI1 = NI(1) .* eye(3); NI2 = NI(2) .* eye(3); NI3 = NI(3) .* eye(3);

% COLLOC2
eta = +1/sqrt(3);
NIIp = [(0.5 * (2 * eta-1)),   (-2 * eta), (0.5 * (2 * eta+1))]; % shape functions derivative evaluated in COLLOC2
NII  = [(0.5 * eta * (eta-1)), (1-eta^2),  (0.5 * eta * (eta+1))];

NII1 = NII(1) .* eye(3); NII2 = NII(2) .* eye(3); NII3 = NII(3) .* eye(3);
NZero = zeros(3,3);

NShapeElem = [NI1,  NZero,NI2,  NZero,NI3,  NZero;...
    NZero,NI1,  NZero,NI2,  NZero,NI3;...
    NII1, NZero,NII2, NZero,NII3, NZero;...
    NZero,NII1, NZero,NII2, NZero,NII3];

for n = 1:nbar
    
    n1 = Bar.Conn(n, 1);
    n2 = Bar.Conn(n, 2);
    n3 = Bar.Conn(n, 3);
    
    NodeCoord = [NodeCoord;beam_model.Node.Coord(n1,:)';zeros(3,1);beam_model.Node.Coord(n2,:)';zeros(3,1);...
        beam_model.Node.Coord(n3,:)';zeros(3,1)];
    
    for DOF = 1:3
        if beam_model.Node.DOF2(n1,DOF)
            rdofs1 = [rdofs1,DOF + 18*(n-1)];
            cdofs1 = [cdofs1,double(beam_model.Node.DOF2(n1,DOF))];
            vdofs1 = [vdofs1,1];
        end
    end
    
    for DOF = 4:6
        if beam_model.Node.DOF2(n1,DOF)
            rdofs2 = [rdofs2,DOF + 18*(n-1)];
            cdofs2 = [cdofs2,double(beam_model.Node.DOF2(n1,DOF))];
            vdofs2 = [vdofs2,1];
        end
    end
    
    for DOF = 1:3
        if beam_model.Node.DOF2(n2,DOF)
            rdofs1 = [rdofs1,DOF + 18*(n-1) + 6];
            cdofs1 = [cdofs1,double(beam_model.Node.DOF2(n2,DOF))];
            vdofs1 = [vdofs1,1];
        end
    end
    
    for DOF = 4:6
        if beam_model.Node.DOF2(n2,DOF)
            rdofs2 = [rdofs2,DOF + 18*(n-1) + 6];
            cdofs2 = [cdofs2,double(beam_model.Node.DOF2(n2,DOF))];
            vdofs2 = [vdofs2,1];
        end
    end
    
    for DOF = 1:3
        if beam_model.Node.DOF2(n3,DOF)
            rdofs1 = [rdofs1,DOF + 18*(n-1) + 12];
            cdofs1 = [cdofs1,double(beam_model.Node.DOF2(n3,DOF))];
            vdofs1 = [vdofs1,1];
        end
    end
    
    for DOF = 4:6
        if beam_model.Node.DOF2(n3,DOF)
            rdofs2 = [rdofs2,DOF + 18*(n-1) + 12];
            cdofs2 = [cdofs2,double(beam_model.Node.DOF2(n3,DOF))];
            vdofs2 = [vdofs2,1];
        end
    end
    
    % offset global coords
    f1 = Node.R(:,:, n1) * Bar.Offset(n, 1:3)';
    f2 = Node.R(:,:, n2) * Bar.Offset(n, 4:6)';
    f3 = Node.R(:,:, n3) * Bar.Offset(n, 7:9)';
    
    x = [Node.Coord(n1,1)+f1(1) Node.Coord(n2,1)+f2(1) Node.Coord(n3,1)+f3(1)];
    y = [Node.Coord(n1,2)+f1(2) Node.Coord(n2,2)+f2(2) Node.Coord(n3,2)+f3(2)];
    z = [Node.Coord(n1,3)+f1(3) Node.Coord(n2,3)+f2(3) Node.Coord(n3,3)+f3(3)];
    
    JI(1) = (NIp* x');
    JI(2) = (NIp* y');
    JI(3) = (NIp* z');
    NIp = NIp ./ norm(JI);
    
    JII(1) = (NIIp* x');
    JII(2) = (NIIp* y');
    JII(3) = (NIIp* z');
    NIIp = NIIp ./ norm(JII);
    
    POv1 = [NIp* x',NIp* y',NIp* z'];
    POv2 = [NIIp* x',NIIp* y',NIIp* z'];
    
    PO = [PO;Bar.R(:,:,4,n)'*POv1';zeros(3,1);Bar.R(:,:,5,n)'*POv2';zeros(3,1)];
    
    % determine local initial PO
    
    NID1  = NIp(1) .* eye(3); NID2  = NIp(2) .* eye(3); NID3  = NIp(3) .* eye(3);
    NIID1 = NIIp(1).* eye(3); NIID2 = NIIp(2).* eye(3); NIID3 = NIIp(3).* eye(3);
    
    NShapeElemD = [NID1, NZero,NID2, NZero,NID3, NZero;...
        NZero,NID1, NZero,NID2, NZero,NID3;...
        NIID1,NZero,NIID2,NZero,NIID3,NZero;...
        NZero,NIID1,NZero,NIID2,NZero,NIID3];
    
    D1 = beam_model.Beam.D(:,:,1,n);
    D2 = beam_model.Beam.D(:,:,2,n);
    
    BR1 = beam_model.Beam.R(:,:,4,n);
    BR2 = beam_model.Beam.R(:,:,5,n);
    
    % Node rotations and offsets are repeated for the end nodes here!!
    NR0      = blkdiag(NR0,beam_model.Node.R(:,:,n1),beam_model.Node.R(:,:,n1),...
        beam_model.Node.R(:,:,n2),beam_model.Node.R(:,:,n2),...
        beam_model.Node.R(:,:,n3),beam_model.Node.R(:,:,n3));
    D       = blkdiag(D,D1,D2);
    BARR0    = blkdiag(BARR0,BR1,BR1,BR2,BR2);
    NShape  = blkdiag(NShape,NShapeElem);
    NShapeD = blkdiag(NShapeD,NShapeElemD);
    Offsets = [Offsets;Bar.Offset(n, 1:3)';zeros(3,1);Bar.Offset(n, 4:6)';zeros(3,1);Bar.Offset(n, 7:9)';zeros(3,1)];
end


% Associate a DOF to the RBE0 nodes in order to carry the state space over
% to the aerodynamic panels correctly
if beam_model.SolParam.Aeroelastic == 1
    beam_model.Node.DOFAero = beam_model.Node.DOF;
    beam_model.Info.ndofAero = beam_model.Info.ndof;
    beam_model.Aero.Interp.Ic;
    for gridn = 1:ngrid
        if ~isempty(beam_model.Node.Aero.Index(gridn).data)
            for i = 1:length(beam_model.Node.Aero.Index(gridn).data)
                beam_model.Node.DOFAero(beam_model.Node.Aero.Index(gridn).data(i),:) = [beam_model.Info.ndofAero+1 : beam_model.Info.ndofAero+6];
                beam_model.Info.ndofAero = beam_model.Info.ndofAero + 6;
            end
        end
    end
end


if beam_model.SolParam.Aeroelastic == 1
    
    if (beam_model.Info.spline_type > 1)
        IntCol = kron(beam_model.Aero.Interp.Ic,eye(6));
        IntFor = kron(beam_model.Aero.Interp.Imv,[eye(3),zeros(3)]);
    end
end

% Generate a matrix that carries over the displacement states over to the
% aerodynamic splining states:
% THIS COPIES THE DISPLACEMENTS ACROSS, FURTHER WORK NEEDS TO BE DONE TO
% THE DISPLACEMENT DUE TO ROTATION
DOF2STR = zeros(6*ngrid,beam_model.Info.ndof);
for n = 1:ngrid
    dof = beam_model.Node.DOF(n, 1:6);
    index = find(dof);
    for i = index
        DOF2STR(6*(n-1)+i,dof(i)) = 1;
        if beam_model.SolParam.Aeroelastic ==1
            if ~isempty(beam_model.Node.Aero.Index(n).data)
                for isp = 1:length(beam_model.Node.Aero.Index(n).data)
                    AeroIndex = beam_model.Node.Aero.Index(n).data(isp);
                    DOF2STR(6*(AeroIndex-1)+i,dof(i)) = 1;
                end
            end
        end
    end
end


ROT2STR = zeros(6*ngrid,beam_model.Info.ndof);
ROT2TRA = zeros(6*ngrid,beam_model.Info.ndof);
for n = 1:ngrid
    dof   = beam_model.Node.DOF(n, 1:6);
    index = find(dof(4:6));
    for i = index + 3
        ROT2STR(6*(n-1)+i,dof(i)) = 1;
        ROT2TRA(6*(n-1)+i-3,dof(i)) = 1;
        if beam_model.SolParam.Aeroelastic == 1
            if ~isempty(beam_model.Node.Aero.Index(n).data)
                for isp = 1:length(beam_model.Node.Aero.Index(n).data)
                    AeroIndex = beam_model.Node.Aero.Index(n).data(isp);
                    ROT2STR(6*(AeroIndex-1)+i,dof(i)) = 1;
                    ROT2TRA(6*(AeroIndex-1)+i-3,dof(i)) = 1;
                end
            end
        end
    end
end

if beam_model.SolParam.Aeroelastic == 1
    
    IntSpline = zeros(beam_model.Info.ndof,6*ngrid);
    SkAeroArm = zeros(6*ngrid,6*ngrid);
    r = [];c = [];v = [];
    for n = 1:ngrid
        
        dof = beam_model.Node.DOF(n,1:6);
        index = find(dof);
        for i = index
            IntSpline(dof(i),6*(n-1)+i) = 1;
            if ~isempty(beam_model.Node.Aero.Index(n).data)
                for isp = 1:length(beam_model.Node.Aero.Index(n).data)
                    AeroIndex = beam_model.Node.Aero.Index(n).data(isp);
                    IntSpline(dof(i),6*(AeroIndex-1)+i) = 1;
                end
            end
        end
        
        if ~isempty(beam_model.Node.Aero.Index(n).data)
            for isp = 1:length(beam_model.Node.Aero.Index(n).data)
                AeroIndex = beam_model.Node.Aero.Index(n).data(isp);
                SkAeroArm(6*(AeroIndex-1)+1:6*(AeroIndex-1)+3,6*(AeroIndex-1)+1:6*(AeroIndex-1)+3) = crossm(beam_model.Node.Aero.Coord(n).data(:,isp));
            end
        end
        
        r = [r,6*(n-1)+4:6*(n-1)+6];
        c = [c,6*(n-1)+1:6*(n-1)+3];
        v = [v,ones(1,3)];
    end
    
    Extpitrat = zeros(npanels,6*npanels);
    L2Vector = zeros(3*npanels,npanels);
    for n = 1:npanels
        L2Vector(3*(n-1)+1:3*(n-1)+3,n) = 1;
        Extpitrat(n,6*(n-1)+5) = 1;
    end
        
    SumAeroArm = sparse(r,c,v,6*ngrid,6*ngrid);
    
    % Generate a moment arm matrix from the aerodynamic moment about the center
    % of gravity
    if size(beam_model.Aero.lattice_vlm.VORTEX,2) == 8
        b1 = 4;
    elseif size(beam_model.Aero.lattice_vlm.VORTEX,2) == 5
        b1 = 1;
    else
        error('check lattice size')
        b1 = size(beam_model.Aero.lattice_vlm.VORTEX,2)/2;
    end
    
    p1(:,:) = beam_model.Aero.lattice_vlm.VORTEX(:,b1,:);		%Calculating panel vortex midpoint
    p2(:,:) = beam_model.Aero.lattice_vlm.VORTEX(:,b1+1,:);	%to use as a force locus
    QuarterCOLLOC(:,:) = (p1+p2)./2;	    % LOCAL control point, vortex midpoint.
    CGMomentArm = QuarterCOLLOC-ones(size(QuarterCOLLOC,1),1)*beam_model.WB.CG;
    
    WindVel = sparse(3*npanels,6*npanels);
    CrossCGArm = sparse(3*npanels,3*npanels);
    
    for n = 1:npanels
        CrossCGArm(3*(n-1)+1:3*(n-1)+3,3*(n-1)+1:3*(n-1)+3) = crossm(CGMomentArm(n,:));
        WindVel(3*(n-1)+1:3*(n-1)+3,6*(n-1)+1:6*(n-1)+3) = eye(3);
    end
    
    % Extract
    
    ExtLift = sparse(3,3*npanels);
    for n = 1:3
        ExtLift(n,n:3:end) = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    states.L2Vector   = sparse(L2Vector);
    states.SkAeroArm  = sparse(SkAeroArm);
    states.SumAeroArm = sparse(SumAeroArm);
    states.WindVel    = WindVel;
    states.CrossCGArm = CrossCGArm;
    states.ExtLift    = ExtLift;
    
    % Creates a large sparse matrix that interpolates a set of states
    % (including the aeronodes) to the aeropanels
    states.InterpCol    = sparse(IntCol);
    states.InterpFor    = sparse(IntFor);
    states.InterpSpline = sparse(IntSpline);
end

% Create a large sparse matrix that copies the displacements from the master
% nodes to the slave aerodynamic splining nodes
states.DOF2Str = sparse(DOF2STR);
states.ROT2Str = sparse(ROT2STR);
states.ROT2TRA = sparse(ROT2TRA);
% Block Diagonal of the sectional stiffness matrices
states.D         = sparse(D);

states.BARR0     = sparse(BARR0);
states.NR0       = sparse(NR0);
states.NR        = states.NR0;
states.NShape    = sparse(NShape);
states.NShapeD   = sparse(NShapeD);
states.Offsets   = Offsets;
states.NodeCoord = NodeCoord;
states.Extpitrat = Extpitrat;
%states.DOFstates = sparse([rdofs1,rdofs2],[cdofs1,cdofs2],[vdofs1,vdofs2],nbar*18,ndof);

%states.DOF       = sparse(rdof,cdof,vdof,6*ngrid,ndof); % recovers coordinates from DOF
%states.Cangle    = sparse(rangle,cangle,vangle,3*ngrid,ndof);
%states.PO        = PO;
%states.Gmat      = sparse(eye(3*2*2*nbar));
%states.C_tr      = sparse(rdoftr,cdoftr,vdoftr,ndof,ndof);
%states.C_rot     = sparse(rdofrot,cdofrot,vdofrot,ndof,ndof);

%states.DOFstatestr = sparse(rdofs1,cdofs1,vdofs1,nbar*18,ndof);
%states.DOFstatesrot = sparse(rdofs2,cdofs2,vdofs2,nbar*18,ndof);

%states.C_tr(states.C_tr>1) = 1;
%states.C_rot(states.C_rot>1) = 1;

beam_model.StateMatrices = states;

end