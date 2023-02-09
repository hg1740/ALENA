% No Control or support cards and only angle of attack trim
% 
% 
%   Author: Dario Calderon 

function [Static,beam_model] = solve_linear_trim_customisedV2(Static,beam_model,TrimIdx)

if nargin == 2
    TrimIdx = 1;
end

IndDrag = 1;

% Define the perturbation value used in the finite difference
EPS = D2R(0.001);

% Initialise the deformed lattice from the input lattice
beam_model.Aero.lattice = beam_model.Aero.lattice_vlm;

% fid for printing information
fid = beam_model.Param.FID;

% Define the centre of gravity as the reference point for the moment
% calculations
beam_model.Aero.geo.ref_point = beam_model.WB.CG;

% Trim Index to to be used to populate the output structure (change this?)
%TrimIdx = beam_model.SolParam.trimI;

beam_model.Res.SOL = 'Static linear constrained aeroelastic';
beam_model.Res.Aero = [];

% Define some useful counts
ngrid = beam_model.Info.ngrid;
nbar  = beam_model.Info.nbar;
nbeam = beam_model.Info.nbeam;
ndof  = beam_model.Info.ndof;
npa   = length(beam_model.Aero.lattice_vlm.COLLOC);
nc    = beam_model.Aero.geo.nc;

% Recover the number of
nr = [];
if nc
    nr = find(beam_model.Aero.Trim.CS.MPC == 0);
end

% Recover the number or fixed and free control surface dofs
if beam_model.Aero.geo.nc
    ncs_fixed = sum(beam_model.Res.CS.Fixed);
    ncs_free  = length(beam_model.Res.CS.Fixed) - ncs_fixed;
else
    ncs_free = 0;
end

% store bar internal forces
beam_model.Res.Bar.CForces = zeros(2, 6, nbar);
beam_model.Res.Beam.CForces = zeros(2, 6, nbeam);

% store bar internal strains and curvatures
beam_model.Res.Bar.CStrains = zeros(2, 6, nbar);
beam_model.Res.Beam.CStrains = zeros(2, 6, nbeam);

% store bar stresses
beam_model.Res.Bar.CStresses = zeros(2, 4, nbar);
beam_model.Res.Beam.CStresses = zeros(2, 4, nbeam);
beam_model.Res.Bar.CSM  = [];
beam_model.Res.Beam.CSM = [];

% store nodal displacement
beam_model.Res.NDispl = zeros(ngrid, 6);
NODEPOS = beam_model.Node.Coord;

% store updated bar rotations
beam_model.Res.Bar.R = beam_model.Bar.R;
beam_model.Res.Bar.Colloc = beam_model.Bar.Colloc;

beam_model.Res.Beam.R = beam_model.Beam.R;
beam_model.Res.Beam.Colloc = beam_model.Beam.Colloc;

% store updated node rotation
beam_model.Res.NRd = beam_model.Node.R;

% Fudge Variables to run
beam_model.Celas.ID = [];
%beam_model.Info.nbeam = 0;

%beam_model.Info.nbaero = 0;
beam_model.Param.FUSE_DP = [];

% Assemble STIFFNESS matrix
K = st_lin_matrix(beam_model.Info, beam_model.Node.DOF, beam_model.Node.R, beam_model.Node.Coord, beam_model.Bar, beam_model.Beam, beam_model.Celas,beam_model.CBush);

if ~isempty(beam_model.RBE2.ID)
    K = RBE2Assembly(beam_model.RBE2,K);
end

[Kll,ldof] = get_free_stiff(K, beam_model.Node, beam_model.Param.SUPORT, beam_model.Param.EPS);
rdof       = [];

% Assemble MASS matrix
M = ms_matrix(beam_model.Info, beam_model.Node.DOF, beam_model.Node.R, beam_model.ConM, beam_model.Bar, beam_model.Beam);

% Assemble Deformabilitty Influence Matrix
if (beam_model.Info.spline_type == 1)
    [Qaa, CPaeroDef, ~] = DeformCoefficient(beam_model.Aero.Interp.Ic, beam_model.Aero.Interp.Imv, beam_model.Aero.geo, beam_model.Aero.lattice, beam_model.Aero.state);
    CPaeroMDef = beam_model.StateMatrices.CrossCGArm*CPaeroDef;
else
    [Qaa, CPaeroDef, CPaeroMDef] = VLM_Gradient(beam_model.Aero.state,beam_model.StateMatrices,beam_model,beam_model.Aero.lattice);
end

% Define the translation acceleration vector for the rigid body
% ACC    = zeros(6,1);
% ACC(1) = 0;
% ACC(3) = - (beam_model.Param.G);

ACC = -beam_model.Aero.Trim.FM.Value(7:12);


% Determine the inertial loads
Fi = gf_iner_nodal(ndof, beam_model.Node.DOF, M, ACC);

Precurvature = 0;

if Precurvature == 1
    
    load('Precurvature.mat')
    
    F_pre = ...
        get_precurvature_load(beam_model.Info.nbeam, beam_model.Beam, beam_model.Node, ...
        beam_model.Res.NDispl, beam_model.Info.ngrid,Precurvature);
    
    F_curv = get_internal_forces(ndof, beam_model.Node.DOF, F_pre);
    
end
% store stability derivatives
CREF = 1;

BREF = beam_model.Aero.ref.b_ref;
SREF = beam_model.Aero.ref.S_ref;
VREF = beam_model.Aero.state.AS;

RHOREF = beam_model.Aero.state.rho;
QINF   = 0.5 * RHOREF * VREF^2;
QINFS  = QINF * SREF;

% ONLY RIGID TRIM
beam_model.Res.Aero.RStab_Der = [];
beam_model.Res.Aero.RIntercept = [];
beam_model.Res.Aero.RTrim_sol = [];

%% CASE 1: Static Rigid Load at null condition
dummy_aero = beam_model.Aero;
dummy_aero.state.alpha = 0; dummy_aero.state.betha = 0;
dummy_aero.state.P = 0;     dummy_aero.state.Q = 0;      dummy_aero.state.R = 0;

beam_model.Aero.AIC = [];

% If the twist or camber is corrected then update the lattice to account
% for this. This should not be done too often ... and it can be costly
if beam_model.twist_corr == 1
    
    if beam_model.Info.spline_type == 1
        lattice_defo = update_vlm_mesh1(beam_model.Node, zeros(ndof,1), beam_model.Aero, beam_model);
    else
        lattice_defo = update_vlm_mesh(beam_model.Node, NODEPOS, ...
            beam_model.Node.R, beam_model.Aero, beam_model);
    end
else
    lattice_defo = beam_model.Aero.lattice;
end

% Determine the global rotation matrix due to new angles of attack
GlobalRotationMatrix = [cos(dummy_aero.state.betha)*cos(dummy_aero.state.alpha),        -sin(dummy_aero.state.betha),          cos(dummy_aero.state.betha)*sin(dummy_aero.state.alpha) ;...
    cos(dummy_aero.state.alpha)*sin(dummy_aero.state.betha),         cos(dummy_aero.state.betha),          sin(dummy_aero.state.betha)*sin(dummy_aero.state.alpha) ;...
    -sin(dummy_aero.state.alpha),                                              0,                           cos(dummy_aero.state.alpha)]';

% Determine the new flow vector
FlowVector = (beam_model.Aero.state.AS*GlobalRotationMatrix*[1,0,0]')';

% Run the vlm for the new parameters
results0 = VLM_HS(FlowVector,dummy_aero.state.rho,lattice_defo,beam_model.Aero.state.SIMXZ,[],[],[],0,[]);

% Extract some useful forces
results0 = Aero_process(results0,lattice_defo,beam_model.Aero.state.rho,FlowVector,GlobalRotationMatrix,beam_model.Aero.ref.S_ref,beam_model.Aero.geo,beam_model.Res.WB.CG);

% Update the results.F as they are the local values ... needed for that
% static analysis
results0.F = results0.Fglobal;

% Transfer the results.F over to the aerodynamic mesh
Fa0 = gf_transfer_aero_nodal(beam_model.Info, beam_model.Node.DOF, beam_model.Node, dummy_aero, results0);

%   Total external forces
if Precurvature == 1
    F = Fi + Fa0 - F_curv;
else
    F = Fi + Fa0;
end
if ~isempty(beam_model.RBE2.ID)
    F = RBE2Assembly2(beam_model.RBE2,F);
end

% Setup the Fa_state_aero variable - aerodynamic derivatives
Fa_State_aero = zeros(npa*3, 5);

%% CASE 2: Rigid Body Variations
% - Alpha
state_count = 1;
dummy_aero  = beam_model.Aero;
dummy_aero.state.alpha = EPS;
dummy_aero.state.betha = 0;
dummy_aero.state.P     = 0;
dummy_aero.state.Q     = 0;
dummy_aero.state.R     = 0;

% NOTE: lattice_defo is not recalculated here ... no need!

% Determine the global rotation matrix due to new angle of attack
GlobalRotationMatrix = [cos(dummy_aero.state.betha)*cos(dummy_aero.state.alpha),        -sin(dummy_aero.state.betha),          cos(dummy_aero.state.betha)*sin(dummy_aero.state.alpha) ;...
    cos(dummy_aero.state.alpha)*sin(dummy_aero.state.betha),         cos(dummy_aero.state.betha),          sin(dummy_aero.state.betha)*sin(dummy_aero.state.alpha) ;...
    -sin(dummy_aero.state.alpha),                                              0,                           cos(dummy_aero.state.alpha)]';

% Determine the new flow vector
FlowVector = (beam_model.Aero.state.AS*GlobalRotationMatrix*[1,0,0]')';

% Run the vlm for the new parameters
results = VLM_HS(FlowVector,dummy_aero.state.rho,lattice_defo,beam_model.Aero.state.SIMXZ,[],results0.AIC,[],0,results0.IndDownwash);

% Extract some useful forces
results = Aero_process(results,lattice_defo,beam_model.Aero.state.rho,FlowVector,GlobalRotationMatrix,beam_model.Aero.ref.S_ref,beam_model.Aero.geo,beam_model.Res.WB.CG);

% Calculate the derivative of the aerodynamic forces on the structural
% nodes to perturbations of the angle of attack

% Update the results.F as they are the local values ... needed for that
% static analysis
results.F = results.Fglobal;

Fa_ALPHA = gf_transfer_aero_nodal(beam_model.Info, beam_model.Node.DOF, beam_model.Node, dummy_aero, results);
Fa_ALPHA = (Fa_ALPHA - Fa0)./ EPS; % get force variation

Fa_State_aero([1:3:end], state_count) = [results.F(:,1)-results0.F(:,1)]./ EPS;
Fa_State_aero([2:3:end], state_count) = [results.F(:,2)-results0.F(:,2)]./ EPS;
Fa_State_aero([3:3:end], state_count) = [results.F(:,3)-results0.F(:,3)]./ EPS;

Fa_Rigid_aero(:,state_count) = (results.FORCES - results0.FORCES)'/EPS;
Fa_RigidM_aero(:,state_count) = (results.MOMENT - results0.MOMENT)'/EPS;
%%
% - Sideslip
Fa_BETA = zeros(ndof,1);
Fa_P    = zeros(ndof,1);
Fa_Q    = zeros(ndof,1);
Fa_R    = zeros(ndof,1);

%% CASE 3: Control Surface Variations
dummy_aero = beam_model.Aero;
dummy_aero.state.alpha = 0;
dummy_aero.state.betha = 0;
dummy_aero.state.P     = 0;
dummy_aero.state.Q     = 0;
dummy_aero.state.R     = 0;

FCl = [];FCr = [];
FCtot = [];FCtot_aero = [];FCtot_aeroM = [];FCtot_Rigid_aero = [];FCtot_RigidM_aero = [];

if (~isempty(nc))
    %print_controls(outp, nr, beam_model.Aero.lattice_vlm.Control.Name, beam_model.Res.CS);
    FCtot = zeros(ndof, length(nr));
    FCtot_aero = zeros(npa*3, length(nr));
    FCtot_aeroM = zeros(npa*3, length(nr));
    FCl = zeros(length(ldof), length(nr));
    FCr = zeros(length(rdof), length(nr));
    %
    LUMP_DOF = [];
    LUMP_COEFF = [];
    %
    for k=1:length(nr)
        
        % erase all rotations
        beam_model.Res.CS.Value(1:length(beam_model.Aero.Trim.CS.MPC)) = 0.0;
        % set master rotation
        %         EPS = pi*20/180;
        beam_model.Res.CS.Value(nr(k)) = EPS;
        % look for slave
        cdof = find(beam_model.Aero.Trim.CS.MPC == nr(k));
        LUMP_DOF(k).data = [nr(k)];
        LUMP_COEFF(k).data = [1];
        if (~isempty(cdof))
            beam_model.Res.CS.Value(cdof) = EPS*beam_model.Aero.Trim.CS.Coeff(cdof);
            LUMP_DOF(k).data = [nr(k), cdof];
            LUMP_COEFF(k).data = [1, beam_model.Aero.Trim.CS.Coeff(cdof)];
        end
        
        NODEPOS = beam_model.Node.Coord;
        
        lattice_defo_control = rotate_control_norm(beam_model.Aero.ref, dummy_aero.state, beam_model.Aero.geo, ...
            lattice_defo, beam_model.Res.CS.Value, ...
            beam_model.Aero.lattice.Control.Hinge);
        
        
        % Zero rigid body matrices simplify the rotation matrix and the
        % inflow vector
        GlobalRotationMatrix = eye(3);
        
        FlowVector = beam_model.Aero.state.AS*[1,0,0];
        
        results = VLM_HS(FlowVector,dummy_aero.state.rho,lattice_defo_control,beam_model.Aero.state.SIMXZ,[],results0.AIC,[],0,results0.IndDownwash);
        results = Aero_process(results,lattice_defo_control,beam_model.Aero.state.rho,FlowVector,GlobalRotationMatrix,beam_model.Aero.ref.S_ref,beam_model.Aero.geo,beam_model.Res.WB.CG);
        
        % Update the results.F as they are the local values ... needed for that
        % static analysis
        results.F = results.Fglobal;
        
        Fa_C = gf_transfer_aero_nodal(beam_model.Info, beam_model.Node.DOF, beam_model.Node, dummy_aero, ...
            results);
        
        Fa_C = (Fa_C - Fa0) ./ EPS;
        FCtot_aero([1:3:end], k) = [results.F(:,1)-results0.F(:,1)]./ EPS;
        FCtot_aero([2:3:end], k) = [results.F(:,2)-results0.F(:,2)]./ EPS;
        FCtot_aero([3:3:end], k) = [results.F(:,3)-results0.F(:,3)]./ EPS;
        
        FCtot_aeroM([1:3:end], k) = [results.M(:,1)-results0.M(:,1)]./ EPS;
        FCtot_aeroM([2:3:end], k) = [results.M(:,2)-results0.M(:,2)]./ EPS;
        FCtot_aeroM([3:3:end], k) = [results.M(:,3)-results0.M(:,3)]./ EPS;
        
        FCl(:,k) = Fa_C(ldof,1);
        FCr(:,k) = Fa_C(rdof,1);
        FCtot(:,k) = Fa_C;
        FCtot_Rigid_aero(:,k) = (results.FORCES - results0.FORCES)'/EPS;
        FCtot_RigidM_aero(:,k) = (results.MOMENT - results0.MOMENT)'/EPS;
    end
    beam_model.Res.CS.Value(1:length(beam_model.Aero.Trim.CS.MPC)) = 0.0;
    
end

%%
Kax    = -[Fa_ALPHA,Fa_BETA,Fa_P,Fa_Q,Fa_R,FCtot];

KaxDOF = Kax;
if ~isempty(beam_model.RBE2.ID)
    Kax = RBE2Assembly2(beam_model.RBE2,Kax);
end

Kaxl = Kax(ldof,:);

QaaDOF = Qaa;
if ~isempty(beam_model.RBE2.ID)
    Qaa = RBE2Assembly(beam_model.RBE2,Qaa);
end
Qaall = Qaa(ldof, ldof);

beam_model.Res.NDispl = zeros(ngrid, 6);
beam_model.Res.NRd    = beam_model.Node.R;

%% CHECK THIS!!!
acc_dof = beam_model.Aero.Trim.FM.Fixed(7:12);
indexa = find(acc_dof); % fixed acc dofs
index_acc = setdiff([1:6], indexa); % free acc dofs
fm_dof = beam_model.Aero.Trim.FM.Fixed(2:6);
indexfm = find(fm_dof); % fixed fm dofs
index_fm = setdiff([1:5], indexfm); % free fm dofs
index_cs = [];
indexcs = [];
if nc
    nr = find(beam_model.Aero.Trim.CS.MPC == 0); % get master rotations
end
if (~isempty(nr))
    cs_dof = beam_model.Aero.Trim.CS.Fixed(nr); % check if master deflections are fixed
    indexcs = find(cs_dof); % fixed cs dofs
    index_cs = setdiff([1:length(nr)], indexcs); % free cs dofs
end

%%
beam_model.Res.Struct.M = M;
beam_model.Res.Struct.K = K;

if ~isempty(beam_model.RBE2)
    Fa02 = RBE2Assembly2(beam_model.RBE2,Fa0);
else
    Fa02 = Fa0;
end

beam_model.Res.Struct.F = F - Fa02;

beam_model.Res.Struct.D = [];

beam_model.Res.Struct.ldof = ldof;
beam_model.Res.Struct.rdof = rdof;

% aero
beam_model.Res.Aero.Qaa    = Qaa;
beam_model.Res.Aero.Kax    = Kax;
beam_model.Res.Aero.QaaDOF = QaaDOF;
beam_model.Res.Aero.KaxDOF = KaxDOF;
beam_model.Res.Aero.Fa0    = Fa02;
beam_model.Res.Aero.Fa0DOF = Fa0;

% matrices along aero panels
beam_model.Res.CPaero.State         = Fa_State_aero;
beam_model.Res.CPaero.Control       = FCtot_aero;
beam_model.Res.CPaero.Defo          = CPaeroDef;
beam_model.Res.CPaero.DefoM         = CPaeroMDef;
beam_model.Res.CPaero.F0            = zeros(npa*3,1);
beam_model.Res.CPaero.F0(1:3:end,1) = results0.F(:,1);
beam_model.Res.CPaero.F0(2:3:end,1) = results0.F(:,2);
beam_model.Res.CPaero.F0(3:3:end,1) = results0.F(:,3);

%%
Kall    = Kll - Qaall;
invKall = inv(Kall);

%%
ALX   = invKall * Kaxl;
UINTL = invKall * F(ldof,1);
UDD(:,1) = ACC;

UCS = [];
if (~isempty(nr))
    defl = beam_model.Aero.Trim.CS.Value(nr)';
    if (~isempty(indexcs))
        UCS(:,1) = defl(indexcs);
    end
    free_cs = intersect(nr,  find(beam_model.Res.CS.Fixed==0));
    NAME = beam_model.Aero.lattice.Control.Name;
end

ADOF = [];
aero_data = [];
for k=1:length(beam_model.Aero.ID)
    n1 = beam_model.Aero.lattice_vlm.DOF(k,1,1);
    n2 = beam_model.Aero.lattice_vlm.DOF(k,1,2);
    ADOF = [ADOF, [n1:n2]];
    for j=n1:n2
        coord = (beam_model.Aero.lattice_vlm.VORTEX([j],[4],:) + beam_model.Aero.lattice_vlm.VORTEX([j],[5],:)).*0.5;
        aero_data = [aero_data; [coord(:,:,1), coord(:,:,2), coord(:,:,3)]];
    end
end

beam_model.Res.Aero.RStab_Der.Control.Name = {};
beam_model.Res.Aero.DStab_Der.Control.Name = {};

% Solve for twist first!
if (~isempty(nr))
    beam_model.Res.Aero.RStab_Der.Control.Name = beam_model.Aero.lattice_vlm.Control.Name(nr);
    beam_model.Res.Aero.DStab_Der.Control.Name = beam_model.Aero.lattice_vlm.Control.Name(nr);
end

%% Generate the Trim Jacobian Matrix
% I need to fix this for multiple RBE2s
% Fa_Alpha_RBE = -Kaxl(:,1);
% Fc_tot_RBE = -Kaxl(:,6:end);

dF_dalpha = Fa_Rigid_aero(3) + beam_model.StateMatrices.ExtLift(3,:)*CPaeroDef*invKall*Fa_ALPHA;
dM_dalpha = Fa_RigidM_aero(2) + beam_model.StateMatrices.ExtLift(2,:)*CPaeroMDef*invKall*Fa_ALPHA;

if ncs_free
    dF_dbeta  = FCtot_Rigid_aero(3,index_cs) + beam_model.StateMatrices.ExtLift(3,:)*CPaeroDef*invKall*FCtot(:,index_cs);
    dM_dbeta  = FCtot_RigidM_aero(2,index_cs) + beam_model.StateMatrices.ExtLift(2,:)*CPaeroMDef*invKall*FCtot(:,index_cs);
    Trim_J_R = [Fa_Rigid_aero(3),FCtot_Rigid_aero(3,index_cs);Fa_RigidM_aero(2),FCtot_RigidM_aero(2,index_cs)];
    Trim_J_F = [dF_dalpha,dF_dbeta;dM_dalpha,dM_dbeta];
else
    Trim_J_R = Fa_Rigid_aero(3);
    Trim_J_F = dF_dalpha;
end

if beam_model.SolParam.Trim == 1
    
    %% Solve for the trim variables
    
    % I need to fix this so that it takes into account fixed control surfaces!!
    % Need to calculate the difference in the required Rigid Forces (All 6
    % DOF?)
    
    % Calculate the required change in trim forces
    if (~isempty(nr))
        if isempty(indexcs)
            DeltaLiftR = -ACC(3)*(beam_model.WB.MCG(1,1)) - (results0.FORCES(3));
            DeltaLiftF = -ACC(3)*(beam_model.WB.MCG(1,1)) - (results0.FORCES(3) + beam_model.StateMatrices.ExtLift(3,:)*CPaeroDef*invKall*Fa0);
            
            DeltaMomentR = -(results0.MOMENT(2));
            DeltaMomentF = -(results0.MOMENT(2)+ beam_model.StateMatrices.ExtLift(2,:)*CPaeroMDef*invKall*Fa0);
        else
            DeltaLiftR = -ACC(3)*(beam_model.WB.MCG(1,1)) - (results0.FORCES(3) + sum(FCtot_Rigid_aero(3,indexcs).*UCS'));
            DeltaLiftF = -ACC(3)*(beam_model.WB.MCG(1,1)) - (results0.FORCES(3) + sum(FCtot_Rigid_aero(3,indexcs).*UCS') + beam_model.StateMatrices.ExtLift(3,:)*CPaeroDef*invKall*(Fa0+sum(FCl(:,indexcs)*UCS,2)));
            
            DeltaMomentR = -(results0.MOMENT(2) + sum(FCtot_RigidM_aero(2,indexcs).*UCS'));
            DeltaMomentF = -(results0.MOMENT(2) + sum(FCtot_RigidM_aero(2,indexcs).*UCS' + beam_model.StateMatrices.ExtLift(2,:)*CPaeroMDef*invKall*(Fa0+sum(FCl(:,indexcs)*UCS,2)) ));
            
        end
    else
        DeltaLiftR = -ACC(3)*(beam_model.WB.MCG(1,1)) - (results0.FORCES(3));
        DeltaLiftF = -ACC(3)*(beam_model.WB.MCG(1,1)) - (results0.FORCES(3) + beam_model.StateMatrices.ExtLift(3,:)*CPaeroDef*invKall*Fa0);
    end
    
    if ncs_free
        if size(Trim_J_R,2)>2
            Trim_SOL_R = pinv(Trim_J_R)*[DeltaLiftR;DeltaMomentR]; % Add initial moment value here!
            Trim_SOL_F = pinv(Trim_J_F)*[DeltaLiftF;DeltaMomentF];
        else
            Trim_SOL_R = Trim_J_R\[DeltaLiftR;DeltaMomentR]; % Add initial moment value here!
            Trim_SOL_F = Trim_J_F\[DeltaLiftF;DeltaMomentF];
        end
    else
        Trim_SOL_R = DeltaLiftR/Trim_J_R;
        Trim_SOL_F = DeltaLiftF/Trim_J_F;
    end
%     
%     fprintf('\n\n\t... Trim Angle for Rigid Wing:         %f deg.',   180*Trim_SOL_R(1)/pi);
%     if ncs_free
%         fprintf('\n\t... Elevator Angle for Rigid Wing:     %f deg.',   180*Trim_SOL_R(2)/pi);
%     end
%     fprintf('\n\n\t... Trim Angle for Flexible Wing:      %f deg.',   180*Trim_SOL_F(1)/pi);
%     if ncs_free
%         fprintf('\n\t... Elevator Angle for Flexible Wing:  %f deg.',   180*Trim_SOL_F(2)/pi);
%     end
else
    
    Trim_SOL_F(1) = pi*beam_model.SolParam.WingAngle(1)/180;
    Trim_SOL_R(1) = pi*beam_model.SolParam.WingAngle(1)/180;
    
    %fprintf('\n\n\t... Fixed Angle for Wing:         %f deg.',   180*Trim_SOL_R(1)/pi);
    
end

if nc
    ControlDefR = zeros(1,length(nr));
    if ncs_free
        ControlDefR(index_cs) = Trim_SOL_R(2:end);
    end
    if ~isempty(indexcs)
        ControlDefR(indexcs) = UCS';
    end
    ControlDefF = zeros(1,length(nr));
    if ncs_free
        ControlDefF(index_cs) = Trim_SOL_F(2:end);
    end
    if ~isempty(indexcs)
        ControlDefF(indexcs) = UCS';
    end
end

% I need to fix this so that it takes into account fixed control surfaces!!
if nc
    UX(:,1) = [Trim_SOL_R(1),0,0,0,0,ControlDefR];
    UX(:,2) = [Trim_SOL_F(1),0,0,0,0,ControlDefF];
else
    UX(:,1) = [Trim_SOL_R(1),0,0,0,0];
    UX(:,2) = [Trim_SOL_F(1),0,0,0,0];
end

UR = [];

UL = - ALX * UX(:,2) + UINTL;

SOL =  zeros(size(K,1),1); SOL(rdof,1) = UR; SOL(ldof,1) = UL;
if ~isempty(beam_model.RBE2.ID)
    SOL = RBE2disp(beam_model.RBE2,SOL,ndof);
end

%Fext = F(ldof,1) - Kaxl * UX(:,2);

Fext = -Kaxl * UX(:,2);

%%%% REMOVE THE rdof from the SOL!!!!!!!!!
gdef = zeros(beam_model.Info.ngrid, 6);
Fgrid = zeros(beam_model.Info.ngrid, 6);
for n = 1:beam_model.Info.ngrid
    dof = beam_model.Node.DOF(n, 1:6);
    index = find(dof);
    if ~isempty(index)
        gdef(n, index) = SOL(dof(index));
        Fgrid(n, index) = Fext(dof(index));
    end
end

% store nodal displacement
beam_model.Res.NDispl(:,:,1) = gdef;

% set delta Rot
for n = 1:beam_model.Info.ngrid
    beam_model.Res.NRd(:,:,n,1) = Rmat(gdef(n, 4:6)');
end

[DS, DR] = get_nodal_displ(beam_model.Info.ngrid, beam_model.Node.DOF, SOL); %

% update BAR rotations
beam_model.Res.Bar.R = update_bar_rot(nbar, DR, beam_model.Bar.Conn, beam_model.Res.Bar.R , DS); % ok

% update BEAM rotations
beam_model.Res.Beam.R = update_bar_rot(nbeam, DR, beam_model.Beam.Conn, beam_model.Res.Beam.R , DS);
%
COORD = beam_model.Node.Coord + beam_model.Res.NDispl(:,1:3);
beam_model.Res.Bar.Colloc = bar_defo_colloc(nbar, beam_model.Bar, beam_model.Node.DOF, COORD, beam_model.Res.NRd);
beam_model.Res.Beam.Colloc = bar_defo_colloc(nbeam, beam_model.Beam, beam_model.Node.DOF, COORD, beam_model.Res.NRd);
beam_model.Res.WB = [];

% calculate new mass matrix
[beam_model.Res.WB.CG, beam_model.Res.WB.MCG, beam_model.Res.WB.MRP] = ...
    wb_set_conm_mass(beam_model.Info.nconm, beam_model.Node.Index, COORD, beam_model.Res.NRd, beam_model.Param.GRDPNT, beam_model.ConM);

% set bar mass CG
[beam_model.Res.WB.CG, beam_model.Res.WB.MCG, beam_model.Res.WB.MRP] =...
    wb_add_bar_mass(nbar, COORD, beam_model.Res.NRd, beam_model.Res.WB.CG, beam_model.Param.GRDPNT, beam_model.Res.WB.MCG, beam_model.Res.WB.MRP, beam_model.Bar);

% set beam mass CG
[beam_model.Res.WB.CG, beam_model.Res.WB.MCG, beam_model.Res.WB.MRP] =...
    wb_add_bar_mass(nbeam, COORD, beam_model.Res.NRd, beam_model.Res.WB.CG, beam_model.Param.GRDPNT, beam_model.Res.WB.MCG, beam_model.Res.WB.MRP, beam_model.Beam);

% get principal axes
%[beam_model.Res.WB.MCG_pa, beam_model.Res.WB.R_pa] = wb_principal_axis(beam_model.Res.WB.MCG);

%
% update aerobeam nodes (if any)
%
if (beam_model.Info.nrbe0 > 0)
    AERO_POS = update_aerobeam_node(beam_model.Info.ngrid, beam_model.Node, beam_model.Res.NDispl(:,1:3,1),...
        beam_model.Res.NRd(:,:,:,1)-repmat(eye(3,3),[1,1,beam_model.Info.ngrid]));
    % update coord database with slave nodes position
    for n=1:beam_model.Info.ngrid
        ne = length(beam_model.Node.Aero.Index(n).data);
        if ne
            beam_model.Res.NDispl(beam_model.Node.Aero.Index(n).data, 1:3, 1) = AERO_POS(n).data';
        end
    end
    clear AERO_POS;
end

[beam_model.Res.Bar.CForces, beam_model.Res.Bar.CStrains, beam_model.Res.Bar.CStresses, beam_model.Res.Bar.CSM] = ...
    get_bar_force_strain(beam_model.Info.nbar, beam_model.Bar, beam_model.PBar, beam_model.Mat, beam_model.Node, ...
    beam_model.Res.NDispl, beam_model.Param.FUSE_DP);

[beam_model.Res.Beam.CForces, beam_model.Res.Beam.CStrains, beam_model.Res.Beam.CStresses, beam_model.Res.Beam.CSM] = ...
    get_bar_force_strain(beam_model.Info.nbeam, beam_model.Beam, beam_model.PBeam, beam_model.Mat, beam_model.Node, ...
    beam_model.Res.NDispl, beam_model.Param.FUSE_DP);

BARR = beam_model.Res.Bar.R;
BARF =  store_global_elem_force(nbar, BARR, beam_model.Res.Bar.CForces);

BEAMR = beam_model.Res.Beam.R;
BEAMF =  store_global_elem_force(nbeam, BEAMR, beam_model.Res.Beam.CForces);

%beam_model.Aero.state.alpha=(pi*angle/180);

beam_model.Res.Aero.XYZ = aero_data;

% Store Rigid Loads

beam_model.Aero.lattice = beam_model.Aero.lattice_vlm;

if nc
    % Update the Master rotations
    beam_model.Res.CS.Value(nr) = ControlDefF;
    
    index = find(beam_model.Aero.Trim.CS.MPC);
    for k=1:length(index)
        beam_model.Res.CS.Value(1,index(k)) = beam_model.Aero.Trim.CS.Coeff(index(k)) * beam_model.Res.CS.Value(1,beam_model.Aero.Trim.CS.MPC(index(k)));
    end
    
    beam_model.Aero.lattice = update_hinge_lines(beam_model.Aero.geo, ...
        beam_model.Aero.lattice);
    
    beam_model.Aero.lattice = rotate_control_norm(beam_model.Aero.ref, beam_model.Aero.state, beam_model.Aero.geo, ...
        beam_model.Aero.lattice, beam_model.Res.CS.Value, ...
        beam_model.Aero.lattice.Control.Hinge);
end

[MAC, MAC_LE, MAC_AC] = CalculateMAC(beam_model.Aero);

% Store Flexible Loads

beam_model.Res.Aero.DStab_Der.Alpha.Cm = dM_dalpha/(0.5*beam_model.Aero.state.rho*beam_model.Aero.state.AS^2*beam_model.Aero.ref.S_ref*MAC);

beam_model.Res.Aero.DStab_Der.Alpha.Cl = dF_dalpha/(0.5*beam_model.Aero.state.rho*beam_model.Aero.state.AS^2*beam_model.Aero.ref.S_ref);

beam_model.Res.Aero.DTrim_sol        = store_trim_sol(fid, UDD(:,1), UX(:,2), beam_model.Res.Aero.DStab_Der.Control.Name);
beam_model.Res.Aero.F0_DTrim         = recover_trim_vlm_forces_defo(beam_model.Res.CPaero.F0, beam_model.Res.Aero.DTrim_sol, beam_model.Res.CPaero, SOL, CREF, BREF, SREF, VREF, ADOF);

% % Iterate through Parts
% if isfield(beam_model,'PartList')
%     fprintf('\n\n');
%     for i = 1:length(beam_model.PartList)
%         
%         PartName = beam_model.PartList(i).Name;
%         
%         PartLift = sum(beam_model.Res.Aero.F0_DTrim(beam_model.PartList(i).AeroIdx,3));
%         
%         fprintf('\t\t -%s: Lift = %8.4f N\n',PartName,PartLift);
%         
%     end
%     
% end
%beam_model.Res.Aero.Fgrid            = Fa_ALPHA* + FCl(:,k)
beam_model.Res.Aero.results.F        = [];
beam_model.Res.Aero.results.F(:,1) = beam_model.Res.Aero.F0_DTrim(:,1);
beam_model.Res.Aero.results.F(:,2) = beam_model.Res.Aero.F0_DTrim(:,2);
beam_model.Res.Aero.results.F(:,3) = beam_model.Res.Aero.F0_DTrim(:,3);
beam_model.Res.Aero.results  = Aero_process(beam_model.Res.Aero.results,beam_model.Aero.lattice,beam_model.Aero.state.rho,[beam_model.Aero.state.AS,0,0],eye(3),beam_model.Aero.ref.S_ref,beam_model.Aero.geo,beam_model.Res.WB.CG);

NR = beam_model.Res.NRd;
%BARR = update_bar_rot(beam_model.Info.nbar, DR, beam_model.Bar.Conn, BARR, gdef);

%% Calculate Drag for Flexible Wing

NODEPOS = beam_model.Node.Coord + beam_model.Res.NDispl(:,1:3);

% Update the deflection and unit normals for the flexible solution
if beam_model.Info.spline_type == 1
    lattice_defo = update_vlm_mesh1(beam_model.Node, SOL, beam_model.Aero,beam_model);
else
    lattice_defo = update_vlm_mesh(beam_model.Node, NODEPOS, ...
        beam_model.Res.NRd, beam_model.Aero,beam_model);
end

if nc
    lattice_defo = rotate_control_norm(beam_model.Aero.ref, beam_model.Aero.state, beam_model.Aero.geo, ...
        lattice_defo, beam_model.Res.CS.Value, ...
        lattice_defo.Control.Hinge);
end

% Re-initialise the case:
dummy_aero.state.alpha = Trim_SOL_F(1);
dummy_aero.state.betha = 0;
dummy_aero.state.P     = 0;
dummy_aero.state.Q     = 0;
dummy_aero.state.R     = 0;
dummy_aero.lattice_vlm.N      = lattice_defo.N;
%dummy_aero.lattice_vlm.N(:,2) = 0.*dummy_aero.lattice_vlm.N(:,2);
dummy_aero.lattice_vlm.N(:,1) = dummy_aero.lattice_vlm.N(:,1)./sqrt(sum(dummy_aero.lattice.N.^2,2));
dummy_aero.lattice_vlm.N(:,2) = dummy_aero.lattice_vlm.N(:,2)./sqrt(sum(dummy_aero.lattice.N.^2,2));
dummy_aero.lattice_vlm.N(:,3) = dummy_aero.lattice_vlm.N(:,3)./sqrt(sum(dummy_aero.lattice.N.^2,2));

% Zero rigid body matrices simplify the rotation matrix and the
% inflow vector
% Determine the global rotation matrix due to new angle of attack
GlobalRotationMatrix = [cos(dummy_aero.state.betha)*cos(dummy_aero.state.alpha),        -sin(dummy_aero.state.betha),          cos(dummy_aero.state.betha)*sin(dummy_aero.state.alpha) ;...
    cos(dummy_aero.state.alpha)*sin(dummy_aero.state.betha),         cos(dummy_aero.state.betha),          sin(dummy_aero.state.betha)*sin(dummy_aero.state.alpha) ;...
    -sin(dummy_aero.state.alpha),                                              0,                           cos(dummy_aero.state.alpha)]';

% Determine the new flow vector
FlowVector = (beam_model.Aero.state.AS*GlobalRotationMatrix*[1,0,0]')';

results = VLM_HS(FlowVector,dummy_aero.state.rho,dummy_aero.lattice_vlm,beam_model.Aero.state.SIMXZ,[],[],[],1,[]);

results = Aero_process(results,dummy_aero.lattice_vlm,beam_model.Aero.state.rho,FlowVector,GlobalRotationMatrix,beam_model.Aero.ref.S_ref,beam_model.Aero.geo,beam_model.Res.WB.CG);

% Update the results.F as they are the local values ... needed for that
% static analysis
results.F = results.Fglobal;


% %[results,beam_model] = solver9(dummy_aero.state, beam_model.Aero.geo,dummy_aero.lattice_vlm,beam_model);
% [results,beam_model] = solver9(dummy_aero.state, beam_model.Aero.geo,lattice_defo,beam_model);
% results = Aero_coeff(results,beam_model.Aero.lattice_vlm,beam_model.Aero.lattice,dummy_aero.state,beam_model.Aero.ref,beam_model.Aero.geo,beam_model.Linear,beam_model.SolParam.AeroLinear);

snapnow;
beam_model.Res.Aero.results.gamma      = results.gamma;
beam_model.Res.Aero.results.gammapanel = results.gammapanel;
beam_model.Res.Aero.results.LPost      = results.L;
beam_model.Res.Aero.results.CLPost     = results.CL;
beam_model.Res.Aero.results.CD         = results.CD;
beam_model.Res.Aero.results.D          = results.D;
beam_model.Res.Aero.results.spanside   = results.spanside;
beam_model.Res.Aero.results.spancd     = results.spancd;
%beam_model.Res.Aero.results.CD = results.Trefftz_drag_Coeff;
%beam_model.Res.Aero.results.treffts = results.treffts;
beam_model.Res.SOL = SOL;

beam_model.Res.state.alpha = Trim_SOL_F(1);

%results = TrefftzDrag(results,beam_model.Aero.state,dummy_aero.lattice_vlm,beam_model.Aero.ref);

beam_model.Aero.lattice = lattice_defo;
%%
beam_model.Lattice.flexible{TrimIdx} = beam_model.Aero.lattice;
%beam_model.Loads.flexible{TRIM_INDEX}   = storeloads(beam_model);
beam_model.Results.flexible{TrimIdx} = beam_model.Res;

%[~,beam_model.Res.Beam.CForces] = BeamRotationMatrix(beam_model);

if ~isempty(beam_model.Res.Bar.CForces)
    IntLoads = storeinternalloads(beam_model.Res.Bar.CForces);
    BarForces.Fx(:,1) = IntLoads.Fx;BarForces.Fy(:,1) = IntLoads.Fy;BarForces.Fz(:,1) = IntLoads.Fz;
    BarForces.Mx(:,1) = IntLoads.Mx;BarForces.My(:,1) = IntLoads.My;BarForces.Mz(:,1) = IntLoads.Mz;
else
    IntLoads = storeinternalloads(beam_model.Res.Beam.CForces);
    BarForces.Fx(:,1) = IntLoads.Fx;BarForces.Fy(:,1) = IntLoads.Fy;BarForces.Fz(:,1) = IntLoads.Fz;
    BarForces.Mx(:,1) = IntLoads.Mx;BarForces.My(:,1) = IntLoads.My;BarForces.Mz(:,1) = IntLoads.Mz;
end

if ~isempty(beam_model.Res.Bar.CStrains)
    IntLoads = storeinternalloads(beam_model.Res.Bar.CStrains);
    BarStrains.epsx(:,1) = IntLoads.Fx;BarStrains.epsx(:,1) = IntLoads.Fy;BarStrains.epsx(:,1) = IntLoads.Fz;
    BarStrains.kapx(:,1) = IntLoads.Mx;BarStrains.kapy(:,1) = IntLoads.My;BarStrains.kapz(:,1) = IntLoads.Mz;
else
    IntLoads = storeinternalloads(beam_model.Res.Beam.CStrains);
    BarStrains.epsx(:,1) = IntLoads.Fx;BarStrains.epsy(:,1) = IntLoads.Fy;BarStrains.epsz(:,1) = IntLoads.Fz;
    BarStrains.kapx(:,1) = IntLoads.Mx;BarStrains.kapy(:,1) = IntLoads.My;BarStrains.kapz(:,1) = IntLoads.Mz;
end


if ~isempty(beam_model.Res.Bar.CForces)
    IntLoads = storeinternalloads(BARF);
    BarForces_Global.Fx(:,1) = IntLoads.Fx;BarForces_Global.Fy(:,1) = IntLoads.Fy;BarForces_Global.Fz(:,1) = IntLoads.Fz;
    BarForces_Global.Mx(:,1) = IntLoads.Mx;BarForces_Global.My(:,1) = IntLoads.My;BarForces_Global.Mz(:,1) = IntLoads.Mz;
else
    IntLoads = storeinternalloads(BEAMF);
    BarForces_Global.Fx(:,1) = IntLoads.Fx;BarForces_Global.Fy(:,1) = IntLoads.Fy;BarForces_Global.Fz(:,1) = IntLoads.Fz;
    BarForces_Global.Mx(:,1) = IntLoads.Mx;BarForces_Global.My(:,1) = IntLoads.My;BarForces_Global.Mz(:,1) = IntLoads.Mz;
end

ex_node = zeros(ngrid,3);
ey_node = zeros(ngrid,3);
ez_node = zeros(ngrid,3);

for k = 1:ngrid
    ex_node(k,1)=NR(1,1,k);
    ey_node(k,1)=NR(1,2,k);
    ez_node(k,1)=NR(1,3,k);
    ex_node(k,2)=NR(2,1,k);
    ey_node(k,2)=NR(2,2,k);
    ez_node(k,2)=NR(2,3,k);
    ex_node(k,3)=NR(3,1,k);
    ey_node(k,3)=NR(3,2,k);
    ez_node(k,3)=NR(3,3,k);
end

Displacements = zeros(ngrid,6);
Displacements(:,1:3) = beam_model.Res.NDispl(:,1:3);
Displacements(:,4) = -atan(ez_node(:,2)./ez_node(:,3));
Displacements(:,5) = atan(ez_node(:,1)./ez_node(:,3));
Displacements(:,6) = atan(ex_node(:,2)./ex_node(:,1));

Static.Displacements{TrimIdx} = Displacements;
Static.Positions{TrimIdx}     = NODEPOS;
Static.NonDimAeroForces{TrimIdx} = [beam_model.Res.Aero.results.spancd,...
    beam_model.Res.Aero.results.spancs,beam_model.Res.Aero.results.spancl];
Static.p_mid_r{TrimIdx}     = beam_model.Res.Aero.results.p_mid_r;
Static.BarForces{TrimIdx}   = BarForces;
Static.BarStrains{TrimIdx}  = BarStrains;
Static.BarForces_Global{TrimIdx} = BarForces_Global;
Static.State{TrimIdx}       = dummy_aero.state;
Static.Gamma{TrimIdx}       = results.gamma;
Static.Fext{TrimIdx}        = Fgrid;
Static.Lattice{TrimIdx}     = beam_model.Aero.lattice;
Static.spanarea{TrimIdx} = beam_model.Res.Aero.results.spanarea;

Static.spandrag{TrimIdx} = beam_model.Res.Aero.results.spandrag;
Static.spanside{TrimIdx} = beam_model.Res.Aero.results.spanside;
Static.spanlift{TrimIdx} = beam_model.Res.Aero.results.spanlift;
%beam_model.Results.flexible.NR      = NR;
%beam_model.Results.flexible.NODEPOS = NODEPOS;
%beam_model.Results.flexible.Res     = beam_model.Res;
%beam_model.Results.flexible.CL_alpha   = beam_model.Aero.CL_alpha;
%beam_model.Results.flexible.CL_control = beam_model.Aero.CL_control;

%fprintf('\n\t\t -%s: Tip Displacement = %8.4f m\n','Wing',Displacements(beam_model.WingNodes(end),3));


end
