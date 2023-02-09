% RIGID_TRIM  Rigid Trim solver
%
%   Summary: The following function performs a rigid trim for a given
%   aerodynamic model etc.
%
%   Aero solvers: 
%          VLM - HS
%          VLM - VR
%          Strip Theory
%
%   Inputs:
%      - beam_model
%      - Results    - FM
%                   - CS
%      - Aero
%
%   Outputs:
%      - beam_model
%
%   Authors: -Originally found in NeoCASS
%            -Heavily modified by Dario Calderon (01/12/2017) - University of
%           Bristol
%
function Results = rigid_trim(beam_model,Results,AERO)

% Reset the AIC matrix at the beginning as this may need to be recalculated
Results.Aero.AIC = [];

% Reset the induced downwash at the beginning
Results.Aero.IndDownwash = [];

% Determine number of free trim variables
nfm_fixed = sum(Results.FM.Fixed);
nfm_free  = length(Results.FM.Fixed) - nfm_fixed;

% Determine number of free control surfaces
if AERO.geo.nc
    ncs_fixed = sum(Results.CS.Fixed);
    ncs_free  = length(Results.CS.Fixed) - ncs_fixed;
else
    ncs_free = 0;
end

% Determine the total number of free DOFs for trim
NDOF = nfm_free + ncs_free;

% Intialise the trim solution (based on NDOF)
sol  = zeros(NDOF,1);

% Define delta for DOF
EPS  = 1e-5;

% Define the tolerance on the norm of the forces
TOLL = 1e-3;

% Define the maximum number of steps 
NSTEP = 30;

% Initialise the trim Jacobian
J = [];

% Begin the trim iteration loop
for N = 1:NSTEP

    [RES0,Results] = ...
        get_rigid_vlm_res(beam_model,Results,beam_model.Aero,sol,NDOF);

    NORM = norm(RES0);
    
    % Check convergence
    if (NORM < TOLL)
        break;
    end
    
    % determine numeric jacobian
    if mod(N,10) == 1 
        
        %J = zeros(NDOF, NDOF);
        
        for k = 1:NDOF
            dsol = sol;
            dsol(k) = sol(k) + EPS;
            
            [RES,Results] =...
                get_rigid_vlm_res(beam_model,Results,beam_model.Aero,dsol,NDOF);
            
            J(:,k) = (RES-RES0)./EPS;
        end
    end
    
    if size(J,1)< NDOF
        sol = sol -pinv(J)*RES0;
    else
        sol = sol - J\RES0;
    end
end

if (NORM > TOLL)
    error(['Trim was not possible ... Norm = ' num2str(NORM)]);
end

Results.Aero.J = J;

end

%% Get the trim residual in vertical force
function [RES,TrimResults] = get_rigid_vlm_res(beam_model,TrimResults,Aero,sol,NDOF)

TOTFM = 15;

%RES = zeros(NDOF,1);

count = 0;

nr = find(Aero.Trim.CS.MPC == 0);

for n = 1:TOTFM % flight mechanics parameters
    if ~TrimResults.FM.Fixed(1,n)
        count = count +1;
        TrimResults.FM.Value(1,n) = sol(count);
    end
end

% Get and set the rigid body trim variables
TrimResults.state.alpha = TrimResults.FM.Value(2);
TrimResults.state.betha = TrimResults.FM.Value(3);

[TrimResults.state.alpha, TrimResults.state.betha,TrimResults. state.P, TrimResults.state.Q, TrimResults.state.R] = ...
    get_state_trim_vars(TrimResults.FM);

if Aero.geo.nc
    
    % Find the Master CS and substitute with the NR variables
    [~,J] = find(TrimResults.CS.Fixed == 0);
    TrimResults.CS.Value(J) = sol(count+1:end);
    
    % Get and set the control surface deflections
    TrimResults.CS.Value = apply_constr_eq(Aero.geo.nc, TrimResults.CS.Fixed, Aero.Trim.CS.MPC, Aero.Trim.CS.Coeff, TrimResults.CS.Value);
    
    % Update the lattice normals to account for the change in control
    % surface deflection
    
    if beam_model.MeshType == 3 % VLM - HS
        beam_model.Aero.lattice_defo = rotate_control_norm(Aero.ref, TrimResults.state, Aero.geo, ...
            Aero.lattice, TrimResults.CS.Value, ...
            Aero.lattice.Control.Hinge);
        
    elseif beam_model.MeshType == 4 % VLM - VR

        beam_model.Aero.lattice_defo = rotate_VR_norm(beam_model.Aero.ref,beam_model.Res.state,beam_model.Aero.geo,...
            beam_model.Aero.lattice,beam_model.Res.CS.Value,...
            beam_model.Aero.lattice.Control.Hinge,beam_model.Aero.lattice.Control.HingeVort);
    
    else % Don't bother with the strip theory
        beam_model.Aero.lattice_defo = beam_model.Aero.lattice;
    end
    
else
    beam_model.Aero.lattice_defo = Aero.lattice;
end


% TrimResults.state.alpha = pi*10/180;

GlobalRotationMatrix = [cos(TrimResults.state.betha)*cos(TrimResults.state.alpha),        -sin(TrimResults.state.betha),          cos(TrimResults.state.betha)*sin(TrimResults.state.alpha) ;...
                        cos(TrimResults.state.alpha)*sin(TrimResults.state.betha),         cos(TrimResults.state.betha),          sin(TrimResults.state.betha)*sin(TrimResults.state.alpha) ;...
                              -sin(TrimResults.state.alpha),                                              0,                           cos(TrimResults.state.alpha)]';

FlowVector = (Aero.state.AS*GlobalRotationMatrix*[1,0,0]')';


switch beam_model.MeshType
    
    case {1,2} % Strip Theory
        [results.la,results.Wind,ez_aero,ey_aero] = fLocal_AoA2(beam_model.Res.NRd,beam_model,beam_model.Aero.lattice_defo,TrimResults.state,beam_model.Gust,[],[],0);
        
        results = stripAero(TrimResults.state, beam_model.Aero.geo, beam_model.Aero.lattice_defo,results,...
            beam_model.SolParam.AeroLinear,beam_model.Aero.CL_alpha,beam_model.Aero.CL_control,TrimResults.CS,nr,...
            beam_model.Aero.geo.panelarea,ez_aero,ey_aero,beam_model.Res.AeroRd);
        
    case 3 % Horseshoe vortex
        results = VLM_HS(FlowVector,beam_model.Aero.state.rho,beam_model.Aero.lattice_defo,Aero.state.SIMXZ,[],TrimResults.Aero.AIC,[],0,TrimResults.Aero.IndDownwash);

        
    case 4 % Vortex Ring Elements
        [results,beam_model] = VLM_VR(state, beam_model.Aero.geo, beam_model.Aero.lattice_defo,beam_model.Gust,beam_model.Aero.Interp.Ic,[],beam_model.SolParam.Gust,beam_model);
end

% simplify this step
results = Aero_process(results,beam_model.Aero.lattice_defo,beam_model.Aero.state.rho,FlowVector,GlobalRotationMatrix,Aero.ref.S_ref,beam_model.Aero.geo,beam_model.Res.WB.CG);

% Fa = gf_nl_aero_nodal(beam_model.Info, beam_model.Node.DOF, ...
%                 beam_model.Node, beam_model.Aero, results,beam_model.SolParam.AeroLinear,beam_model.Res.NRd);
% 
% for n = 1:beam_model.Info.ngrid
%     dof = beam_model.Node.DOF(n, 1:6);
%     index = find(dof);
%     if ~isempty(index)
%         Fgrid(n, index) = Fa(dof(index));
%     end
% end

%results = Aero_coeff(results,beam_model.Aero.lattice_vlm,beam_model.Aero.lattice_defo,beam_model.Res.state,beam_model.Aero.ref,beam_model.Aero.geo,beam_model.SolParam.Linear,beam_model.SolParam.AeroLinear);

% Work out the difference in forces and moments
G = beam_model.Param.G;

M  = beam_model.WB.MCG(1);

if NDOF == 1
    % Lift Equilibrium
    RES = results.L - M * G;
else
    % Lift Equilibrium
    RES(1,1) = results.L - M * G;
    
    % Pitching Moment Equilibrium (neglects angular velocities)
    RES(2,1) = results.TotalMoment(2);
end

TrimResults.Aero = results;

end