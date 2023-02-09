%%OPTIMISETWIST: The following function solves the twist required for
%%a desired aero load distribution

function [Static,beam_model] = optimiseTwist(Static,beam_model,Aircraft_param,Sizing_param,load_idx,linear)

fp = 1;

AeroCount = 0;

EPS = 1e-3;

mesh.DW = [];
mesh.DW_sym = [];

% Initialise some counts
countaero    = 0;
countnpanels = 0;
countpanels  = 0;

q          = 0.5*beam_model.Aero.state.rho*beam_model.Aero.state.AS^2;
VLMspan    = zeros(sum(beam_model.Aero.Wing.geo.ny),sum(beam_model.Aero.geo.nx + beam_model.Aero.geo.fnx));
panel_area = tarea(beam_model.Aero.lattice_vlm.XYZ);
sp_mid     = (beam_model.Aero.lattice_vlm.XYZ(:,1,:)+beam_model.Aero.lattice_vlm.XYZ(:,2,:))/2;

% Generate a recovery matrix for the circulation and force
for ii = 1:size(beam_model.Aero.Wing.geo.ny,1)
    for m = 1:beam_model.Aero.Wing.geo.ny(ii,1)
        countaero = countaero + 1;
        countnpanels = countnpanels + 1; % groups main surface and elevator panels?
        numpanels = beam_model.Aero.Wing.geo.nx(ii,1) + beam_model.Aero.Wing.geo.fnx(ii,1);
        VLMspan(countaero,countpanels + 1:countpanels + numpanels) = 1;
        countpanels = countpanels + numpanels;
    end
end

span_area = VLMspan*panel_area';


% Perform an iteration loop where the goal is to drive down the required
% change in the twist to achieve our goal
while 1
    
    % Get the trim deflections so that we can get the gradient about the
    % correct deflection point
    [Static,beam_model]     = performTrim(Static,beam_model,load_idx,linear,fp);
    
    % Get the circulation and lift for the null condition
    [gamma,lift,mesh]  = getIQ(beam_model, Static.Positions{load_idx}, beam_model.Aero.geo.TWIST, linear, mesh);
    
    spangamma0 = VLMspan*gamma;
    spanlift0  = VLMspan*lift;
    
    nspan = sum(beam_model.Aero.Wing.geo.ny);
    
    % Initialise the gradient matrices
    dGamma_dTwist = zeros(nspan,nspan);
    dCl_dTwist    = zeros(nspan,nspan);
    
    %% TWIST PERTURBATION
    for sc = 1:nspan
        TWIST        = beam_model.Aero.geo.TWIST;
        TWIST(sc)    = beam_model.Aero.geo.TWIST(sc) + EPS;
        [gamma,lift,mesh] = getIQ(beam_model, Static.Positions{load_idx}, TWIST, linear, mesh);
        
        spangamma = VLMspan*gamma;
        spanlift  = VLMspan*lift;
        
        dGamma_dTwist(:,sc) = (spangamma-spangamma0)/EPS;
        dCl_dTwist(:,sc)    = (spanlift-spanlift0)./(q*span_area)/EPS;
    end
    
    % Recover the spanwise distance of each panel
    p_mid_y = VLMspan*sp_mid(:,1,2)./sum(VLMspan,2);
    
    % Calculate the inverse of the gradient matrix to save time later on
    invGamma_TWIST = inv(dGamma_dTwist);
    
    % Calculate the required twist to achieve an elliptical distribution
    Req_angle = beam_model.Aero.geo.TW(1,:,1);
    
    dTwist = zeros(sum(beam_model.Aero.geo.ny),1);
    
    % Solve the root circulation for the required root twist
    %     if isfield(Aircraft_param.Wing.Planform,'Cl')
    %
    %         Target_Cl = interp1(beam_model.Aero.ref.b_ref*Aircraft_param.Wing.Planform.Cl(:,1),Aircraft_param.Wing.Planform.Cl(:,2),p_mid_y,'linear','extrap');
    %
    %         Cl_0 = spanlift0./(q*span_area);
    %
    %         dTwist(1:nspan) = dCl_dTwist\((Target_Cl-Cl_0));
    %         dTwist(1:nspan) = dTwist(1:nspan) - dTwist(1);
    %
    %     else
    
    % This loophole has been added to take into account the fact that
    % target Cl has been selected without adding a distribution to target.
    if strcmp(Sizing_param.AeroDistribution,'TargetCl')
        if ~isfield(Aircraft_param.Wing.Planform,'Cl')
            Sizing_param.AeroDistribution = 'Elliptical';
        end
    end
    
    if strcmp(Sizing_param.AeroDistribution,'TargetCl')
        
        Target_Cl = interp1(beam_model.Aero.ref.b_ref*Aircraft_param.Wing.Planform.Cl(:,1),Aircraft_param.Wing.Planform.Cl(:,2),p_mid_y,'linear','extrap');
        
        Cl_0 = spanlift0./(q*span_area);
        
        dTwist(1:nspan) = dCl_dTwist\((Target_Cl-Cl_0));
        dTwist(1:nspan) = dTwist(1:nspan) - dTwist(1);
        
    elseif strcmp(Sizing_param.AeroDistribution,'Elliptical')
        Gamma_ROOT = (Req_angle - beam_model.Aero.geo.TWIST(1) + invGamma_TWIST(1,:) * spangamma0)/...
            (invGamma_TWIST(1,:) * sqrt(1-(p_mid_y/beam_model.Aero.ref.b_ref).^2));
        
        dTwist(1:nspan) = dGamma_dTwist\(Gamma_ROOT*(sqrt(1-(p_mid_y/beam_model.Aero.ref.b_ref).^2)) - spangamma0);
        
    elseif strcmp(Sizing_param.AeroDistribution,'Triangular')
        Gamma_ROOT = (Req_angle - beam_model.Aero.geo.TWIST(1) + invGamma_TWIST(1,:) * spangamma0)/...
            (invGamma_TWIST(1,:) * (1-(p_mid_y/beam_model.Aero.ref.b_ref)));
        
        dTwist(1:nspan) = dGamma_dTwist\(Gamma_ROOT*(1-(p_mid_y/beam_model.Aero.ref.b_ref)) - spangamma0);
        
    elseif strcmp(Sizing_param.AeroDistribution,'Bell-Shaped')
        
        span_theta = acos(p_mid_y/beam_model.Aero.ref.b_ref);
        
        Gdist = sin(span_theta).^3;
        
        Gamma_ROOT = (Req_angle - beam_model.Aero.geo.TWIST(1) + invGamma_TWIST(1,:) * spangamma0)/...
            (invGamma_TWIST(1,:) * Gdist);
        
        dTwist(1:nspan) = dGamma_dTwist\(Gamma_ROOT*Gdist - spangamma0);
    end
    %     end
    
    beam_model.Aero.geo.TWIST = beam_model.Aero.geo.TWIST + dTwist;
    
    beam_model.Aero.geo.TW(:,:,1) = beam_model.Aero.geo.TWIST;
    beam_model.Aero.geo.TW(:,:,2) = beam_model.Aero.geo.TWIST;
    
    Res_twist = norm(dTwist);
    %fprintf(fp,['\n... delta twist = ' num2str(Res_twist) ' \n']);
    
    if Res_twist < 10^(-4) || AeroCount > 5
        break
    end
    
    AeroCount = AeroCount +1;
end

[Static,beam_model] = performTrim(Static,beam_model,load_idx,Sizing_param.Linear,fp);

end

function [Static,beam_model] = performTrim(Static,beam_model,load_idx,Linear,fp)

if Linear == 1
    [Static,beam_model] = solve_linear_trim_customisedV2(Static,beam_model,load_idx);
else
    [Static,beam_model] = solve_nonlinear_trim_customised(Static,beam_model,load_idx);
end

end

function [gamma,lift,mesh] = getIQ(beam_model, NODEPOS, TWIST, Linear, mesh)

%% NULL CONDITION
dummy_aero = beam_model.Aero;
dummy_aero.state.alpha = beam_model.Res.state.alpha;
dummy_aero.state.betha = 0;
dummy_aero.state.P     = 0;
dummy_aero.state.Q     = 0;
dummy_aero.state.R     = 0;

beam_model.Aero.geo.TW(:,:,1) = TWIST;
beam_model.Aero.geo.TW(:,:,2) = TWIST;

% Update the lattice to take account for the twist and camber as well as
% the new deflection from the trim solution
if beam_model.Info.spline_type == 1
    lattice_defo = update_vlm_mesh1(beam_model.Node, beam_model.Res.SOL, beam_model.Aero,beam_model);
elseif beam_model.Info.spline_type == 2
    lattice_defo = update_vlm_mesh(beam_model.Node, NODEPOS, ...
        beam_model.Res.NRd, beam_model.Aero,beam_model);
else
    lattice_defo = update_vlm_mesh(beam_model.Node, NODEPOS, ...
        beam_model.Res.NRd, beam_model.Aero,beam_model);
end

% Modify the mesh depending on whether it is linear/nonlinear
if Linear == 1
    dummy_aero.lattice_vlm.N = lattice_defo.N;
    dummy_aero.lattice_vlm.N(:,2) = 0.*dummy_aero.lattice_vlm.N(:,2);
    dummy_aero.lattice_vlm.N(:,1) = dummy_aero.lattice_vlm.N(:,1)./sqrt(sum(dummy_aero.lattice.N.^2,2));
    dummy_aero.lattice_vlm.N(:,2) = dummy_aero.lattice_vlm.N(:,2)./sqrt(sum(dummy_aero.lattice.N.^2,2));
    dummy_aero.lattice_vlm.N(:,3) = dummy_aero.lattice_vlm.N(:,3)./sqrt(sum(dummy_aero.lattice.N.^2,2));
else
    dummy_aero.lattice_vlm = lattice_defo;
end

% Determine the global rotation matrix due to new angles of attack - THIS
% MAKES IT A NONLINEAR PROBLEM!!!
GlobalRotationMatrix = [cos(dummy_aero.state.betha)*cos(dummy_aero.state.alpha),        -sin(dummy_aero.state.betha),          cos(dummy_aero.state.betha)*sin(dummy_aero.state.alpha) ;...
    cos(dummy_aero.state.alpha)*sin(dummy_aero.state.betha),         cos(dummy_aero.state.betha),          sin(dummy_aero.state.betha)*sin(dummy_aero.state.alpha) ;...
    -sin(dummy_aero.state.alpha),                                              0,                           cos(dummy_aero.state.alpha)]';

% Determine the new flow vector
FlowVector = (beam_model.Aero.state.AS*GlobalRotationMatrix*[1,0,0]')';

% Run the vlm for the new parameters
results = VLM_HS(FlowVector,dummy_aero.state.rho,dummy_aero.lattice_vlm,beam_model.Aero.state.SIMXZ,[],[],mesh.DW,0,[],mesh.DW_sym);

% Post-process some of the results
results = Aero_process(results,dummy_aero.lattice_vlm,dummy_aero.state.rho,FlowVector,GlobalRotationMatrix,dummy_aero.ref.S_ref,beam_model.Aero.geo,beam_model.Res.WB.CG);

gamma  = results.gamma;
lift   = results.Fglobal(:,3);

mesh.DW     = results.DW;
mesh.DW_sym = results.DW_sym;

end