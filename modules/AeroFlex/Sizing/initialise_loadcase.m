function [beam_model,Aircraft_param] = initialise_loadcase(beam_model,LoadCases,load_idx,AddFuel,Optim,TotalFuel,Aircraft_param)

PlotFuel = 0;


beam_model.Param.G        = LoadCases(load_idx).URDD3;
beam_model.SolParam.trimI = load_idx;
beam_model                = setflightcond(beam_model,load_idx);
beam_model.Aero.lattice   = beam_model.Aero.lattice_vlm;
beam_model.Aero.AIC       = [];

% --------------------------
% Set the fuel configuration
% --------------------------
if AddFuel == 1
    
    fuelmass    = zeros(1,Optim.Wing.Nsec);
    
    FuelConf = LoadCases(load_idx).Fuel*(TotalFuel); % Half wing
    
    for j = 1:Optim.Wing.Nsec
        fuelmass(j) = (Optim.Wing.V_int(j)/sum(Optim.Wing.V_int))*FuelConf;
        Fuel_Idx = find(beam_model.Conm2.ID == beam_model.PartIDs.Fuel.Conm2IDs(j));
        beam_model.Conm2.M(1:3,1:3,Fuel_Idx) = (fuelmass(j))*eye(3);
    end
    
    if PlotFuel == 1
        figure(50);plot(Optim.Wing.y_ends(2:end),fuelmass./Optim.Wing.y_lbox','ro');
        title(['Fuel mass distribution for ' num2str(100*LoadCases(load_idx).Fuel) '% fuel']);
        xlabel('y [m]');ylabel('Fuel Mass per meter [Kg/m]')
    end
    
    % --------------------------
    % Temporarily remove payload properties
    % --------------------------    
    Payload_Idx = find(beam_model.Conm2.ID == beam_model.PartIDs.Payload.Conm2IDs);
    
    % Reassign the payload properties to the correct conm2 entry
    beam_model.Conm2.M(1:6,1:6,Payload_Idx)  = zeros(6);
    
    [beam_model.ConM,beam_model.Info] = ...
        sortCONMArray(beam_model.Conm1,beam_model.Conm2,beam_model.Node,beam_model.Coord,beam_model.Param,beam_model.Info);
    
    beam_model.WB = WBnCG(beam_model.Node,beam_model.Param,beam_model.ConM,beam_model.WB,beam_model.Bar,beam_model.Beam,beam_model.Info);
    
    % --------------------------
    % Allocate payload properties
    % --------------------------
    
    if isempty(LoadCases(load_idx).Payload)
        payloadmass = LoadCases(load_idx).AC_Mass/2 - beam_model.WB.MCG(1);
    else
        payloadmass = LoadCases(load_idx).Payload/2;
    end
    % Reassign the payload properties to the correct conm2 entry
    beam_model.Conm2.M(1:3,1:3,Payload_Idx)  = (payloadmass)*eye(3);
    
    % Calculate the
    Aircraft_cg_x   = LoadCases(load_idx).PayloadCG*beam_model.Aero.ref.MAC + beam_model.Aero.ref.MAC_LE_x;
    Payload_cg = [(Aircraft_cg_x*(beam_model.WB.MCG(1) + payloadmass) - beam_model.WB.CG(1)*beam_model.WB.MCG(1))/payloadmass,0,0]; 

    % Recalculate the offset of the lumped mass
    Parent          = Aircraft_param.Payload.Parent;
    BluffNodeIdx    = Aircraft_param.(Parent).Planform.NodeIdx + Aircraft_param.Payload.Planform.ParentIdx;
    NodeOffset      = Payload_cg - beam_model.Node.Coord((beam_model.Node.ID == BluffNodeIdx),:);
    
    radius = max(Aircraft_param.Payload.Planform.radius);
    
    % Calculate the inertial properties of the bluff body
    [~,Ixx,Iyy,Izz] = GuessFuselageMassMatrix(payloadmass,Aircraft_param.Payload.Planform.Length,radius,0.75);
    
    beam_model.Conm2.M(4:6,4:6,Payload_Idx)  = diag([Ixx,Iyy,Izz]);
    beam_model.Conm2.Offset(Payload_Idx,1:3) = NodeOffset;
    
    % --------------------------
    % Recalculate model mass properties
    % --------------------------
    [beam_model.ConM,beam_model.Info] = ...
        sortCONMArray(beam_model.Conm1,beam_model.Conm2,beam_model.Node,beam_model.Coord,beam_model.Param,beam_model.Info);
    
    beam_model.WB = WBnCG(beam_model.Node,beam_model.Param,beam_model.ConM,beam_model.WB,beam_model.Bar,beam_model.Beam,beam_model.Info);
    
    beam_model.Res.WB.CG = beam_model.WB.CG;
end

end