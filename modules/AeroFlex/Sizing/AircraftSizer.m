%AIRCRAFTSIZER: A routine that outputs a fully sized aircraft
%
%
%   Author: Dario Calderon

function [Static,Optim_Variables,beam_model] = AircraftSizer(GUI_Input,hp, logfcn)

%Log function not supplied ?
if nargin < 3
    
    %Go with good old disp
    logfcn = @disp;
    
end

%% Extract the required structures from the main input structure
LoadCases       = GUI_Input.LoadCases;
beam_model      = GUI_Input.beam_model;
Aircraft_param  = GUI_Input.Aircraft_param;
Sizing_param    = GUI_Input.Parameters;
Optim_Variables = GUI_Input.Optim_Variables;

if nargin < 2
    hp = [];
end

%Sizing_param.Linear = 1;

% This is currently needed so that the trim analysis is treated properly
beam_model.SolParam.Trim = 1;

%% Predefine the static and dynamic outputs
Static  = [];
Dynamic = [];

%% Clear variable to free up space
clear GUI_Input

%% Extract the aircraft parts that can be optimised
LiftingParts  = fieldnames(Optim_Variables);
AircraftParts = fieldnames(Aircraft_param);

%% Setup the load cases [beam_model.Aero.Trim(:)]
n_loads  = length(LoadCases);
for idx = 1:n_loads
    beam_model = SetupLoadCases(beam_model,idx,LoadCases(idx));
end

%% Calculate Section Properties
for i = 1:length(LiftingParts)   
    Optim_Variables.(LiftingParts{i}) = BeamReduction(Optim_Variables.(LiftingParts{i}));
end

% Set the flight conditions to the last load case (cruise)
beam_model = setflightcond(beam_model,n_loads);

%% Parasitic Drag Estimation
Cd0_local = [];
Cd0       = 0;
for i = 1:length(AircraftParts)
    if isfield(Aircraft_param.(AircraftParts{i}),'LiftingSurface')
        if Aircraft_param.(AircraftParts{i}).LiftingSurface == 1
            
            % Calculate some necessary geometries of the wing profile
            Aircraft_param.(AircraftParts{i}).AeroSurface = AeroSurfaceProperties(beam_model.Aero.(AircraftParts{i}).geo,...
                Aircraft_param.(AircraftParts{i}).Planform);
            
            % Calculate the parasitic drag
            Aircraft_param.(AircraftParts{i}).AeroSurface = ParasiteDragWing(beam_model.Aero.(AircraftParts{i}).geo,...
                beam_model.Aero.state,beam_model.Aero.ref,Aircraft_param.(AircraftParts{i}).AeroSurface);
            
            Cd0_local = [Cd0_local;Aircraft_param.(AircraftParts{i}).AeroSurface.cd0]; %#ok<AGROW>
        else
            Aircraft_param.(AircraftParts{i}).AeroSurface = ParasiteDragFuse( Aircraft_param.(AircraftParts{i}).Planform,beam_model.Aero.state,beam_model.Aero.ref);
        end
        
        Cd0  = Cd0 + Aircraft_param.(AircraftParts{i}).AeroSurface.Cd0;
    end
end

%%
logfcn('Generating variables for nonlinear beam solver');
% -----------------------------------------
% Create Grids, Bars, RBE0s,SET1 and SPLINE1 for the lifting sections
% -----------------------------------------
for i = 1:length(LiftingParts)
    beam_model = GenerateLiftStructure(beam_model,Optim_Variables.(LiftingParts{i}),...
        Aircraft_param.(LiftingParts{i}),(LiftingParts{i}),Sizing_param);
end

% -----------------------------------------
% Add SPCs for static analysis
% -----------------------------------------
% Define the boundary conditions for the aircraft. Including rigid
% Connections
beam_model =  GenerateBCs(beam_model,Aircraft_param);

%% Allocate masses to Conm2 entries
[beam_model,SecondaryMass] = AllocateMasses(beam_model,Aircraft_param,Sizing_param,Optim_Variables);

beam_model = sortArrays(beam_model);

%% Extract Ids? This should not be necessary
StartNode = find(beam_model.Node.ID == beam_model.PartIDs.Wing.BeamNodes(1));
EndNode   = find(beam_model.Node.ID == beam_model.PartIDs.Wing.BeamNodes(end));
if Sizing_param.BeamType == 0
    StartBar = find(beam_model.Bar.Conn(:,1) == StartNode);
    EndBar   = find(beam_model.Bar.Conn(:,3) == EndNode);
else
    StartBar = find(beam_model.Beam.Conn(:,1) == StartNode);
    EndBar   = find(beam_model.Beam.Conn(:,3) == EndNode);
end

% Store the bars and nodes for plotting purposes
beam_model.WingNodes   = (StartNode:EndNode);
beam_model.WingMeanBar = (StartBar:EndBar);
beam_model.WingBar     = (StartBar:2*EndBar);

[beam_model.Aero,beam_model.Info] = ...
    sortINTERPArray(beam_model.Aero,beam_model.Param,beam_model.Info);

%% I need to divide

beam_model.WB = WBnCG(beam_model.Node,beam_model.Param,beam_model.ConM,beam_model.WB,beam_model.Bar,beam_model.Beam,beam_model.Info);

beam_model = PrintMasses(beam_model);

beam_model.Aero.Interp.Temp = 0;

% Choose which aerodynamic solver to use
beam_model.MeshType = 3;

%Sets the spline matrix??
beam_model = aeroelastic_interface(beam_model);

%%  Initialise the optimisation variables
for i = 1:length(LiftingParts)
    
    % Initialise the SOL Variable
    Optim_Variables.(LiftingParts{i}) = InitialiseOptimVariables(Optim_Variables.(LiftingParts{i}));
    
    % Initialise the Mass
    Optim_Variables.(LiftingParts{i}).Mass(1) = ObjFunc(Optim_Variables.(LiftingParts{i}).SOL, Optim_Variables.(LiftingParts{i}));
    
    Optim_Variables.(LiftingParts{i}).Iteration = 0;
end

% Find the Wing CAERO cards
wingidx = [];
for i = 1:length(beam_model.PartIDs.Wing.CAEROIDs)
    wingidx = [wingidx,find(beam_model.Aero.ID == beam_model.PartIDs.Wing.CAEROIDs(i))];
end

% Calculate the aerodynamic centres and store somewhere.
[mac, mac_LE, mac_AC] = CalculateMAC(beam_model.Aero,wingidx);

beam_model.Aero.ref.MAC       = mac;
beam_model.Aero.ref.MAC_LE_x  = mac_LE;
beam_model.Aero.ref.MAC_ac    = mac_AC;

% beam_model.Aero.geo.TWIST = 0.5*(squeeze(beam_model.Aero.geo.TW(:,:,1))+squeeze(beam_model.Aero.geo.TW(:,:,2)));
beam_model.Res.WB.CG = beam_model.WB.CG;

beam_model.SolParam.Aeroelastic = 1;
beam_model = set_state_matrices(beam_model);

% -------------------------------
% Store sizing history
% -------------------------------
for i = 1:length(LiftingParts)
    SizingHistory.(LiftingParts{i}) = store_sizing_history(Optim_Variables.(LiftingParts{i}));
end

Static = [];

iter_idx = 0;

while 1
    
    iter_idx = iter_idx + 1;
    logfcn(' ');
    logfcn(['----------- Iteration # ' num2str(iter_idx) ' -----------']);
    
    % --------------------
    % Run Jig Twist Optimisation
    % --------------------------
    %Sizing_param.Aero_optimisation = 1;
    if Sizing_param.Aero_optimisation == 1
        
        % Perform a jig twist optimisation here:
        logfcn(' ');
        logfcn('... Performing a jig twist correction');
        
        % Assign a mass configuration:
        [beam_model,Aircraft_param] = initialise_loadcase(beam_model,LoadCases,n_loads,...
            Sizing_param.AddFuel,Optim_Variables,0.5*Aircraft_param.Weights.FB,Aircraft_param);
        
        beam_model.Res = rigid_trim(beam_model,beam_model.Res,beam_model.Aero);
        
        [Static,beam_model] = optimiseTwist(Static,beam_model,Aircraft_param,Sizing_param,n_loads,Sizing_param.Linear);
        
    end
    
    
    % --------------------------
    % Run the Sizing Loads
    % --------------------------
    logfcn('Extracting Sizing Loads...');
    
    for iload = 1:n_loads-1
        
        [beam_model,Aircraft_param] = initialise_loadcase(beam_model,LoadCases,iload,...
            Sizing_param.AddFuel,Optim_Variables,0.5*Aircraft_param.Weights.FB,Aircraft_param);
        
        beam_model.Res = rigid_trim(beam_model,beam_model.Res,beam_model.Aero);
        
        % Perform a nonlinear trim analysis
        logfcn(sprintf('Performing flexible %-.2f g trim analysis', LoadCases(iload).URDD3/9.81));
        
        if Sizing_param.Linear == 1
            [Static,beam_model] = solve_linear_trim_customisedV2(Static,beam_model,iload);
        else
            [Static,beam_model] = solve_nonlinear_trim_customised(Static,beam_model,iload,hp,logfcn);
        end
        
        % Plot and store loads that are fed into the optimiser
        Optim_Variables = store_internal_loads(beam_model,Optim_Variables,LiftingParts,iload);
        
    end
    
    % --------------------------
    % Create an envelope of the internal loads
    % --------------------------
    Optim_Variables = store_envelope_internal_loads(Optim_Variables,LiftingParts);
    
    if Sizing_param.Struc_optimisation == 1
        
        % --------------------------
        % Run the optimisation on the different parts of the aircraft
        % --------------------------
        logfcn('Performing a Structural Sizing');
        
        for i = 1:length(LiftingParts)
            if Aircraft_param.(LiftingParts{i}).Optimisation == 1
                
                [Optim_Variables.(LiftingParts{i}),~,~] = run_optimisation(Optim_Variables.(LiftingParts{i}));
                
                logfcn(['Structural ' LiftingParts{i} ' mass (half model) = ' sprintf('%.2f',Optim_Variables.Wing.fval) ' Kg']);
                
                % Recover the section properites
                Optim_Variables.(LiftingParts{i}).Mass(iter_idx+1) = Optim_Variables.(LiftingParts{i}).fval;
                Optim_Variables.(LiftingParts{i}).X0 = Optim_Variables.(LiftingParts{i}).SOL;
                Optim_Variables.(LiftingParts{i}) = RecoverBoxProp(Optim_Variables.(LiftingParts{i}));
            end
        end
        
        
        % --------------------------
        % Update the beam model to account for the new stiffness and mass
        % properties
        % --------------------------
        for i = 1:length(LiftingParts)
            if Aircraft_param.(LiftingParts{i}).Optimisation == 1
                beam_model = UpdateBeamProperties(beam_model,Optim_Variables.(LiftingParts{i}),Aircraft_param.(LiftingParts{i}),LiftingParts{i},Sizing_param);
            end
        end
        
        % Update the element stiffness and mass matrices
        beam_model.Beam = set_cbeam_database(beam_model.Info.nbeam, beam_model.Beam, beam_model.PBeam, beam_model.Mat, beam_model.Node, []);
        
        for i = 1:length(LiftingParts)
            if Aircraft_param.(LiftingParts{i}).Optimisation == 1
                SizingHistory.(LiftingParts{i}) = UpdateSizingHistory(SizingHistory.(LiftingParts{i}),...
                    Optim_Variables.(LiftingParts{i}));
            end
        end
        
        if min(abs(diff(Optim_Variables.Wing.Mass))) < 0.1 || iter_idx > 20
            break;
        end
    else
        break;
    end
end

% Perform 1g maneouver for Breguet Range Calculation
logfcn('Extracting Design Cruise Data ...');

[beam_model,Aircraft_param] = initialise_loadcase(beam_model,LoadCases,n_loads,...
    Sizing_param.AddFuel,Optim_Variables,0.5*Aircraft_param.Weights.FB,Aircraft_param);

beam_model.Res = rigid_trim(beam_model,beam_model.Res,beam_model.Aero);


% Perform trim for Breguet Range Design Case
if Sizing_param.Linear == 1
    [Static,beam_model] = solve_linear_trim_customisedV2(Static,beam_model,n_loads);
else
    [Static,beam_model] = solve_nonlinear_trim_customised(Static,beam_model,n_loads,hp, logfcn);
end

Optim_Variables = store_internal_loads(beam_model,Optim_Variables,LiftingParts,n_loads);

% Calculate Parasitic drag as a funciton of Cl
beam_model.Results.flexible{1,n_loads}.Aero.results.Cdp = Cd0 + 0.38*Cd0.*(beam_model.Res.Aero.results.CN.^2);

logfcn(sprintf('Final Full Aircraft Mass at (0.5 fuel) = %-8.4g Kg',2*beam_model.WB.MCG(1)));

% ----------------------------
%
% ----------------------------

beam_model.WB = WBnCG(beam_model.Node,beam_model.Param,beam_model.ConM,beam_model.WB,beam_model.Bar,beam_model.Beam,beam_model.Info);
beam_model = PrintMasses(beam_model);

% Recover the twist for the different wing parts for reflecting purposes
if ~isempty(Aircraft_param.Wing)
    aero_index = [];
    for j = 1:length(beam_model.PartIDs.Wing.CAEROIDs)
        index = find(beam_model.Aero.ID == beam_model.PartIDs.Wing.CAEROIDs(j)) ;
        aero_index = [aero_index,index];
    end
    
    beam_model.Aero.Wing.geo.TW(:,1,1) = beam_model.Aero.geo.TWIST(aero_index);
    beam_model.Aero.Wing.geo.TW(:,1,2) = beam_model.Aero.geo.TWIST(aero_index);
end

if isfield(Aircraft_param,'StbdHTP')
    if ~isempty(Aircraft_param.StbdHTP)
        aero_index = [];
        for j = 1:length(beam_model.PartIDs.StbdHTP.CAEROIDs)
            index = find(beam_model.Aero.ID == beam_model.PartIDs.StbdHTP.CAEROIDs(j)) ;
            aero_index = [aero_index,index];
        end
        
        beam_model.Aero.StbdHTP.geo.TW(:,1,1) = beam_model.Aero.geo.TWIST(aero_index);
        beam_model.Aero.StbdHTP.geo.TW(:,1,2) = beam_model.Aero.geo.TWIST(aero_index);
    end
end

% --------------------------------- %
% Sort out the VTP Splines
% --------------------------------- %
% There should be flag that ignores the aero if it has a symmetry flag
if isfield(Aircraft_param,'VTP')
    if ~isempty(Aircraft_param.VTP)
        
        beam_model = generateVTPAero(beam_model,Optim_Variables.VTP,Aircraft_param,Sizing_param,n_loads);
        
    end
end

if Sizing_param.Ref_model == 1
    beam_model = reflectmodel(beam_model,Aircraft_param,Sizing_param);
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%% ReGenerate the VORTEX LATTICE
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
[beam_model.Aero.lattice_vlm, beam_model.Aero.ref] = ...
    vlm_setup(1, beam_model.Aero.geo, beam_model.Aero.state, beam_model.Aero.ref);

beam_model.Info.ncaero = numel(beam_model.Aero.geo.ny);

beam_model = genPartIds(beam_model,LiftingParts);
logfcn('Sizing Finished!');

end