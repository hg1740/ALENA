%%ALLOCATEMASSES:
% 
% 
%   Author: Dario Calderon 

function [beam_model,SecondaryMass] = AllocateMasses(beam_model,Aircraft_param,Parameters,Optim_Variables)

Plot_fuel = 0;

AircraftParts = fieldnames(Aircraft_param);

for i = 1:length(AircraftParts)
    if isfield(Aircraft_param.(AircraftParts{i}),'LiftingSurface')
        if Aircraft_param.(AircraftParts{i}).LiftingSurface == 0
            
            if isempty(Aircraft_param.(AircraftParts{i}).Parent)
                BluffMass = Aircraft_param.(AircraftParts{i}).Properties.Mass / 2;
            else
                BluffMass = Aircraft_param.(AircraftParts{i}).Properties.Mass;
            end
            
            radius = max(Aircraft_param.(AircraftParts{i}).Planform.radius);
            
            % Calculate the inertial properties of the bluff body
            [~,Bluff_Ixx,Bluff_Iyy,Bluff_Izz] = GuessFuselageMassMatrix(BluffMass,Aircraft_param.(AircraftParts{i}).Planform.Length,radius,0.75);
            
            if ~isempty(Aircraft_param.(AircraftParts{i}).Parent)
                Parent          = Aircraft_param.(AircraftParts{i}).Parent;
                BluffNodeIdx    = Aircraft_param.(Parent).Planform.NodeIdx + Aircraft_param.(AircraftParts{i}).Planform.ParentIdx;
                NodeOffset      = Aircraft_param.(AircraftParts{i}).Properties.CG - beam_model.Node.Coord((beam_model.Node.ID == BluffNodeIdx),:);
            else
                BluffNodeIdx = Aircraft_param.(AircraftParts{i}).Planform.NodeIdx+1;
                NodeOffset = [0,0,0];
            end
            
            % Populate the conm2 arrays with the new lumped masses
            [beam_model.Conm2,beam_model.Info] = ...
                CONM2Creator(beam_model.Conm2,beam_model.Info,Aircraft_param.(AircraftParts{i}).Planform.CONM2Idx+1,BluffNodeIdx,...
                BluffMass,NodeOffset,Bluff_Ixx,Bluff_Iyy,Bluff_Izz);
            
            beam_model.PartIDs.(AircraftParts{i}).Conm2IDs = Aircraft_param.(AircraftParts{i}).Planform.CONM2Idx+1;
            
        end
    end
end

% ---------------------------
% Add Secondary Mass (CONM2)
% ---------------------------
if Parameters.AddSecondMass == 1
    [SecondaryMass,beam_model] = SecMass(beam_model,Aircraft_param,Optim_Variables);
else
    SecondaryMass.Total = 0;
end

% Post-process arrays
if Parameters.AddFuel == 1
    fuelmass = zeros(1,Optim_Variables.Wing.Nsec);
    beam_model.PartIDs.Fuel.Conm2IDs = [];
    
    FuelConf = 0.5*(0.5*Aircraft_param.Weights.FB); % 50% fuel case for twist optimisation
    for j = 1:Optim_Variables.Wing.Nsec
        fuelmass(j) = (Optim_Variables.Wing.V_int(j)/sum(Optim_Variables.Wing.V_int))*FuelConf;
        fuelnode1 = beam_model.PartIDs.Wing.BeamNodes(j+1);
        fuelnode2 = 610000 + j;
        [beam_model.Conm2,beam_model.Info] = ...
            CONM2Creator(beam_model.Conm2,beam_model.Info,fuelnode2,fuelnode1,...
            fuelmass(j),[0,0,0]);
        beam_model.PartIDs.Fuel.Conm2IDs = [beam_model.PartIDs.Fuel.Conm2IDs,fuelnode2];
    end
    if Plot_fuel == 1
        figure(50);plot(Optim_Variables.Wing.y_ends(2:end),fuelmass./Optim_Variables.Wing.y_lbox','ro');
        title('Fuel mass distribution for 50% fuel');
        xlabel('y [m]');ylabel('Fuel Mass per meter [Kg/m]')
    end
end


end