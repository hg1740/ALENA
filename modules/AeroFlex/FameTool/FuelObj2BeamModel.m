%% FUELOBJ2BEAMMODEL:
%   Summary:
%       This funtion interrogates the fuel distribution objects and adds
%       the lumped mass to the beam model structure.

function beam_model = FuelObj2BeamModel(FuelDistribution,beam_model)

% Get the fuel distributions for the current loadcase
nfuel = numel(FuelDistribution);

for i = 1:nfuel
    
    % Get a list of new IDs, Grids, and coordinate frames
    newid    = [FuelDistribution(i).AeroFlexConm2.id]';
    newnodes = [FuelDistribution(i).AeroFlexConm2.grid]';
    newcid   = [FuelDistribution(i).AeroFlexConm2.cid]';
        
    %Check for row/column mismatch
    %   - TODO : This should be removed eventually
    if isrow(beam_model.Conm2.ID)
        beam_model.Conm2.ID = beam_model.Conm2.ID';
    end
    if isrow(beam_model.Conm2.Node)
        beam_model.Conm2.Node = beam_model.Conm2.Node';
    end
    if isrow(beam_model.Conm2.CID)
        beam_model.Conm2.CID = beam_model.Conm2.CID';
    end
    
    % Aggregate the exisiting fields
    beam_model.Conm2.ID   = [beam_model.Conm2.ID;newid];
    beam_model.Conm2.Node = [beam_model.Conm2.Node;newnodes];
    beam_model.Conm2.CID  = [beam_model.Conm2.CID;newcid];    
    
    % Add the 6x6 mass matrices to the structure
    for j = 1:numel(newid)
        M = [FuelDistribution(i).AeroFlexConm2(j).m*eye(3),zeros(3);zeros(3),FuelDistribution(i).AeroFlexConm2(j).i];
        beam_model.Conm2.M(:,:,end+1) = M;
    end
    
    % Pass the fuel offsets into the Conm2 structure
    beam_model.Conm2.Offset = [beam_model.Conm2.Offset;vertcat(FuelDistribution(i).AeroFlexConm2.offset)];
    
end

% Update the conm2 count
beam_model.Info.ncom2 = length(beam_model.Conm2.ID);

end