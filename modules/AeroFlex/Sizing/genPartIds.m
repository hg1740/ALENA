function beam_model = genPartIds(beam_model,LiftingParts)

PartIDs = fieldnames(beam_model.PartIDs);
beam_model.PartId = [];

for i = 1:numel(LiftingParts)
    
    % Identify if it is a lifting part
    idx = find(ismember(PartIDs,LiftingParts{i}));
    
    % If it finds a lifting part
    if ~isempty(idx)
        
        if ~strcmpi(LiftingParts{i},'vtp')
            
            % Find the right/left hand side nodes
            [~,nodeidx] = intersect(beam_model.Node.ID,beam_model.PartIDs.(LiftingParts{i}).BeamNodes);
            
            NodeCoord = beam_model.Node.Coord(nodeidx,2);
            
            Node_r = (NodeCoord>=0);
            Node_l = (NodeCoord<=0);
            
            nId = numel(beam_model.PartId);
            beam_model.PartId(nId + 1).Id = 1;
            beam_model.PartId(nId + 1).Part = [LiftingParts{i},'_R'];
            beam_model.PartId(nId + 1).Type = 'GRID';
            beam_model.PartId(nId + 1).data = beam_model.PartIDs.(LiftingParts{i}).BeamNodes(Node_r);
            beam_model.PartId(nId + 1).index = nodeidx(Node_r);
            
            nId = numel(beam_model.PartId);
            beam_model.PartId(nId + 1).Id = 1;
            beam_model.PartId(nId + 1).Part = [LiftingParts{i},'_L'];
            beam_model.PartId(nId + 1).Type = 'GRID';
            beam_model.PartId(nId + 1).data = beam_model.PartIDs.(LiftingParts{i}).BeamNodes(Node_l);
            beam_model.PartId(nId + 1).index = nodeidx(Node_l);
            
            % Find the right/left hand side beam elements
            [~,beamidx] = intersect(beam_model.Beam.ID,beam_model.PartIDs.(LiftingParts{i}).CBeam);
            
            BeamCoord = beam_model.Node.Coord(beam_model.Beam.Conn(beamidx,3),2);
            
            Beam_r = (BeamCoord > 0);
            Beam_l = (BeamCoord < 0);
            
            nId = numel(beam_model.PartId);
            beam_model.PartId(nId + 1).Id = 1;
            beam_model.PartId(nId + 1).Part = [LiftingParts{i},'_R'];
            beam_model.PartId(nId + 1).Type = 'CBEAM';
            beam_model.PartId(nId + 1).data = beam_model.PartIDs.(LiftingParts{i}).CBeam(Beam_r);
            beam_model.PartId(nId + 1).index = beamidx(Beam_r);
            
            nId = numel(beam_model.PartId);
            beam_model.PartId(nId + 1).Id = 1;
            beam_model.PartId(nId + 1).Part = [LiftingParts{i},'_L'];
            beam_model.PartId(nId + 1).Type = 'CBEAM';
            beam_model.PartId(nId + 1).data = beam_model.PartIDs.(LiftingParts{i}).CBeam(Beam_l);
            beam_model.PartId(nId + 1).index = beamidx(Beam_l);
            
            % Find the right/left hand side aero panels
            [~,aeroidx] = intersect(beam_model.Aero.ID,beam_model.PartIDs.(LiftingParts{i}).CAEROIDs);
            
            AeroCoord = beam_model.Aero.geo.starty(aeroidx) + beam_model.Aero.geo.b(aeroidx);
            
            Aero_r = (AeroCoord > 0);
            Aero_l = (AeroCoord < 0);
            
            nId = numel(beam_model.PartId);
            beam_model.PartId(nId + 1).Id = 1;
            beam_model.PartId(nId + 1).Part = [LiftingParts{i},'_R'];
            beam_model.PartId(nId + 1).Type = 'CAERO';
            beam_model.PartId(nId + 1).data = beam_model.PartIDs.(LiftingParts{i}).CAEROIDs(Aero_r);
            beam_model.PartId(nId + 1).index = aeroidx(Aero_r);
            
            nId = numel(beam_model.PartId);
            beam_model.PartId(nId + 1).Id = 1;
            beam_model.PartId(nId + 1).Part = [LiftingParts{i},'_L'];
            beam_model.PartId(nId + 1).Type = 'CAERO';
            beam_model.PartId(nId + 1).data = beam_model.PartIDs.(LiftingParts{i}).CAEROIDs(Aero_l);
            beam_model.PartId(nId + 1).index = aeroidx(Aero_l);
        else
                       
            nId = numel(beam_model.PartId);
            beam_model.PartId(nId + 1).Id = 1;
            beam_model.PartId(nId + 1).Part = [LiftingParts{i}];
            beam_model.PartId(nId + 1).Type = 'GRID';
            beam_model.PartId(nId + 1).data = beam_model.PartIDs.(LiftingParts{i}).BeamNodes;
            
            nId = numel(beam_model.PartId);
            beam_model.PartId(nId + 1).Id = 1;
            beam_model.PartId(nId + 1).Part = [LiftingParts{i}];
            beam_model.PartId(nId + 1).Type = 'CBEAM';
            beam_model.PartId(nId + 1).data = beam_model.PartIDs.(LiftingParts{i}).CBeam;
            
            nId = numel(beam_model.PartId);
            beam_model.PartId(nId + 1).Id = 1;
            beam_model.PartId(nId + 1).Part = [LiftingParts{i}];
            beam_model.PartId(nId + 1).Type = 'CAERO';
            beam_model.PartId(nId + 1).data = beam_model.PartIDs.(LiftingParts{i}).CAEROIDs;
            
        end
    end
end


end