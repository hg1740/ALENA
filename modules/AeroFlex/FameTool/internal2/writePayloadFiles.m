function writePayloadFiles(homeFolder,obj)

% Findout how much fuel is in each file
for i = 1:numel(obj.Inp.NeofuelMassFiles)
    
    % Read the file into matlab
    blk = textread(obj.Inp.NeofuelMassFiles{i},'%s','delimiter','\n');
    
    % Carry out the regexp function on each line and check whether it is not empty
    idx    = find(~cellfun(@isempty,regexp(blk,'Total Fuel Mass:')));
    tab = blk(idx+1);
    FuelMassFile(i) = str2double(tab{1,1});
end

% Recover the mean aerodynamic chord
WingCAERO = findobj(obj.Mdl.Caero,'part','Wing');
[mac,mac_le,mac_ac] = recoverMAC(WingCAERO);

CG_ADJUSTMENT = 22.2 - (0.15*mac+mac_le); % 15% MAC
CG_ADJUSTMENT = 0;
% Run through the trim cases and generate an input file for each of them.

for trimI = 1:numel(obj.Fame.LoadCases)
    
    % Fuel Conm2 File
    FuelMass     = obj.Fame.LoadCases(trimI).fuel;
    
    if FuelMass > 0
        [~,FuelFileIdx]    = min(abs(FuelMass - FuelMassFile));
        [~,fueltype,~] = fileparts(obj.Inp.NeofuelMassFiles{FuelFileIdx});
        
    end
    
    fp = fopen(fullfile(homeFolder,['\Neo_model\Payload_Neo_LC' num2str(trimI) '.dat']), 'w');
    
    fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
    fprintf(fp,'$       CONM CARDS \n');
    fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
    
    PayloadIdx = find(~cellfun(@isempty,regexp({obj.Mdl.Conm2.type},'Payload')));
    
    % Get the aircraft conm2 indices for this loadcase
    AircraftConm2 = findobj(obj.Mdl.Conm2,'Type','Structural','-or','Type',...
        'Secondary','-or','Type','Engine','-or','Type',fueltype);
    
    FuseNode = findobj(obj.Mdl.Grid,'part','Fuselage');
    
    % Get the cg of the aircraft
    
    NodeId = [obj.Mdl.Grid.id];
    MCg = 0;
    M   = 0;
    for i = 1:numel(AircraftConm2)
       
        % Get the nodal coordinate
        ConmIdx = find(NodeId == AircraftConm2(i).grid);
        Conmcg = obj.Mdl.Grid(ConmIdx).coord;
        
        % Add the offset to get the correct cg
        Conmcg = Conmcg + AircraftConm2(i).offset;
        
        % Multiply the Cg by the mass and add
        MCg = MCg + Conmcg*AircraftConm2(i).m;
        M = M + AircraftConm2(i).m;
    end

    % Calculate the cg of the aircraft
    Cg = MCg/M;
    
    for i = PayloadIdx
        
        % This is a BODGE, because I don't necessarily have the correct
        % fuel files.
        
        %payloadmass = obj.Fame.LoadCases(trimI).AC_mass - obj.Fame.LoadCases(trimI).fuel - obj.Fame.Geometry.Weights.oew;
        if ~obj.Opts.Geom.addTail
            payloadmass = obj.Fame.LoadCases(trimI).Wing_lift/(9.81*obj.Fame.LoadCases(trimI).Loadfactor) - FuelMassFile(FuelFileIdx) - obj.Fame.Geometry.Weights.oew;
        else
            payloadmass = obj.Fame.LoadCases(trimI).AC_mass - FuelMassFile(FuelFileIdx) - obj.Fame.Geometry.Weights.oew;
        end
        
        Cg_payload = obj.Fame.LoadCases(trimI).AC_cg;
        
        % Calculate mac
        Aircraft_cg_x = Cg_payload*mac + mac_le + CG_ADJUSTMENT;
        Payload_cg    = [(Aircraft_cg_x*(M + payloadmass) - Cg(1)*M)/payloadmass,0,0];
        PayOffset     = Payload_cg - FuseNode.coord;
        
        % Update the inertial properties of the payload mass
        [~,Ixx,Iyy,Izz] = GuessFuselageMassMatrix(payloadmass,16,1.5,0.75);
        Inertia = diag([Ixx,Iyy,Izz]);
        
        % Update the payload mass object
        obj.Mdl.Conm2(i).m      = payloadmass;
        obj.Mdl.Conm2(i).offset = PayOffset;
        obj.Mdl.Conm2(i).i      = Inertia;
        
        % Edit the payload to account for the load case
        str = genNeocassDeck(obj.Mdl.Conm2(i));
        
        for j = 1:numel(str)
            fprintf(fp,'%s\n',str{j});
        end
    end
    
    fclose(fp);
end

end
