function writeMainInputfile(homeFolder,obj)

% Findout how much fuel is in each file
for i = 1:numel(obj.Inp.NeofuelMassFiles)
    
    % Read the file into matlab
    blk = textread(obj.Inp.NeofuelMassFiles{i},'%s','delimiter','\n');
    
    % Carry out the regexp function on each line and check whether it is not empty
    idx    = find(~cellfun(@isempty,regexp(blk,'Total Fuel Mass:')));
    tab = blk(idx+1);
    FuelMassFile(i) = str2double(tab{1,1});
end

for trimI = 1:numel(obj.Fame.LoadCases)
    fp = fopen(fullfile(homeFolder,['Neo_model\MainFile_LC' num2str(trimI) '.dat']), 'w');
    
    fprintf(fp,'$---------- MAIN INPUT FILE --------------$\n');
    % Trim input file
    fprintf(fp,'INCLUDE~ TrimInput_LC%i.dat\n',trimI);
    
    % Structural and aerodynamic model
    fprintf(fp,'INCLUDE~ FameModel_Neo.dat\n');
    
    % Payload File
    fprintf(fp,'INCLUDE~ Payload_Neo_LC%i.dat\n',trimI);
    
    % Fuel Conm2 File
    FuelMass     = obj.Fame.LoadCases(trimI).fuel;
    if FuelMass > 0 
        [~,FuelIdx] = min(abs(FuelMass - FuelMassFile));
        [~,n,e] = fileparts(obj.Inp.NeofuelMassFiles{FuelIdx});
        fprintf(fp,'INCLUDE~ %s\n',[n,e]);
    end
    
    fclose(fp);
end
end
