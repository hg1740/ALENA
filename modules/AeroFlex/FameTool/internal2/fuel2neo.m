%% FUEL2NEO Generate fuel Conm2 cards
%
% Use fuel data from FAME to create Conm2 cards in NASTRAN format. Each
% Conm2 object has a type field in the object. All the fuel related objects
% are found and then written into separate files, based on their
% classification. The output is written into the output directory.
function fuel2neo(writeFolder,obj)

lType = {obj.Mdl.Conm2.type}';
uType = unique(lType);

idxType = find(~cellfun(@isempty,regexpi(uType,'fuel')));

for j = 1:numel(idxType)
    
    idx = find(~cellfun(@isempty,regexp(lType,uType{idxType(j)})));
    
    fileName = fullfile(writeFolder,['Neo_model' filesep uType{idxType(j)},'.dat']);
    
    
    INPUTFILES.NeofuelMassFiles{j} = fileName;
    
    fp = fopen(fileName, 'w');
    
    fuelmass = sum([obj.Mdl.Conm2(idx).m]);
    
    fprintf(fp,'\n');
    fprintf(fp,'Total Fuel Mass: \n');
    fprintf(fp,'%8.4f\n',fuelmass);
    fprintf(fp,'\n');
    
    fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
    fprintf(fp,'$  CONM CARDS \n');
    fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
    
    for i = 1:numel(idx)
        
        str = genNeocassDeck(obj.Mdl.Conm2(idx(i)));
        
        for k = 1:numel(str)
            fprintf(fp,'%s\n',str{k});
        end
    end
    
    fclose all;
end

obj.Inp = setFields(obj.Inp,INPUTFILES);


end
