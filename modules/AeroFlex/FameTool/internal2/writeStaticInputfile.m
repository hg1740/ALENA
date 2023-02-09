function writeStaticInputfile(homeFolder,obj)

fp = fopen(fullfile(homeFolder,'StaticInputFile.inc'), 'w');
fprintf(fp,'SOL 144\n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$------- SPC PROPERTIES --------------------------------------------------\n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');

spcCases = unique([obj.Mdl.Spc.id]);

for i = 1:numel(spcCases)
    fprintf(fp,'SPC= %1.0f\n',spcCases(i));
end

for i = 1:numel(obj.Mdl.Spc)
    
    str = genNeocassDeck(obj.Mdl.Spc(i));
    
    for j = 1:numel(str)
        fprintf(fp,'%s\n',str{j});
    end
end

fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$------- DLOAD PROPERTIES ------------------------------------------------\n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'LOAD= 1\n'); % TODO - Removed space Dario -> Etienne

idx = find([obj.Mdl.Spc.dof] == 135);

for i = 1:numel(idx)
    fprintf(fp,'DLOAD   1       %-8d  3       0       0       0       0       0  \n',obj.Mdl.Spc(idx(i)).grids);
    fprintf(fp,'        0       0       0 \n');
end

fclose(fp);
end
