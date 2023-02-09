function data = splitfuns(varargin)

data = readfile('fameConverter.m');
idx = find(strncmp(data,'function ',9));

idx(1) = 1;
idx = cat(1,idx,numel(data) + 1);

fprintf('\n\n');

for i = 2:numel(idx) - 1
    
    str = data(idx(i):idx(i + 1) - 1);
    
    eqidx = strfind(str{1,:},' = ');
    bridx = strfind(str{1,:},'(');
    fcnName = strtrim(str{1}(eqidx + 1:bridx - 1));
    
    fid = fopen([fcnName,'.m'],'w');
    fprintf('Creating function: %s...\n',[fcnName,'.m']);
    
    for j = 1:numel(str)
        fprintf(fid,'%s\n',str{j});
    end
    
    fclose(fid);
end

end

function data = readfile(fileName)

fid = fopen(fileName);
data = cell(0,1);
i = 1;

while feof(fid) == 0
    data{i,1} = fgetl(fid);
    i = i + 1;
end

fclose(fid);

end
