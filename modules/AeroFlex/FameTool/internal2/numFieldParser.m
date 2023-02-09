function num = numFieldParser(line, index)

FIELD = 16;

if length(line) < 8 + FIELD * (index - 2)
    num = 0;
else
    minc  = min(length(line), 8 + (index - 1) * FIELD);
    field = strtok(line(8 + (index - 2) * FIELD + 1:minc));
    
    if ~isempty(field)
        num = str2num(field);
    else
        num = 0;
    end
    
end

end

