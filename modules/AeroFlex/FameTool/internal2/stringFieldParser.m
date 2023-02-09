function str = stringFieldParser(line, index)

FIELD = 16;

if length(line) < 8 + FIELD * (index - 2)
    
    str = [];
    
else
    
    field = line(8 + (index - 2) * FIELD + 1:end);
    
    if ~isempty(field)
        
        if length(field) > FIELD
            str = strtok(char(field(1:FIELD)));
        else
            str = strtok(char(field));
        end
        
    else
        
        str = [];
        
    end
end

end
