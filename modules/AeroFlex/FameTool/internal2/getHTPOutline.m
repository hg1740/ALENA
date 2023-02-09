function Outline = getHTPOutline(fname)

if nargin < 1
    
   error('HTP Outline file not found');
    
end

keyword = { '%      x.outline       y.outline       z.outline'};

noutline = 0;

if ~exist(fname,'file')
    error('Unable to find file %s.\n', fname);
else
    fp       = fopen(fname, 'r');
    skipLine = false;
    
    while ~feof(fp)
        if ~skipLine
            tline = fgetl(fp);
        else
            skipLine = false;
        end
        CARD = tline;
        
        switch CARD
            % Outline card
            case keyword{1}
                fgetl(fp);
                fgetl(fp);
                tline    = fgetl(fp);
                noutline = noutline + 1;
                Outline(noutline,1:3) = strread(tline);
                tline    = fgetl(fp);
                
                % continuation detected
                while ~isequal(tline(1),'%')
                    noutline = noutline + 1;
                    Outline(noutline,1:3) = strread(tline);
                    tline    = fgetl(fp);
                end
                
                Outline  = Outline / 1000;
                skipLine = false;
        end
    end
    
    fclose(fp);
    
end

fclose all;



end