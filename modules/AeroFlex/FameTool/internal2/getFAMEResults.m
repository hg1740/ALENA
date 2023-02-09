function results = getFAMEResults(FolderName)

PathName='/export/fmwpp01/flexible';


%% FETCH AERODYNAMICS
FileName='/wng_ads_lc001.fmwpp01';
fid = 1;
PARAM.INCLUDE{1}=[FolderName PathName FileName];

READ_INCLUDE = true;
NFILE = 0;
neta=0;

if ~exist([FolderName PathName FileName])
    
    error('\nUnable to find file %s.', FileName);
else
    
    fp = fopen(PARAM.INCLUDE{1}, 'r');
    
    skip_line = false;
    fprintf(fid,'\nReading %s file...', FileName);
    while ~feof(fp)
        
        if ~skip_line
            tline = fgetl(fp);
        else
            skip_line = false;
        end
        
        stringdata=textscan(tline,'%f');
        if ~isempty(stringdata{1,1})
            neta=neta+1;
            Aerodynamic(neta,:)=stringdata{1,1}';
        end
        
    end % end of while
    fclose(fp);
    
end
%% FETCH DEFORMATIONS
FileName='/wng_deform_lc001.fmwpp01';
fid = 1;
PARAM.INCLUDE{1}=[FolderName PathName FileName];

READ_INCLUDE = true;
NFILE = 0;
neta=0;

if ~exist([FolderName PathName FileName])
    
    error('\nUnable to find file %s.', FileName);
else
    
    fp = fopen(PARAM.INCLUDE{1}, 'r');
    
    skip_line = false;
    fprintf(fid,'\nReading %s file...', FileName);
    while ~feof(fp)
        
        if ~skip_line
            tline = fgetl(fp);
        else
            skip_line = false;
        end
        
        stringdata=textscan(tline,'%f');
        if ~isempty(stringdata{1,1})
            neta=neta+1;
            Deformation(neta,:)=stringdata{1,1}';
        end
        
    end % end of while
    fclose(fp);
    
end

%% FETCH BOX PROPERTIES
FileName='/wng_box_properties.fmwpp01';
fid = 1;
PARAM.INCLUDE{1}=[FolderName PathName FileName];

READ_INCLUDE = true;
NFILE = 0;
neta=0;

if ~exist([FolderName PathName FileName])
    
    error('\nUnable to find file %s.', FileName);
else
    
    fp = fopen(PARAM.INCLUDE{1}, 'r');
    
    skip_line = false;
    fprintf(fid,'\nReading %s file...', FileName);
    while ~feof(fp)
        
        if ~skip_line
            tline = fgetl(fp);
        else
            skip_line = false;
        end
        
        stringdata=textscan(tline,'%f');
        if ~isempty(stringdata{1,1})
            neta=neta+1;
            BoxProperties(neta,:)=stringdata{1,1}';
        end
        
    end % end of while
    fclose(fp);
    %
    
end

%% FETCH THICKNESS AND CHORD PROPERTIES
FileName='/wng_addgeodata.fmwpp01';
fid = 1;
PARAM.INCLUDE{1}=[FolderName PathName FileName];

READ_INCLUDE = true;
NFILE = 0;
neta=0;

if ~exist([FolderName PathName FileName])
    
    error('\nUnable to find file %s.', FileName);
else
    
    fp = fopen(PARAM.INCLUDE{1}, 'r');
    
    skip_line = false;
    fprintf(fid,'\nReading %s file...', FileName);
    while ~feof(fp)
        
        if ~skip_line
            tline = fgetl(fp);
        else
            skip_line = false;
        end
        
        stringdata=textscan(tline,'%f');
        if ~isempty(stringdata{1,1})
            neta=neta+1;
            AirfoilProperties(neta,:)=stringdata{1,1}';
        end
        
    end 
    fclose(fp);
end
results.LoadCases=LoadCases;
results.Aerodynamic=Aerodynamic;
results.Deformation=Deformation;
results.BoxProperties=BoxProperties;
end