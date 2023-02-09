%GETFAMEBOXPROPERTIES: A function that extracts the box properties of the
%   wing.

function [results,InputFiles] = getFameBoxGeometry(FolderName)

PathName = ['export',filesep,'fmwpp01',filesep,'flexible'];

FileName = 'wng_box_properties.fmwpp01';

InputFiles.boxPropertiesFiles{1} = fullfile(FolderName,PathName,FileName);

FilePath = [FolderName,filesep, PathName,filesep, FileName];

neta = 0;

if ~exist(FilePath)
    error('\nUnable to find file %s.', FileName);
else
    
    fp = fopen(FilePath, 'r');
    
    skip_line = false;
    
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
    
end

results = BoxProperties;

end