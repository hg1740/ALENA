function [results,InputFiles] = getFameDeformations(FolderName)

% The pathname is essentially hardcoded 
PathName = ['export',filesep,'fmwpp01',filesep,'flexible'];

list = dir(fullfile(FolderName,PathName,'wng_deform_lc*.fmwpp01'));

for i = 1:length(list)

    FileName = list(i).name;
    
    FilePath = [FolderName,filesep, PathName,filesep, FileName];
    InputFiles.wingDeformationsFiles{i} = fullfile(FilePath);
      
    Deformation = [];
    
    if ~exist(InputFiles.wingDeformationsFiles{i},'file')
        error('Unable to find file %s.\n', InputFiles.wingDeformationsFiles{i});
    else
        
        s      = textread(InputFiles.wingDeformationsFiles{i} ,'%s','delimiter','\n');
        idx    = find(strncmp(s,'%',1));
        s(idx) = [];
        s      = str2num(char(s));
        
        Deformation.eta      = s(:,1);
        Deformation.Uxd2     = s(:,2);
        Deformation.Uxd1     = s(:,3);
        Deformation.Ux       = s(:,4)/1000;
        Deformation.Uzd2     = s(:,5);
        Deformation.Uzd1     = s(:,6);
        Deformation.Uz       = s(:,7)/1000;
        Deformation.Thetad1  = s(:,8);
        Deformation.Theta    = s(:,9);
 
    end
    results(i) = Deformation;
end

end