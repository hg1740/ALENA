%GETFAMEBOXPROPERTIES: A function that extracts the box properties of the
%   wing.

function [results,InputFiles] = getFameBoxGeometry(FolderName)

PathName = ['export',filesep,'fmwpp01',filesep,'flexible'];

FileName = 'wng_skin_thickness.fmwpp01';

FilePath = [FolderName,filesep,PathName,filesep,FileName];

InputFiles.thicknessFiles{1} = fullfile(FilePath);

Thickness = [];

if ~exist(InputFiles.thicknessFiles{1},'file')
    
    warning('Unable to find file %s.\n', InputFiles.thicknessFiles{1});

else
    
    s      = textread(InputFiles.thicknessFiles{1} ,'%s','delimiter','\n');
    idx    = find(strncmp(s,'%',1));
    s(idx) = [];
    s      = str2num(char(s));
    
    Thickness.eta      = s(:,1);
    Thickness.Uskn     = s(:,4);
    Thickness.Ush      = s(:,8);
    Thickness.Lskn     = s(:,12);
    Thickness.Lsh      = s(:,16);
    Thickness.Fsp      = s(:,20);
    Thickness.Rsp      = s(:,24);
end

results = Thickness;

end