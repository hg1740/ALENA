%GETFAMEINTERNALLOADS: A function that extracts the internal loads from the
%   FAME output files.
%
% Date of Creation: 01/2018
%
% Authors:  Dario Calderon - University of Bristol
%
%
function [results,InputFiles] = getFameInternalLoads(FolderName)

PathName = ['export',filesep,'fmwpp01',filesep,'flexible'];

list = dir(fullfile(FolderName,PathName,'wng_cld_lc*.fmwpp01'));

% Ignore the first file
list(1) = [];

for i = 1:length(list)
    
    FileName = list(i).name;
    
    FilePath = [FolderName,filesep, PathName,filesep, FileName];
    
    InternalLoads  = [];
    
    InputFiles.internalLoadsFiles{1} = fullfile(FilePath);
    
    if ~exist(FilePath)
        error('\nUnable to find file %s.', FileName);
    else
             
        s      = textread(InputFiles.internalLoadsFiles{1} ,'%s','delimiter','\n');
        idx    = find(strncmp(s,'%',1));
        s(idx) = [];
        s      = str2num(char(s));
        
        InternalLoads.eta     = s(:,1);
        InternalLoads.Fx      = s(:,2);
        InternalLoads.Fy      = s(:,3);
        InternalLoads.Fz      = s(:,4);
        InternalLoads.Mx      = s(:,5);
        InternalLoads.My      = s(:,6);
        InternalLoads.Mz      = s(:,7);
        
    end
    
    results(i) = InternalLoads;
end


end