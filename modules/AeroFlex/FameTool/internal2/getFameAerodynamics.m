%GETFAMEAERODYNAMICS: A function that extracts the aerodynamic loading and
%   twist from the FAME output files.
%
% Date of Creation: 12/2015
%
% Authors:  Dario Calderon - University of Bristol
%           Etienne Coetzee - Airbus (FPO)
%
%
function [results,InputFiles] = getFameAerodynamics(FolderName)

PathName = ['export',filesep,'fmwpp01',filesep,'flexible'];

list = dir(fullfile(FolderName,PathName,'wng_ads_lc*.fmwpp01'));

for i = 1:length(list)
    
    FileName = list(i).name;
    
    FilePath = [FolderName,filesep, PathName,filesep, FileName];
    
    Aerodynamic  = [];
    neta         = 0;
    
    InputFiles.aeroCasesFiles{1} = fullfile(FilePath);
    
    if ~exist(FilePath)
        error('\nUnable to find file %s.', FileName);
    else
             
        s      = textread(InputFiles.aeroCasesFiles{1} ,'%s','delimiter','\n');
        idx    = find(strncmp(s,'%',1));
        s(idx) = [];
        s      = str2num(char(s));
        
        Aerodynamic.eta         = s(:,1);
        Aerodynamic.Cl          = s(:,2);
        Aerodynamic.Cl_l        = s(:,3);
        Aerodynamic.dCl         = s(:,4);
        Aerodynamic.Cm          = s(:,5);
        Aerodynamic.Cm_l        = s(:,6);
        Aerodynamic.dCm         = s(:,7);
        Aerodynamic.Cdi         = s(:,8);
        Aerodynamic.Cdf         = s(:,9);
        Aerodynamic.Cd          = s(:,10);
        Aerodynamic.jig_tw      = s(:,11);
        Aerodynamic.dTw_b       = s(:,12);
        Aerodynamic.dTw_t       = s(:,13);
        Aerodynamic.dTw_bt      = s(:,14);
        Aerodynamic.dTw_man     = s(:,15);
        Aerodynamic.act_tw      = s(:,16);
    end
    
    results(i) = Aerodynamic;
end


end