%% GETFAMEBOXSTIFFFCN Obtain box stiffness along span from FAME

function [Stiffness,Inputfiles] = getFameStiffness(folderName)


pathName = ['export',filesep,'fmwpp01',filesep,'flexible'];
fileName = 'wng_box_properties.fmwpp01';

Inputfiles.boxPropertiesFiles{1} = fullfile(folderName,pathName,fileName);

if ~exist(Inputfiles.boxPropertiesFiles{1},'file')
    error('Unable to find file %s.\n', Inputfiles.boxPropertiesFiles{1});
else
    
    s      = textread(Inputfiles.boxPropertiesFiles{1} ,'%s','delimiter','\n');
    idx    = find(strncmp(s,'%',1));
    s(idx) = [];
    s      = str2num(char(s));
    
    Stiffness.eta  = s(:,1);
    Stiffness.EIxx = s(:,2)/1000000;
    Stiffness.EIxz = s(:,3)/1000000;
    Stiffness.EIzz = s(:,4)/1000000;
    Stiffness.EKxx = s(:,5)/1000000;
    Stiffness.EKzz = s(:,6)/1000000;
    Stiffness.GJ   = s(:,7)/1000000;
    Stiffness.EA   = s(:,8);
    Stiffness.Abox = s(:,9);
    Stiffness.Xec  = s(:,10)/1000;
    Stiffness.Zec  = s(:,11)/1000;
    Stiffness.Xsc  = s(:,12)/1000;
    Stiffness.Zsc  = s(:,13)/1000;
    
end

end

