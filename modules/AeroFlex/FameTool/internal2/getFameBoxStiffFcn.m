%% GETFAMEBOXSTIFFFCN Obtain box stiffness along span from FAME
%
% Source FAME file : 'export\fmwpp01\flexible\wng_box_properties.fmwpp01'

%   Copyright 2016 University of Bristol
%   Private function.
function [Results,Inputfiles] = getFameBoxStiffFcn(folderName)

datPathName = ['export',filesep,'fmwpp01',filesep,'flexible'];
boxStiffnessFile = fullfile(folderName,datPathName,'wng_box_properties.fmwpp01');

Inputfiles.boxStiffnessFiles = boxStiffnessFile;

s      = textread(boxStiffnessFile,'%s','delimiter','\n');
idx    = find(strncmp(s,'%',1));
s(idx) = [];
s      = str2num(char(s));

Results.eta  = s(:,1);
Results.EIxx = s(:,2);
Results.EIxz = s(:,3);
Results.EIzz = s(:,4);
Results.EKxx = s(:,5);
Results.EKzz = s(:,6);
Results.GJ   = s(:,7);
Results.EA   = s(:,8);
Results.Abox = s(:,9);
Results.Xec  = s(:,10);
Results.Zec  = s(:,11);
Results.Xsc  = s(:,12);
Results.Zsc  = s(:,13);

end

