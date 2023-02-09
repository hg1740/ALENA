%% GETFAMEBOXAREASFCN Obtain area of strigers and skins along the span from FAME
%
% Source FAME file : '1_structure\1_10_wing\l3_fame-w\flexible\box1_thickness+area.txt'

%   Copyright 2016 University of Bristol
%   Private function.
function [Results,Inputfiles] = getFameBoxAreasFcn(folderName)

datPathName = ['1_structure',filesep,'1_10_wing',filesep,'l3_fame-w',filesep,'flexible',filesep,'box1_after_sizing'];

boxAreaFile = fullfile(folderName,datPathName,'box1_thickness+area.txt');

Inputfiles.boxAreaFiles = boxAreaFile;

s = textread(boxAreaFile,'%s','delimiter','\n');

str      = 'bm  by bcs    ksi    eta   eta* n_ustr';
idx      = find(strncmp(s,str,numel(str)));
blk      = s(idx + 3:end);

str      = '============';
idx      = find(strncmp(blk,str,numel(str)));
blk(idx) = [];

s = str2num(char(blk));

Results.eta     = s(:,5);
Results.aUshTot = s(:,13);
Results.aLshTot = s(:,20);
Results.a       = Results.aUshTot + Results.aLshTot;

end

