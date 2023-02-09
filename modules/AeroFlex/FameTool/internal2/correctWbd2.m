%% CORRECTWBD2   Apply mass corrections
% 
% Make a correction to the mass definition if needed. FAME defaults the CG
% position to the reference position if it is not defined in a weight
% function. You can define your own EXCEL spreadsheet with the following
% columns. The header will be ignored if you define column names.
%
% NAME, MASS, CGX, CGY, CGZ
%
% NAME : Name of part in FAME file
% MASS : Mass of part
% CGX  : CG position in x - axis
% CGY  : CG position in y - axis
% CGZ  : CG position in z - axis
%
% We need to tell the program that this file is present. We do this by
% setting the input flag and defining the path to the corrections file
%
%   Example:
%      obj = fame2mat;
%      obj.Opts.Struct.massCorrection = true;
%      obj.Inp.wbdCorrFiles = fullfile(FOLDERNAME,FILENAME);
%      obj = run(obj);
%
%
%   See also |getMassData| method within Fame2mat object.

%   Copyright 2016 University of Bristol
%   Private function.


function wbd = correctWbd2(wbd,wbdCorrFile)

if ~exist(wbdCorrFile,'file')
    wbdCorrFile = [];
end

fprintf('        Applying structural mass corrections ...\n');

[num,txt,~] = xlsread(wbdCorrFile);

if strcmp(txt{1},'NAME')
    txt(1,:) = [];
end

for i = 1:size(num,1)
    idx = find(strcmp(txt{i,1},{wbd.Name}));
    
    if ~isempty(idx)
        idx                   = idx(1);
        wbd(idx).Pts.x        = num(i,2);
        wbd(idx).Pts.y        = num(i,3);
        wbd(idx).Pts.z        = num(i,4);
        wbd(idx).Pts.Mass     = num(i,1);
        wbd(idx).Part.Mass    = num(i,1);
        wbd(idx).Part.CG      = num(i,2:4);
        wbd(idx).Part.Inertia = zeros(3);
    end
    
end

end

