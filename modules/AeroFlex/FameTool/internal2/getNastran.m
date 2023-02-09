%% GETNASTRAN Extract NASTRAN data into a model object
%
% The NASTRAN file is parsed and then the information is placed into the
% relevant objects in the |Model| object. The numbering in FAME is
% arbitrary, hence the objects are renumbered starting from:
%
%   Structural Nodes Right Wing    : 110001  
%   Structural Nodes Left Wing     : 120001
%   Rbe0 Aero Points LE Right Wing : 211001
%   Rbe0 Aero Points TE Right Wing : 212001
%   Rbe0 Aero Points LE Left Wing  : 221001
%   Rbe0 Aero Points TE Left Wing  : 222001
%   Engine Right                   : 310001
%   Engine Left                    : 320001
%
%   See also |getFameStructure| method within Fame2mat object.

%   Copyright 2016 University of Bristol
%   Private function.
function obj = getNastran(obj, logfcn)

if nargin < 2
    logfcn = @(s) fprintf('%s\n', s);
end
            
folderName = obj.Dirs.fameFolder;

datPathName = ['1_structure',filesep,'1_10_wing',filesep,'l3_fame-w',filesep,'flexible',filesep,'deformations'];
nastranFile = fullfile(folderName,datPathName,'model_nastran__ldcase0001.dat');

obj.Inp.nastranFiles = nastranFile;

% Read nastran file
logfcn(' Parsing NASTRAN file...');
nas = nas2mat(nastranFile);
obj.Mdl = nas;

% Renumber grid and cbar elements here to make it easier to reflect the
% model later. 
startIdOffsetNumber = 110000;
nGrid = numel(obj.Mdl.Grid);

for i = 1:nGrid
    obj.Mdl.Grid(i).id = obj.Mdl.Grid(i).id + startIdOffsetNumber;
    obj.Mdl.Grid(i).part = 'StbdWing';
end

nSpc = numel(obj.Mdl.Spc);
for i = 1:nSpc
    obj.Mdl.Spc(i).id    = obj.Mdl.Spc(i).id + startIdOffsetNumber;
    obj.Mdl.Spc(i).grids = obj.Mdl.Spc(i).grids + startIdOffsetNumber;
    obj.Mdl.Spc(i).part = 'StbdWing';
end

nBar = numel(obj.Mdl.Cbar);
for i = 1:nBar
    obj.Mdl.Cbar(i).id      = obj.Mdl.Cbar(i).id + startIdOffsetNumber;
    obj.Mdl.Cbar(i).pid     = obj.Mdl.Pbar(i).id + startIdOffsetNumber;
    obj.Mdl.Cbar(i).conn(1) = obj.Mdl.Cbar(i).conn(1) + startIdOffsetNumber;
    obj.Mdl.Cbar(i).conn(2) = obj.Mdl.Cbar(i).conn(2) + startIdOffsetNumber;
    obj.Mdl.Cbar(i).part = 'StbdWing';
end

nPbar = numel(obj.Mdl.Pbar);
for i = 1:nPbar
    obj.Mdl.Pbar(i).id = obj.Mdl.Pbar(i).id + startIdOffsetNumber;
    obj.Mdl.Pbar(i).part = 'StbdWing';
end

nForce = numel(obj.Mdl.Force);
for i = 1:nForce
    obj.Mdl.Force(i).gid = obj.Mdl.Force(i).gid + startIdOffsetNumber;
end

nMoment = numel(obj.Mdl.Moment);
for i = 1:nMoment
    obj.Mdl.Moment(i).gid = obj.Mdl.Moment(i).gid + startIdOffsetNumber;
end

logfcn(' NASTRAN file succesfully parsed.');










