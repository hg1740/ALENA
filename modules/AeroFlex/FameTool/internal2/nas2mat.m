%% NAS2MAT
%
% This function creates NEOCASS input decks from the FAME model object.

function Nastran = nas2mat( varargin )

if nargin>0
    fileName = varargin{1};
else
    [fileName,pathName] = uigetfile({'*.dat';'*.bdf'},'NASTRAN Input File');
    
    if fileName == 0
        return
    else
        fileName = fullfile(pathName,fileName);
    end
end

% Initialise
Nastran = [];

% Read file into cell array and remove comments etc
[dwc,d] = readfile(fileName);

Nastran.File.contents = d;

% Parse mat cards
Nastran.Mat = readmat(dwc);

% Parse grid data
Nastran.Grid = readgrid(dwc);

% Parse CBAR data
Nastran.Cbar = readcbar(dwc);

% Parse PBAR data
Nastran.Pbar = readpbar(dwc);

% Parse FORCE data
Nastran.Force = readforce(dwc);

% Parse MOMENT data
Nastran.Moment = readmoment(dwc);

% Parse SPC data
%Nastran.Spc = readspc(dwc);

end

%--------------------------------------------------------------------------
function [dwc,d] = readfile(fileName)

fid      = fopen(fileName,'r');
d        = textscan(fid,'%s','delimiter','\n');
fclose(fid);
d        = d{1};
idx      = strmatch('$',d);
dwc      = d;
dwc(idx) = [];
idx      = find(cellfun(@isempty,dwc));
dwc(idx) = [];

end

%--------------------------------------------------------------------------
function MaterObj = readmat(d)

idx = find(~cellfun(@isempty,regexp(d,'MAT1')));

for i = 1:numel(idx)
    blk = d(idx(i):idx(i) + 1);
    blk = regexprep(blk,'MAT1\*','    ');
    blk = regexprep(blk,'\*',' ');
    blk = cell2mat(blk');
    blk = regexprep(blk,'E\+(\w*)\+','E\+$1 \+');
    blk = regexprep(blk,'E\+(\w*)\-','E\+$1 \-');
    blk = regexprep(blk,'E\-(\w*)\+','E\-$1 \+');
    blk = regexprep(blk,'E\-(\w*)\-','E\-$1 \-');
    blk = strtrim(blk);
    blk = regexprep(blk,'                ','  0.0  ');
    blk = str2num(blk);
    
    % Create mat object. G set by NU in object and stress limit
    % values automatically set to 0 on object creation.
    MaterObj(i)      = Mat;
    MaterObj(i).id   = blk(1);
    MaterObj(i).e    = blk(2);
    MaterObj(i).nu   = blk(4);
    MaterObj(i).rho  = blk(5);
    MaterObj(i).a    = blk(6);
    MaterObj(i).tRef = blk(7);
    
end

end

%--------------------------------------------------------------------------
function GridObj = readgrid(d)

idx = find(~cellfun(@isempty,regexp(d,'GRID\*')));

for i = 1:numel(idx)
    blk = d(idx(i):idx(i) + 1);
    blk = regexprep(blk,'GRID\*','    ');
    blk = regexprep(blk,'\*',' ');
    blk = cell2mat(blk');
    blk = regexprep(blk,'E\+(\w*)\+','E\+$1 \+');
    blk = regexprep(blk,'E\+(\w*)\-','E\+$1 \-');
    blk = regexprep(blk,'E\-(\w*)\+','E\-$1 \+');
    blk = regexprep(blk,'E\-(\w*)\-','E\-$1 \-');
    blk = str2num(blk);
    
    % create grid object
    GridObj(i)       = Grid;
    GridObj(i).id    = blk(1);
    GridObj(i).cs    = blk(2);
    GridObj(i).coord = [blk(3),blk(4),blk(5)];
    GridObj(i).cd    = blk(6);
    GridObj(i).ps    = blk(7);
    GridObj(i).seid  = blk(8);
    
end

end

%--------------------------------------------------------------------------
function PbarObj = readpbar(d)

idx = find(~cellfun(@isempty,regexp(d,'PBAR\*')));

for i = 1:numel(idx)
    blk = d(idx(i):idx(i) + 4);
    blk = regexprep(blk,'PBAR\*','    ');
    blk = regexprep(blk,'\*',' ');
    blk = cell2mat(blk');
    blk = regexprep(blk,'E\+(\w*)\+','E\+$1 \+');
    blk = regexprep(blk,'E\+(\w*)\-','E\+$1 \-');
    blk = regexprep(blk,'E\-(\w*)\+','E\-$1 \+');
    blk = regexprep(blk,'E\-(\w*)\-','E\-$1 \-');
    blk = str2num(blk);
    
    % create pbar object
    PbarObj(i)          = Pbar;
    PbarObj(i).id       = blk(1);
    PbarObj(i).si       = [];
    PbarObj(i).type     = [];
    PbarObj(i).mat      = blk(2);
    PbarObj(i).a        = blk(3);
    PbarObj(i).a        = 1e9;
    PbarObj(i).i        = [blk(4),blk(5),blk(18)];
    PbarObj(i).j        = blk(6);
    PbarObj(i).kShear   = [];
    PbarObj(i).rhoNs    = [];
    PbarObj(i).startPnt = [];
end


end

%--------------------------------------------------------------------------
function CbarObj = readcbar(d)

idx = find(~cellfun(@isempty,regexp(d,'CBAR\*')));

for i = 1:numel(idx)
    blk = d(idx(i):idx(i) + 1);
    blk = regexprep(blk,'CBAR\*','    ');
    blk = regexprep(blk,'\*',' ');
    blk = regexprep(blk,'GGG','');
    blk = cell2mat(blk');
    blk = regexprep(blk,'E\+(\w*)\+','E\+$1 \+');
    blk = regexprep(blk,'E\+(\w*)\-','E\+$1 \-');
    blk = regexprep(blk,'E\-(\w*)\+','E\-$1 \+');
    blk = regexprep(blk,'E\-(\w*)\-','E\-$1 \-');
    blk = str2num(blk);
    
    % create bar object
    CbarObj(i)        = Cbar;
    CbarObj(i).id     = blk(1);
    CbarObj(i).pid    = blk(2);
    CbarObj(i).conn   = [blk(3),blk(4)];
    CbarObj(i).barg0  = [];
    CbarObj(i).orient = [blk(5),blk(6),blk(7)];
    CbarObj(i).offset = 'GGG';
end

end

%--------------------------------------------------------------------------
function ForceObj = readforce(d)

idx = find(~cellfun(@isempty,regexp(d,'FORCE\*')));

for i = 1:numel(idx)
    blk = d(idx(i):idx(i) + 1);
    blk = regexprep(blk,'FORCE\*','    ');
    blk = regexprep(blk,'\*',' ');
    blk = cell2mat(blk');
    blk = regexprep(blk,'E\+(\w*)\+','E\+$1 \+');
    blk = regexprep(blk,'E\+(\w*)\-','E\+$1 \-');
    blk = regexprep(blk,'E\-(\w*)\+','E\-$1 \+');
    blk = regexprep(blk,'E\-(\w*)\-','E\-$1 \-');
    blk = str2num(blk);
    
    % create force object
    ForceObj(i)          = Force;
    ForceObj(i).id       = blk(1);
    ForceObj(i).gid      = blk(2);
    ForceObj(i).cid      = blk(3);
    ForceObj(i).forceMag = blk(4);
    ForceObj(i).x        = blk(5);
    ForceObj(i).y        = blk(6);
    ForceObj(i).z        = blk(7);
end


end

%--------------------------------------------------------------------------
function MomentObj = readmoment(d)

idx = find(~cellfun(@isempty,regexp(d,'MOMENT\*')));

for i = 1:numel(idx)
    blk = d(idx(i):idx(i) + 1);
    blk = regexprep(blk,'MOMENT\*','    ');
    blk = regexprep(blk,'\*',' ');
    blk = cell2mat(blk');
    blk = regexprep(blk,'E\+(\w*)\+','E\+$1 \+');
    blk = regexprep(blk,'E\+(\w*)\-','E\+$1 \-');
    blk = regexprep(blk,'E\-(\w*)\+','E\-$1 \+');
    blk = regexprep(blk,'E\-(\w*)\-','E\-$1 \-');
    blk = str2num(blk);
    
    % create moment object
    MomentObj(i)           = Moment;
    MomentObj(i).id        = blk(1);
    MomentObj(i).gid       = blk(2);
    MomentObj(i).cid       = blk(3);
    MomentObj(i).momentMag = blk(4);
    MomentObj(i).x         = blk(5);
    MomentObj(i).y         = blk(6);
    MomentObj(i).z         = blk(7);
end

end

%--------------------------------------------------------------------------
function SpcObj = readspc(d)

idx = find(~cellfun(@isempty,regexp(d,'SPC  ')));

for i = 1:numel(idx)
    blk = d(idx(i));
    blk = regexprep(blk,'SPC     ','');
    
    % create spc object
    SpcObj(i)       = Spc;
    SpcObj(i).id    = str2num(blk{1}(1:8));
    SpcObj(i).grids = str2num(blk{1}(9:16));
    SpcObj(i).dof   = str2num(blk{1}(17:24));
    SpcObj(i).d     = str2num(blk{1}(25:32));
end

end

