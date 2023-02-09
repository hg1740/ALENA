%% READWBD2SORTED   Read sorted FAME structural mass data file
%
%  READWBD2SORTED(CLASSNAME) reads mass data from FAME WBD2 file. Thousands
%  of small mass elements are created. The coordinates and mass values are
%  defined in the WBD2 file.
%
%  CLASSNAME is a Fame2mat object
%
%  This is a private function and should be accessed from the |Fame2mat|
%  object.
%
%   Example:
%      obj = fame2mat;
%      obj = getFameResults(obj); % calls this function
%
%

%   Copyright 2016 University of Bristol
%   Private function.
function d = readWbd2Sorted( varargin )

if ~isempty(varargin)
    fileName = varargin{1};
    pathName = [];
else
    [fileName,pathName] = uigetfile('*sorted.205','FAME Mass File');
    
    if fileName == 0
        return
    end
end

% Read the file into matlab
blk = textread(fullfile(pathName,fileName),'%s','delimiter','\n');

% Carry out the regexp function on each line and check whether it is not empty
idx    = find(~cellfun(@isempty,regexp(blk,'--> massdist__start_index=')));
idxrem = [];

% Find the headers before and after the text
for i = 1:numel(idx)
    idxrem = cat(1,idxrem,[idx(i)-1,idx(i)+1:idx(i)+4]');
end

% Get everything else other than the input data
idxrem      = cat(1,[1:idx(1)-2]',idxrem);
blk(idxrem) = [];

idx = find(~cellfun(@isempty,regexp(blk,'--> massdist__start_index=')));
r   = numel(idx);
idx = cat(1,idx,numel(blk) + 1);

for i = 1:r
    
    tab = blk(idx(i):idx(i + 1) - 1);
    % Extract name
    tabname = tab{1};
    
    idxtab          = findstr(tabname,'--> massdist__start_index=');
    idxtab          = [1:8,idxtab:length(tabname)];
    tabname(idxtab) = [];
    tabname         = strtrim(tabname);
    
    d(i).Name = tabname;
    tab(1) = [];
    
    % Extract data
    idxtab = find(~cellfun(@isempty,regexp(tab,'^C')));
    
    tab(idxtab)       = [];
    tab               = char(tab);
    idxexl            = findstr(tab(1,:),'!');
    tab(:,idxexl:end) = [];
    tab               = cellstr(tab);
    
    tab = regexprep(tab,'E\+(\w*)\+','E\+$1 \+');
    tab = regexprep(tab,'E\+(\w*)\-','E\+$1 \-');
    tab = regexprep(tab,'E\-(\w*)\+','E\-$1 \+');
    tab = regexprep(tab,'E\-(\w*)\-','E\-$1 \-');
    
    tab = char(tab);
    tab = str2num(tab);
    
    d(i).Pts.x    = tab(:,4);
    d(i).Pts.y    = tab(:,5);
    d(i).Pts.z    = tab(:,6);
    d(i).Pts.Mass = tab(:,7);
    
    % Calculate CG
    cgX = sum(d(i).Pts.x .* d(i).Pts.Mass) / sum(d(i).Pts.Mass);
    cgY = sum(d(i).Pts.y .* d(i).Pts.Mass) / sum(d(i).Pts.Mass);
    cgZ = sum(d(i).Pts.z .* d(i).Pts.Mass) / sum(d(i).Pts.Mass);
    
    d(i).Part.Mass = sum(d(i).Pts.Mass);
    d(i).Part.CG = [cgX,cgY,cgZ];
   
    
    % Calculate Inertia around CG
    dX = d(i).Pts.x - ones(numel(d(i).Pts.x),1) .* cgX;
    dY = d(i).Pts.y - ones(numel(d(i).Pts.y),1) .* cgY;
    dZ = d(i).Pts.z - ones(numel(d(i).Pts.z),1) .* cgZ;
    
    Ixx = sum((dY.^2 + dZ.^2) .* d(i).Pts.Mass);
    Iyy = sum((dX.^2 + dZ.^2) .* d(i).Pts.Mass);
    Izz = sum((dX.^2 + dY.^2) .* d(i).Pts.Mass);
    Ixy = sum(dX .* dY .* d(i).Pts.Mass);
    Ixz = sum(dX .* dZ .* d(i).Pts.Mass);
    Iyz = sum(dY .* dZ .* d(i).Pts.Mass);
    
    d(i).Part.Inertia = [ Ixx,Ixy,Ixz;
        Ixy,Iyy,Iyz;
        Ixz,Iyz,Izz];
end

% Remove any non-zero masses
dPart = [d.Part];
dMass = [dPart.Mass];
d(dMass == 0) = [];

end

