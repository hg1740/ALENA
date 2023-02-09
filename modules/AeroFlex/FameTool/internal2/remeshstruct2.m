function obj = remeshstruct2(obj)

% Determine where the straight sections are. This will be used to
% distribute the beam elements
eta = unique(obj.Fame.Geometry.Wing.Outline(:,2))/max(unique(obj.Fame.Geometry.Wing.Outline(:,2)));

nelem = 32;

% Build vectors of all the grid and element data
x        = struct2mat(obj.Mdl.Grid,'coord',1);
y        = struct2mat(obj.Mdl.Grid,'coord',2);
z        = struct2mat(obj.Mdl.Grid,'coord',3);
gridId   = [obj.Mdl.Grid.id]';
cbarConn = cell2mat({obj.Mdl.Cbar.conn}');
pbarId   = cell2mat({obj.Mdl.Cbar.pid}');
pbarIm   = cell2mat({obj.Mdl.Pbar.i}');
pbarJ    = cell2mat({obj.Mdl.Pbar.j}');
pbarA    = cell2mat({obj.Mdl.Pbar.a}');
cbarOrient    = cell2mat({obj.Mdl.Cbar.orient}');

halfSpan = max(y);

etaOld = y / halfSpan;
etaOldMid = 0.5*(etaOld(1:end-1) + etaOld(2:end));

% Work out the number of beam elements to put in each section:
secFrac   = diff(eta);
nBeams    = secFrac * nelem;
nBeams    = round(nBeams);
SumBeams  = cumsum(nBeams);
ip        = [0;SumBeams(1:end - 1)];

% Determine new eta values
etaNew = [];
tol    = 0;

for i = 1:numel(nBeams)
    ds     = (eta(i + 1) - eta(i)) / nBeams(i);
    etaNew = cat(2,etaNew,eta(i) + tol:ds:eta(i + 1) - tol);
end

% Rebuild grid elements
etaNew      = unique(etaNew);
etaNewMid   = 0.5*(etaNew(1:end-1) + etaNew(2:end));
startNumber = obj.Opts.Struct.WingStartID;

xNew = interp1(etaOld,x,etaNew,'linear')';
yNew = interp1(etaOld,y,etaNew,'linear')';
zNew = interp1(etaOld,z,etaNew,'linear')';

for i = 1:numel(xNew)
    GridObj(i)       = Grid;
    GridObj(i).id    = i + startNumber;
    GridObj(i).cs    = 0;
    GridObj(i).coord = [xNew(i), yNew(i), zNew(i)];
    GridObj(i).cd    = 0;
    GridObj(i).ps    = 0;
    GridObj(i).seid  = 0;
    GridObj(i).part  = 'StbdWing';
    GridObj(i).type  = 'Beam';
end

% Add the ids to the part object
PartObj      = PartId;
PartObj(1).id   = 1;
PartObj(1).part = 'StbdWing';
PartObj(1).type = 'GRID';
PartObj(1).data = [GridObj.id];

% Rebuild the beeam elements
pbarJNew        = interp1(etaOldMid,pbarJ,etaNewMid,'linear')';
pbarImNew = [];
pbarImNew(:,1)  = interp1(etaOldMid,pbarIm(:,1),etaNewMid,'linear')';
pbarImNew(:,2)  = interp1(etaOldMid,pbarIm(:,2),etaNewMid,'linear')';
pbarImNew(:,3)  = interp1(etaOldMid,pbarIm(:,3),etaNewMid,'linear')';
pbarANew        = interp1(etaOldMid,pbarA,etaNewMid,'linear')';
cbarOrientNew = [];
cbarOrientNew(:,1) = interp1(etaOldMid,cbarOrient(:,1),etaNewMid,'linear')';
cbarOrientNew(:,2) = interp1(etaOldMid,cbarOrient(:,2),etaNewMid,'linear')';
cbarOrientNew(:,3) = interp1(etaOldMid,cbarOrient(:,3),etaNewMid,'linear')';

for i = 1:numel(pbarANew)
   
    PbarObj(i)     = Pbar;
    PbarObj(i).id  = i + startNumber;
    PbarObj(i).mat = 1;
    PbarObj(i).a   = pbarANew(i);
    PbarObj(i).i   = pbarImNew(i,:);
    PbarObj(i).j   = pbarJNew(i);
    PbarObj(i).part = 'StbdWing';
    
    CbarObj(i)        = Cbar;
    CbarObj(i).id     = i + startNumber;
    CbarObj(i).pid    = i + startNumber;
    CbarObj(i).conn   = [GridObj(i).id,GridObj(i + 1).id];
    CbarObj(i).orient = cbarOrientNew(i,:);
    CbarObj(i).offset = 'GGG';
    CbarObj(i).part = 'StbdWing';
    
end

% Add the ids to the part object
PartObj(2).id   = 1;
PartObj(2).part = 'StbdWing';
PartObj(2).type = 'CBAR';
PartObj(2).data = [CbarObj.id];

% Rebuild spc elements
Xn = struct2mat(GridObj,'coord',1);
Yn = struct2mat(GridObj,'coord',2);
Zn = struct2mat(GridObj,'coord',3);

for i = 1:numel(obj.Mdl.Spc)
    
    idxSpc = find(obj.Mdl.Spc(i).grids == gridId);
    
    Xo = repmat(x(idxSpc(1)),size(Xn));
    Yo = repmat(y(idxSpc(1)),size(Yn));
    Zo = repmat(z(idxSpc(1)),size(Zn));
    
    distance = sqrt((Xo - Xn).^2 + (Yo - Yn).^2 + (Zo - Zn).^2);
    
    [~,idxMinDist] = min(distance);
    
    obj.Mdl.Spc(i).grids = GridObj(idxMinDist).id;
end


% Assign to model object
obj.Mdl.Grid    = GridObj;
obj.Mdl.Cbar    = CbarObj;
obj.Mdl.Pbar    = PbarObj;
obj.Mdl.PartId  = cat(2,obj.Mdl.PartId,PartObj);

obj.Fame.BeamModel.Wing.UnsampledData.Eta     = etaOld';
obj.Fame.BeamModel.Wing.UnsampledData.Grids   = [x,y,z];
obj.Fame.BeamModel.Wing.UnsampledData.Beam_I11 = pbarIm(:,1);
obj.Fame.BeamModel.Wing.UnsampledData.Beam_I22 = pbarIm(:,2);
obj.Fame.BeamModel.Wing.UnsampledData.Beam_I12 = pbarIm(:,3);
obj.Fame.BeamModel.Wing.UnsampledData.Beam_J   = pbarJ;
obj.Fame.BeamModel.Wing.UnsampledData.Beam_A   = pbarA;
obj.Fame.BeamModel.Wing.UnsampledData.Beam_Orient   = cbarOrient;

obj.Fame.BeamModel.Wing.Eta      = etaNew';
obj.Fame.BeamModel.Wing.Grids    = [xNew,yNew,zNew];
obj.Fame.BeamModel.Wing.Beam_I11 = pbarImNew(:,1);
obj.Fame.BeamModel.Wing.Beam_I22 = pbarImNew(:,2);
obj.Fame.BeamModel.Wing.Beam_I12 = pbarImNew(:,3);
obj.Fame.BeamModel.Wing.Beam_J   = pbarJNew;
obj.Fame.BeamModel.Wing.Beam_A   = pbarANew;
obj.Fame.BeamModel.Wing.Beam_Orient     = cbarOrientNew;
obj.Fame.BeamModel.Wing.YoungsMod  = obj.Mdl.Mat(1).e;
obj.Fame.BeamModel.Wing.ShearMod   = obj.Mdl.Mat(1).g;
obj.Fame.BeamModel.Wing.Poissions  = obj.Mdl.Mat(1).nu;
obj.Fame.BeamModel.Wing.Density    = obj.Mdl.Mat(1).rho;
end

