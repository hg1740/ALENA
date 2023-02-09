%% SYMWINGAERO Reflect wing aero definition around XZ-plane
function obj = symWingAero(obj)

Offsetidx = obj.Opts.Struct.AeroOffsetID;
StructOffsetidx = obj.Opts.Struct.StrucOffsetID;

% Reflect the caero cards
nCaero = numel(obj.Mdl.Caero);

for i = 1:nCaero
    if ~(strcmpi(obj.Mdl.Caero(i).part,'vtp') || strcmpi(obj.Mdl.Caero(i).part,'fuselage'))
        CaeroObj = Caero;
        CaeroObj.id          =  obj.Mdl.Caero(i).id + Offsetidx;
        CaeroObj.startX      =  obj.Mdl.Caero(i).startX;
        CaeroObj.startY      = -obj.Mdl.Caero(i).startY;
        CaeroObj.startZ      =  obj.Mdl.Caero(i).startZ;
        CaeroObj.twist1      =  obj.Mdl.Caero(i).twist1;
        CaeroObj.twist2      =  obj.Mdl.Caero(i).twist2;
        CaeroObj.c           =  obj.Mdl.Caero(i).c;
        CaeroObj.t           =  obj.Mdl.Caero(i).t;
        CaeroObj.b           = -obj.Mdl.Caero(i).b;
        CaeroObj.sw          = -obj.Mdl.Caero(i).sw;
        CaeroObj.rootAirfoil =  obj.Mdl.Caero(i).rootAirfoil;
        CaeroObj.tipAirfoil  =  obj.Mdl.Caero(i).tipAirfoil;
        CaeroObj.dih         =  -obj.Mdl.Caero(i).dih; % TODO changed
        CaeroObj.meshType    =  obj.Mdl.Caero(i).meshType;
        CaeroObj.cid         =  obj.Mdl.Caero(i).cid;
        CaeroObj.ny          =  obj.Mdl.Caero(i).ny;
        CaeroObj.nx          =  obj.Mdl.Caero(i).nx;
        CaeroObj.csType      =  obj.Mdl.Caero(i).csType;
        CaeroObj.part        =  obj.Mdl.Caero(i).part;
        CaeroObj.fc          =  obj.Mdl.Caero(i).fc;
        CaeroObj.fnx         =  obj.Mdl.Caero(i).fnx;
        CaeroObj.symetric    =  obj.Mdl.Caero(i).symetric;
        CaeroObj.flapped     =  obj.Mdl.Caero(i).flapped;
        
        csId = obj.Mdl.Caero(i).csId;
        
        if ~isempty(csId)
            CaeroObj.csId = sprintf('P%s',csId);
        else
            CaeroObj.csId = csId;
        end
        
        CaeroObj.csData = obj.Mdl.Caero(i).csData;
        CaeroObj.csLink = -1;
        obj.Mdl.Caero = cat(2,obj.Mdl.Caero,CaeroObj);
    end
end

% Isolate the Wing grids:
CaeroObj = findobj(obj.Mdl.Caero,'part','StbdWing');

% Isolate the LHS Wing
CaeroObjEndY = [CaeroObj.startY]' + [CaeroObj.b]';

CaeroObjId  = [CaeroObj.id]';
LWingId     = CaeroObjId(CaeroObjEndY < 0);

% Add the ids to the part object
PartObj      = PartId;
PartObj.id   = 1;
PartObj.part = 'PortWing';
PartObj.type = 'CAERO';
PartObj.data = LWingId';

% Add to the objects:
obj.Mdl.PartId  = cat(2,obj.Mdl.PartId,PartObj);

% Isolate the Wing grids:
CaeroObj = findobj(obj.Mdl.Caero,'part','StbdHTP');

% Isolate the LHS Wing
CaeroObjEndY = [CaeroObj.startY]' + [CaeroObj.b]';

CaeroObjId  = [CaeroObj.id]';
LHTPId     = CaeroObjId(CaeroObjEndY < 0);

% Add the ids to the part object
PartObj      = PartId;
PartObj.id   = 1;
PartObj.part = 'PortHTP';
PartObj.type = 'CAERO';
PartObj.data = LHTPId';

% Add to the objects:
obj.Mdl.PartId  = cat(2,obj.Mdl.PartId,PartObj);

% Reflect the set1 cards
gridPts = struct2mat(obj.Mdl.Grid,'id');
gridCoordY = struct2mat(obj.Mdl.Grid,'coord',2);

nSets = numel(obj.Mdl.Sets);

for i = 1:nSets
    
    if ~(strcmpi(obj.Mdl.Sets(i).part,'vtp') || strcmpi(obj.Mdl.Sets(i).part,'fuselage'))
        
        SetsObj = Sets;
        
        SetsObj.id   = obj.Mdl.Sets(i).id + Offsetidx;
        
        % Cannot create duplicates at the centre line
        data = struct2mat(obj.Mdl.Sets(i),'data');
        [~,~,idx] = intersect(data,gridPts);
        gridy  = gridPts(idx);
        idx = find(abs(gridCoordY(idx)) < 1e-6);
        
        ptsLabels = gridy(idx);
        [~,idx,~] = intersect(data,ptsLabels);
        
        SetsObj.data = obj.Mdl.Sets(i).data + StructOffsetidx;
        SetsObj.data(idx) = data(idx);
        SetsObj.part = obj.Mdl.Sets(i).part;
        obj.Mdl.Sets = cat(2,obj.Mdl.Sets,SetsObj);
    end
end

% Reflect the spline cards
nSpline = numel(obj.Mdl.Spline);

for i = 1:nSpline
    if ~(strcmpi(obj.Mdl.Spline(i).part,'vtp') || strcmpi(obj.Mdl.Spline(i).part,'fuselage'))
        
        SplineObj = Spline; % TODO What is this?
        
        SplineObj.id   = obj.Mdl.Spline(i).id + Offsetidx;
        SplineObj.aid  = obj.Mdl.Spline(i).aid + Offsetidx;
        SplineObj.p1   = obj.Mdl.Spline(i).p1;
        SplineObj.p2   = obj.Mdl.Spline(i).p2;
        SplineObj.w    = obj.Mdl.Spline(i).w;
        SplineObj.rmx  = obj.Mdl.Spline(i).rmx;
        SplineObj.cond = obj.Mdl.Spline(i).cond;
        SplineObj.set  = obj.Mdl.Spline(i).set + Offsetidx;
        SplineObj.part = obj.Mdl.Spline(i).part;
        obj.Mdl.Spline = cat(2,obj.Mdl.Spline,SplineObj);
    end
end



end

